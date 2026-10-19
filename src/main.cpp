#include "algorithm/algorithm.h"
#include "algorithm/optimization.h"
#include "algorithm/search_helpers.h"
#include "algorithm/search_types.h"
#include "cellTree.h"
#include "data/cell.h"
#include "data/facet.h"
#include "data/halfspace.h"
#include "data/slab.h"
#include "data/user.h"
#include "geometry/qhull.h"
#include "header.h"
#include "spatial/rtree.h"
#include "utils/utils.h"
#include "utils/anytime.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <climits>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <numeric>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>

int main(const int argc, const char **argv) {
  using std::chrono::steady_clock;
  const auto programStart = steady_clock::now();
  auto finish = [&]() {
    const double totalSeconds =
        std::chrono::duration<double>(steady_clock::now() - programStart)
            .count();
    cout << "TIME_SEC: " << std::fixed << std::setprecision(6)
         << totalSeconds << endl;
    return 0;
  };
  auto invalidWca = [&](const std::string& message) {
    cerr << "[INVALID-WCA] " << message << endl;
    finish();
    return 2;
  };

  bool anytimeDisabled = false;
  if (const char* envValue = std::getenv("ANYTIME_DISABLE")) {
    const string flag(envValue);
    anytimeDisabled =
        !flag.empty() && flag != "0" && flag != "false" && flag != "FALSE";
  }

  double anytimeIntervalSec = 0.0;
  if (!anytimeDisabled) {
    if (const char* envValue = std::getenv("ANYTIME_INTERVAL_SEC")) {
      if (envValue[0] != '\0') {
        char* end = nullptr;
        const double parsed = std::strtod(envValue, &end);
        if (end != envValue && parsed > 0.0 && std::isfinite(parsed)) {
          anytimeIntervalSec = parsed;
        } else {
          cerr << "[ANYTIME] Ignoring invalid ANYTIME_INTERVAL_SEC="
               << envValue << endl;
        }
      }
    }
  }
  double nextAnytimeIntervalSec = anytimeIntervalSec;

  double anytimeStopGapFraction = -1.0;
  auto parseAnytimeGapTarget = [&](const char* raw,
                                   const char* sourceName) {
    if (raw == nullptr || raw[0] == '\0') {
      return;
    }
    char* end = nullptr;
    const double parsed = std::strtod(raw, &end);
    if (end == raw || parsed < 0.0 || !std::isfinite(parsed)) {
      cerr << "[ANYTIME] Ignoring invalid " << sourceName << "="
           << raw << endl;
      return;
    }

    // Accept either a fraction (-agap 0.30) or a percent (-agap 30).
    anytimeStopGapFraction = (parsed > 1.0) ? parsed / 100.0 : parsed;
  };

  if (!anytimeDisabled) {
    parseAnytimeGapTarget(std::getenv("ANYTIME_STOP_GAP_PERCENT"),
                          "ANYTIME_STOP_GAP_PERCENT");
    parseAnytimeGapTarget(Param::read(argc, argv, "-agap", ""), "-agap");
  }

  bool anytimeStopRequested = false;
  // Search incumbents and LP bounds share the reduced slab-index scale.
  // Their comparisons are invariant under adding this common constant.
  long long rankOffset = 0;

  auto requestAnytimeStopIfReady =
      [&](const string& event,
          int incumbentArf,
          const SubRegion* current,
          const SubRegionQueue* innerQueue,
          const SubRegionQueue* globalQueue) {
        if (anytimeDisabled || anytimeStopGapFraction < 0.0 ||
            anytimeStopRequested || incumbentArf == INT_MAX) {
          return false;
        }

        AnytimeGapStatus status;
        if (!anytimeGapTargetReached(incumbentArf, anytimeStopGapFraction,
                                     current, innerQueue, globalQueue,
                                     &status, rankOffset)) {
          return false;
        }

        logAnytimeSnapshot(event, incumbentArf, programStart, current,
                           innerQueue, globalQueue, rankOffset);

        std::ios::fmtflags oldFlags = cout.flags();
        std::streamsize oldPrecision = cout.precision();
        cout << std::fixed << std::setprecision(3)
             << "[ANYTIME-STOP] target_gap_percent="
             << anytimeStopGapFraction * 100.0
             << " achieved_gap_percent=" << status.relativeGap * 100.0
             << " incumbent=" << (long long)incumbentArf + rankOffset
             << " int_bound=" << status.integerLowerBound
             << " certified_gap=" << status.certifiedGap
             << " event=" << event << endl;
        cout.flags(oldFlags);
        cout.precision(oldPrecision);

        anytimeStopRequested = true;
        return true;
      };

  auto maybeLogAnytimeInterval =
      [&](int incumbentArf,
          const SubRegion* current,
          const SubRegionQueue* innerQueue,
          const SubRegionQueue* globalQueue) {
        if (anytimeDisabled || anytimeIntervalSec <= 0.0 ||
            incumbentArf == INT_MAX) {
          return;
        }
        const double elapsed =
            std::chrono::duration<double>(steady_clock::now() - programStart)
                .count();
        if (elapsed + 1e-9 < nextAnytimeIntervalSec) {
          return;
        }
        logAnytimeSnapshot("interval", incumbentArf, programStart, current,
                           innerQueue, globalQueue, rankOffset);
        requestAnytimeStopIfReady("stop_interval", incumbentArf, current,
                                  innerQueue, globalQueue);
        do {
          nextAnytimeIntervalSec += anytimeIntervalSec;
        } while (elapsed + 1e-9 >= nextAnytimeIntervalSec);
      };

  int anytimeBestArf = INT_MAX;
  auto logNewIncumbent =
      [&](const string& event,
          int incumbentArf,
          const SubRegion* current,
          const SubRegionQueue* innerQueue,
          const SubRegionQueue* globalQueue) {
        if (anytimeDisabled) {
          return;
        }
        if (incumbentArf >= anytimeBestArf) {
          return;
        }
        anytimeBestArf = incumbentArf;
        logAnytimeSnapshot(event, incumbentArf, programStart, current,
                           innerQueue, globalQueue, rankOffset);
        requestAnytimeStopIfReady("stop_" + event, incumbentArf, current,
                                  innerQueue, globalQueue);
      };
  // Parse the active CLI used by the current solver. "-N" is kept only so old
  // scripts still run, but it no longer changes the algorithm.
  int dim = atoi(Param::read(argc, argv, "-d", "")); // number of Dimensions
  (void)Param::read(argc, argv, "-N", ""); // Legacy CLI arg, ignored
  const char *productFile = Param::read(argc, argv, "-pf", "");
  const char *userFile = Param::read(argc, argv, "-uf", "");
  const char *verifyComboArg = Param::read(argc, argv, "-vc", "");
  const char *wcaMode = Param::read(argc, argv, "-wca", "fixed");
  const double wcaRadius = std::atof(Param::read(argc, argv, "-wr", "0.4"));
  int wcaHyperplanes = std::atoi(Param::read(argc, argv, "-wn", "0"));
  const uint64_t wcaSeed =
      std::strtoull(Param::read(argc, argv, "-ws", "1"), nullptr, 10);
  vector<int> verifyCombination;
  if (verifyComboArg != nullptr && verifyComboArg[0] != '\0') {
    std::stringstream ss(verifyComboArg);
    std::string token;
    while (std::getline(ss, token, ',')) {
      if (!token.empty()) {
        verifyCombination.push_back(std::stoi(token));
      }
    }
  }

  list<User> usersList;
  vector<User *> users;
  vector<Halfspace> inequalities;
  vector<vector<double>> intersections, TrimIntersections;
  vector<double> solution;
  vector<Facet> facetList;

  // Load users and products, then build the weight-space feasibility region.
  // The default fixed WCA preserves the historical deterministic setup.
  // Benchmark runs can instead request a generated WCA with deterministic
  // hyperplane directions via -wca generated -wr <radius> -wn <N> -ws <seed>.
  // get_ubc_vrep(..., FULL_UNIT_BOX) then adds 0 <= w_i <= 1 and returns the
  // WCA polytope in vertex/facet form for the later search.
  loadUsers(users, usersList, userFile, dim);
  const vector<User*> originalUsers = users;
  Rtree *rtree = loadProducts(productFile, dim);
  if (std::string(wcaMode) == "generated") {
    if (wcaHyperplanes <= 0) {
      wcaHyperplanes = dim;
    }
    CreateWCA(inequalities, dim, wcaHyperplanes, wcaRadius, wcaSeed);
  } else if (std::string(wcaMode) == "fixed") {
    appendFixedWCA(inequalities, dim);
  } else {
    throw std::invalid_argument("Unsupported -wca mode: " +
                                std::string(wcaMode));
  }
  WCAResult ubc = get_ubc_vrep(dim, inequalities, BoxConstraints::FULL_UNIT_BOX);

  if (ubc.vertices.empty()) {
    // Qhull failed to compute WCA vertices (degenerate polytope or infeasible).
    // This is not a valid successful solve because no ARF can be reported.
    return invalidWca("WCA vertex computation returned no vertices");
  }

  intersections = ubc.vertices;

  for (size_t i = 0; i < ubc.facets.size(); ++i) {
    FacetVRep &f = ubc.facets[i];
    facetList.push_back(Facet(f, ubc.vertices));
  }

  vector<User *> tmphullUsers;
  vector<User *> nonProductUsers;
  vector<User *> hullUsers;
  vector<User *> toErase;
  bool start = true;

  // Preprocess the users before the exact search starts:
  // 1) find hull users,
  // 2) trim the WCA edges against those hull users,
  // 3) use the trimmed boundary to prune products/users that can never matter.
  while (nonProductUsers.size() > 0 || start) {
    tmphullUsers.clear();
    nonProductUsers.clear();
    TrimIntersections.clear();

    if (start) {
      HullVRep HconexHull = convex_hull_vrep_from_users(users);
      hull_vertices_users(HconexHull, users, tmphullUsers);
      start = false;
    } else {
      if (hullUsers.size() > dim + 1) {
        HullVRep HconexHull = convex_hull_vrep_from_users(hullUsers);
        hull_vertices_users(HconexHull, hullUsers, tmphullUsers);
      }
    }

    TrimWcaEdges(intersections, tmphullUsers, TrimIntersections);

    if (TrimIntersections.size() == 1) {
      const vector<double>& onlyPoint = TrimIntersections.front();
      const long long exactArf =
          computeExactPointArfFromProductFile(productFile, originalUsers, onlyPoint, dim);
      cout << "[DEGENERATE-WCA] Trimmed WCA has one point; exact ARF computed directly." << endl;
      cout << "[RESULT] Feasible point: ";
      for (double val : onlyPoint) {
        cout << std::setprecision(17) << val << " ";
      }
      cout << endl;
      cout << "\nARF: " << exactArf << endl;
      return finish();
    }

    if (!hullUsers.empty()) {
      eraseUsersInPlace(tmphullUsers, hullUsers);
    }

    auto computePreprocessScoreWindow =
        [&](User* user, double& minscore, double& maxscore) {
          minscore = std::numeric_limits<double>::infinity();
          maxscore = -std::numeric_limits<double>::infinity();

          for (const vector<double>& point : TrimIntersections) {
            const double score = user->CalculateScore(point);
            maxscore = std::max(maxscore, score);
            minscore = std::min(minscore, score);
          }

          if (!std::isfinite(maxscore)) {
            maxscore = minscore;
          }
        };

    for (int i = 0; i < hullUsers.size(); i++) {
      double maxscore;
      double minscore;
      computePreprocessScoreWindow(hullUsers[i], minscore, maxscore);
      hullUsers[i]->checkProduct(maxscore, minscore);
      if (hullUsers[i]->getProductsSize() == 0) {
        nonProductUsers.push_back(hullUsers[i]);
        toErase.push_back(hullUsers[i]);
      } else {
        hullUsers[i]->setEdgeScores(minscore, maxscore);
      }
    }

    if (nonProductUsers.size() > 0) {
      eraseUsersInPlace(hullUsers, nonProductUsers);
    }

    for (int i = 0; i < tmphullUsers.size(); i++) {
      double maxscore;
      double minscore;
      computePreprocessScoreWindow(tmphullUsers[i], minscore, maxscore);
      rtreePruning(tmphullUsers[i], rtree->m_root, dim, minscore, maxscore);
      if (tmphullUsers[i]->getProductsSize() == 0) {
        nonProductUsers.push_back(tmphullUsers[i]);
        toErase.push_back(tmphullUsers[i]);
      } else {
        tmphullUsers[i]->setEdgeScores(minscore, maxscore);
        hullUsers.push_back(tmphullUsers[i]);
      }
    }
  }

  // Remove users that became irrelevant during the iterative pruning stage.
  eraseUsersInPlace(users, toErase);

  // Hull users define the boundary geometry used later by rank lower bounds.
  for (int i = 0; i < hullUsers.size(); i++) {
    hullUsers[i]->createSlabs();
  }

  vector<User *> nonSlabsUsers;
  for (int i = 0; i < users.size(); i++) {
    double maxscore = -std::numeric_limits<double>::infinity();
    double minscore = std::numeric_limits<double>::infinity();
    for (const vector<double>& point : TrimIntersections) {
      const double score = users[i]->CalculateScore(point);
      maxscore = std::max(maxscore, score);
      minscore = std::min(minscore, score);
    }
    if (!std::isfinite(maxscore)) {
      maxscore = minscore;
    }
    if (find(hullUsers.begin(), hullUsers.end(), users[i]) != hullUsers.end())
      continue;
    rtreePruning(users[i], rtree->m_root, dim, minscore, maxscore);
    users[i]->setEdgeScores(minscore, maxscore);
    users[i]->createSlabs();
    if (users[i]->checkSlabEmpty()) {
      nonSlabsUsers.push_back(users[i]);
    }
  }

  eraseUsersInPlace(users, nonSlabsUsers);
  for (const User& user : usersList) rankOffset += user.getRankOffset();
  delete rtree;
  rtree = nullptr;

  // Verification mode: bypass the search and only test the requested slab
  // combination against the current preprocessed user set.
  if (!verifyCombination.empty()) {
    if ((int)verifyCombination.size() != (int)users.size()) {
      throw std::runtime_error(
          "Verify combination size does not match number of users after preprocessing");
    }

    vector<pair<int, User*>> combination;
    int verifyArf = 0;
    for (int u = 0; u < (int)users.size(); ++u) {
      combination.push_back({verifyCombination[u], users[u]});
      verifyArf += verifyCombination[u];
    }

    vector<double> verifySolution;
    bool feasible =
        checkFeasible(combination, users, hullUsers, inequalities, dim,
                      verifySolution, false);
    cout << "[VERIFY-CELL] Feasible: " << (feasible ? "yes" : "no") << endl;
    cout << "[VERIFY-CELL] ARF: " << (long long)verifyArf + rankOffset << endl;
    if (feasible) {
      cout << "[VERIFY-CELL] Witness point: ";
      for (double val : verifySolution) {
        cout << val << " ";
      }
      cout << endl;
    }

    vector<double> exactVerifySolution;
    bool exactFeasible =
        checkFeasibleExact(combination, users, inequalities, dim,
                           exactVerifySolution);
    cout << "[VERIFY-CELL-EXACT] Feasible: "
         << (exactFeasible ? "yes" : "no") << endl;
    if (exactFeasible) {
      cout << "[VERIFY-CELL-EXACT] Witness point: ";
      for (double val : exactVerifySolution) {
        cout << val << " ";
      }
      cout << endl;
    }
    return finish();
  }

  if (TrimIntersections.empty()) {
    return invalidWca("Trimmed WCA has no intersections");
  }
  // Start from the trimmed WCA point with the best exact rank sum. If multiple
  // points tie on ARF, prefer the one with larger raw total score.
  int bestStartIdx = 0;
  int bestStartARF = INT_MAX;
  double bestStartScore = -INFINITY;
  vector<int> startSlabs;
  for (int j = 0; j < (int)TrimIntersections.size(); ++j) {
    const int candidateARF =
        computePointArf(users, TrimIntersections[j], startSlabs);

    double candidateScore = 0.0;
    for (User* user : users) {
      candidateScore += user->CalculateScore(TrimIntersections[j]);
    }

    if (candidateARF < bestStartARF ||
        (candidateARF == bestStartARF && candidateScore > bestStartScore)) {
      bestStartARF = candidateARF;
      bestStartScore = candidateScore;
      bestStartIdx = j;
    }
  }
  vector<double> StartCorner = TrimIntersections[bestStartIdx];

  // Convert the starting point into a concrete slab combination and, at the
  // same time, preload each facet's geometric lower bound and initial
  // per-user border ranges. Phase 1 relies on these facet summaries before any
  // LP tightening happens.
  Cell StartingCell;
  Slab userSlab;
  for (int i = 0; i < users.size(); i++) {
    double score = users[i]->CalculateScore(StartCorner);
    userSlab = users[i]->findClosestSlab(score);
    StartingCell.getCombination().push_back({userSlab, users[i]});

    for (int j = 0; j < facetList.size(); j++) {
      facetList[j].calculateRankFromUser(users[i]);
    }
  }

  StartingCell.setCombinationIndex();

  std::sort(facetList.begin(), facetList.end(),
            [](const Facet& a, const Facet& b) {
              return a.getLowerBound() < b.getLowerBound();
            });

  vector<Cell> infeasibleCells;
  const double EPS = 1e-6; // Tolerance for floating point precision in LP bounds
  solution = StartCorner;
  Cell bestCell = StartingCell;
  int minARF = bestCell.getARF();
  logNewIncumbent("initial_corner", minARF, nullptr, nullptr, nullptr);

  // The canonical search starts from the preprocessed boundary corner.
  maybeLogAnytimeInterval(minARF, nullptr, nullptr, nullptr);
  int numUsers = users.size();
  uint64_t borderLpSolveCount = 0;
  double borderLpSolveSeconds = 0.0;
  uint64_t relaxedLpSolveCount = 0;
  double relaxedLpSolveSeconds = 0.0;

  SubRegionQueue globalQueue;

  // Phase 1: create one root sub-region per facet from the facet's geometric
  // border seed, then keep only facets whose stored selected lower bound can
  // still beat the initial incumbent.
  for (int fIdx = 0;
       fIdx < (int)facetList.size() && !anytimeStopRequested;
       ++fIdx) {
    Facet& facet = facetList[fIdx];

    if (facet.getLowerBound() >= minARF) {
      continue;
    }

    vector<int> minBorders = facet.getMinBorders();
    vector<int> maxBorders = facet.getMaxBorders();

    SubRegion seed;
    seed.facetIdx = fIdx;
    seed.minBorders = minBorders;
    seed.maxBorders = maxBorders;
    seed.splitHalfspaces = {};  // Root of this facet: no splits yet

    solveRelaxedSubRegion(seed, facetList, users, dim,
                          relaxedLpSolveCount, relaxedLpSolveSeconds);

    if (seed.relaxedARF >= minARF) {
      continue;
    }

    globalQueue.push(seed);
  }

  // A global certificate is valid only after every relevant facet has either
  // been pruned or inserted into the root queue. During root construction,
  // an infeasible seed may carry a large sentinel bound and unvisited facets
  // have not contributed bounds yet.
  maybeLogAnytimeInterval(minARF, nullptr, nullptr, &globalQueue);
  if (!anytimeStopRequested) {
    requestAnytimeStopIfReady("stop_roots_ready", minARF, nullptr, nullptr,
                              &globalQueue);
  }

  // Phase 2: exact branch-and-bound inside the surviving facets.
  while (!globalQueue.empty() && !anytimeStopRequested) {
    SubRegion topFacet = globalQueue.top();
    globalQueue.pop();
    maybeLogAnytimeInterval(minARF, &topFacet, nullptr, &globalQueue);
    if (requestAnytimeStopIfReady("stop_global_pop", minARF, &topFacet,
                                  nullptr, &globalQueue)) {
      break;
    }

    if (topFacet.relaxedARF >= minARF) {
      continue;
    }

    Facet& facet = facetList[topFacet.facetIdx];

    SubRegionQueue innerQueue;
    vector<UserHullCache> facetHullCache(numUsers);
    innerQueue.push(topFacet);

    while (!innerQueue.empty() && !anytimeStopRequested) {
      SubRegion current = innerQueue.top();
      innerQueue.pop();
      maybeLogAnytimeInterval(minARF, &current, &innerQueue, &globalQueue);
      if (requestAnytimeStopIfReady("stop_inner_pop", minARF, &current,
                                    &innerQueue, &globalQueue)) {
        break;
      }

      if (current.relaxedARF >= minARF) {
        continue;
      }

      double relaxedMinARF = current.relaxedARF;
      vector<double>& relaxedRanks = current.relaxedRanks;
      vector<double>& wLP = current.wLP;

      if ((int)wLP.size() < dim || (int)relaxedRanks.size() < numUsers) {
        continue;
      }

      // The LP witness gives a real point in w-space. Re-evaluate it against
      // the true slab structure to see whether the relaxed solution itself is
      // already a better exact incumbent.
      int realARF = 0;
      vector<int> integerRanks(numUsers);
      for (int u = 0; u < numUsers; ++u) {
        double score = 0.0;
        const vector<double>& coords = users[u]->getCoordinates();
        for (int d = 0; d < dim; ++d) score += coords[d] * wLP[d];
        Slab& closestSlab = users[u]->findClosestSlab(score);
        integerRanks[u] = closestSlab.getIndex();
        realARF += integerRanks[u];
      }

      if (realARF < minARF) {
        minARF = realARF;
        vector<pair<Slab, User*>> cellCombo;
        for (int u = 0; u < numUsers; ++u)
          cellCombo.push_back({users[u]->getSlab(integerRanks[u]), users[u]});
        bestCell = Cell(cellCombo);
        solution = wLP;
        logNewIncumbent("incumbent_lp_witness", minARF, &current,
                        &innerQueue, &globalQueue);
        if (anytimeStopRequested) {
          break;
        }
      }

      bool hasSplitCandidate = false;
      UserGap bestGap{ -1, 0.0, -1 };
      for (int u = 0; u < numUsers; ++u) {
        int ceilRelaxed = (int)ceil(current.relaxedRanks[u] - EPS);
        int gapInt = integerRanks[u] - ceilRelaxed;
        if (gapInt >= 1) {
          // Split on the slab boundary suggested by the relaxed rank, but
          // canonicalize it so tied thresholds map to a valid user split.
          int splitSlab = (int)floor(current.relaxedRanks[u] + EPS);
          int canonicalSplit =
              canonicalizeSplitSlab(users[u], splitSlab,
                                    current.minBorders[u], current.maxBorders[u]);
          if (canonicalSplit >= current.minBorders[u] &&
              canonicalSplit < current.maxBorders[u]) {
            const double gapValue = (double)gapInt;
            if (!hasSplitCandidate || gapValue > bestGap.gap) {
              bestGap = {u, gapValue, canonicalSplit};
              hasSplitCandidate = true;
            }
          }
        }
      }

      int totalUncertainty = 0;
      for (int u = 0; u < numUsers; ++u) {
          int range = current.maxBorders[u] - current.minBorders[u];
          if (range > 0) {
              totalUncertainty += range; 
          }
      }

      const int MILP_GAP_THRESHOLD = 1;
      const double MILP_AVG_SCORE_WIDTH_THRESHOLD = 5.0;
      double avgScoreWidth =
          averageScoreWidthInSlabUnits(current, users,
                                       facetList[current.facetIdx], dim);
      bool scoreWidthFitsMilp =
          avgScoreWidth <= MILP_AVG_SCORE_WIDTH_THRESHOLD;

      if (totalUncertainty == 0) {
        int fixedARF = 0;
        vector<pair<Slab, User*>> fixedCellCombo;
        vector<pair<int, User*>> fixedCombinationIndex;
        fixedCellCombo.reserve(numUsers);
        fixedCombinationIndex.reserve(numUsers);

        for (int u = 0; u < numUsers; ++u) {
          const int slabIdx = current.minBorders[u];
          fixedARF += slabIdx;
          fixedCellCombo.push_back({users[u]->getSlab(slabIdx), users[u]});
          fixedCombinationIndex.push_back({slabIdx, users[u]});
        }

        if (fixedARF < minARF) {
          vector<double> fixedSolution;
          if (checkFeasibleInFacet(fixedCombinationIndex, users, hullUsers,
                                   facetList[current.facetIdx], dim,
                                   fixedSolution, false)) {
            minARF = fixedARF;
            bestCell = Cell(fixedCellCombo);
            solution = fixedSolution;
            logNewIncumbent("incumbent_fixed_cell", minARF, &current,
                            &innerQueue, &globalQueue);
            if (anytimeStopRequested) {
              break;
            }
          }
        }
        if (anytimeStopRequested) {
          break;
        }
        continue;
      }

      if (!hasSplitCandidate) {
        int lpIntegralLowerBound = (int)ceil(relaxedMinARF);
        if (realARF == lpIntegralLowerBound) {
          if (realARF < minARF) {
            minARF = realARF;
            vector<pair<Slab, User*>> cellCombo;
            for (int u = 0; u < numUsers; ++u)
              cellCombo.push_back({users[u]->getSlab(integerRanks[u]), users[u]});
            bestCell = Cell(cellCombo);
            solution = wLP;
            logNewIncumbent("incumbent_exact_lp_gap", minARF, &current,
                            &innerQueue, &globalQueue);
            if (anytimeStopRequested) {
              break;
            }
          }
          if (anytimeStopRequested) {
            break;
          }
          continue;
        }

        if (!scoreWidthFitsMilp) {
          int fallbackUser = -1;
          int fallbackSplitSlab = -1;
          if (selectMidpointFallbackSplit(
                  users, current.minBorders, current.maxBorders,
                  fallbackUser, fallbackSplitSlab)) {
            bestGap = {fallbackUser, 0.0, fallbackSplitSlab};
            hasSplitCandidate = true;
          }
        }
      }

      // --- Hybrid trigger: switch to MILP for small search spaces ---
      int maxGap = hasSplitCandidate ? (int)bestGap.gap : 0;

      MILPModelEstimate milpEstimate;
      bool milpEligible = false;
      if (maxGap <= MILP_GAP_THRESHOLD &&
          scoreWidthFitsMilp) {
        milpEstimate = estimateMILPModelInFacet(
            users,
            current.minBorders,
            current.maxBorders,
            (int)facetList[current.facetIdx].getVertices().size(),
            (int)current.splitHalfspaces.size(),
            dim);
        milpEligible = true;
      }

      if (milpEligible) {
          // When the residual node is small enough, switch from LP-guided
          // branching to the exact MILP solver for a direct certification.
          Cell cplexResult;
          vector<double> cplexSolution;
          int cplexMin = minARF;
          MILPResult cplexStatus = findMinimumARF_InFacet(
              users,
              current.minBorders,
              current.maxBorders,
              facetList[current.facetIdx].getVertices(),
              current.splitHalfspaces,
              dim,
              cplexResult, cplexSolution, cplexMin, milpEstimate,
              current.nodeScoreMins.empty() ? nullptr : &current.nodeScoreMins,
              current.nodeScoreMaxs.empty() ? nullptr : &current.nodeScoreMaxs);

          if (cplexStatus == MILP_CERTIFIED && (cplexMin < minARF ||
              (cplexMin == minARF && solution.empty()))) {
              minARF = cplexMin;
              bestCell = cplexResult;
              solution = cplexSolution;
              logNewIncumbent("incumbent_milp", minARF, &current,
                              &innerQueue, &globalQueue);
              if (anytimeStopRequested) {
                break;
              }
          }
          if (anytimeStopRequested) {
            break;
          }
          if (cplexStatus == MILP_CERTIFIED || cplexStatus == MILP_INFEASIBLE) {
            continue; // Node is fully handled by MILP.
          }
      }

      if (!hasSplitCandidate) {
        int fallbackUser = -1;
        int fallbackSplitSlab = -1;
        if (selectMidpointFallbackSplit(
                users, current.minBorders, current.maxBorders,
                fallbackUser, fallbackSplitSlab)) {
          bestGap = {fallbackUser, 0.0, fallbackSplitSlab};
          hasSplitCandidate = true;
        }
      }

      if (!hasSplitCandidate) {
        int usersWithRange = 0;
        int totalRange = 0;
        int maxRangeRemaining = 0;
        for (int u = 0; u < numUsers; ++u) {
          int range = current.maxBorders[u] - current.minBorders[u];
          if (range > 0) {
            ++usersWithRange;
            totalRange += range;
            maxRangeRemaining = std::max(maxRangeRemaining, range);
          }
        }

        int lpIntegralLowerBound = (int)ceil(relaxedMinARF);
        bool allUsersFixed = (usersWithRange == 0);
        bool deadEndWouldResolve = (realARF == lpIntegralLowerBound);
        bool fallbackGateTriggered = !scoreWidthFitsMilp;

        std::ostringstream oss;
        oss << "Unresolved sub-region has no split candidate and did not qualify for MILP"
            << " | allUsersFixed=" << (allUsersFixed ? "yes" : "no")
            << " | usersWithRange=" << usersWithRange
            << " | totalRange=" << totalRange
            << " | maxRange=" << maxRangeRemaining
            << " | realARF=" << realARF
            << " | lpIntegralLowerBound=" << lpIntegralLowerBound
            << " | deadEndWouldResolve=" << (deadEndWouldResolve ? "yes" : "no")
            << " | avgScoreWidth=" << avgScoreWidth
            << " | milpEligible=" << (milpEligible ? "yes" : "no")
            << " | milpEntered=no"
            << " | fallbackGateTriggered=" << (fallbackGateTriggered ? "yes" : "no")
            << " | fallbackProducedSplit=no";
        throw std::runtime_error(oss.str());
      }

      vector<SubRegion> children(2);
      for (int branch = 0; branch < 2; ++branch) {
        children[branch].facetIdx = current.facetIdx;
        children[branch].minBorders = current.minBorders;
        children[branch].maxBorders = current.maxBorders;
        children[branch].splitHalfspaces = current.splitHalfspaces;
        children[branch].relaxedBasis = current.relaxedBasis;

        int uIdx = bestGap.userIndex;
        int splitSlab = bestGap.splitSlab;
        double threshold = users[uIdx]->getSlab(splitSlab).getLowerScore();
        const vector<double>& coords = users[uIdx]->getCoordinates();

        if (branch == 1) {
          // Branch "rank_u <= splitSlab".
          vector<double> negCoords(dim);
          for (int j = 0; j < dim; ++j) negCoords[j] = -coords[j];
          children[branch].splitHalfspaces.push_back(Halfspace(negCoords, -threshold));
          children[branch].maxBorders[uIdx] = splitSlab;
        } else {
          // Branch "rank_u >= splitSlab + 1".
          children[branch].splitHalfspaces.push_back(Halfspace(coords, threshold));
          children[branch].minBorders[uIdx] = splitSlab + 1;
        }

        if (!hasValidBorders(children[branch].minBorders, children[branch].maxBorders, users)) {
          children[branch].relaxedARF = 1e18;
          continue;
        }

        // Tighten each child with per-user LP score extrema before paying for
        // the full relaxed ARF solve.
        if (!tightenBordersByNodeScoreLP(
                users,
                children[branch].minBorders,
                children[branch].maxBorders,
                facetList[current.facetIdx].getVertices(),
                children[branch].splitHalfspaces,
                dim,
                borderLpSolveCount,
                borderLpSolveSeconds,
                0.0,
                &children[branch].nodeScoreMins,
                &children[branch].nodeScoreMaxs)) {
          children[branch].relaxedARF = 1e18;
        }
      }

      for (int branch = 0; branch < 2; ++branch) {
        if (children[branch].relaxedARF >= 1e17) {
          continue;
        }

        // Only children whose lower bound still beats the incumbent are worth
        // queuing for deeper search.
        solveRelaxedSubRegion(children[branch], facetList, users, dim,
                              relaxedLpSolveCount, relaxedLpSolveSeconds,
                              &facetHullCache);

        if (children[branch].relaxedARF < minARF) {
          innerQueue.push(children[branch]);
        }
      }
      maybeLogAnytimeInterval(minARF, nullptr, &innerQueue, &globalQueue);
      if (requestAnytimeStopIfReady("stop_after_children", minARF, nullptr,
                                    &innerQueue, &globalQueue)) {
        break;
      }
    }
  }

  if (anytimeStopRequested) {
    cout << "[ANYTIME-STOPPED] Returning incumbent before exact certification "
            "completed." << endl;
  }

  cout << "[RESULT] Slab combination: ";
  for (const pair<int, User *>& idx : bestCell.getCombinationIndex()) {
    cout << idx.first << " ";
  }
  cout << endl;

  cout << "[RESULT] Feasible point: " << std::setprecision(17);
  for (double val : solution) {
    cout << val << " ";
  }
  cout << endl;

  cout << "\nARF: " << (long long)bestCell.getARF() + rankOffset << endl;

  return finish();
}
