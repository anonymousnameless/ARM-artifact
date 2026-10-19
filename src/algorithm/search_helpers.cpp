#include "search_helpers.h"
#include "facet.h"
#include "user.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>

double userScoreAtPoint(User* user, const std::vector<double>& point, int dim) {
  const std::vector<double>& coords = user->getCoordinates();
  double score = 0.0;
  for (int d = 0; d < dim && d < (int)coords.size() &&
                  d < (int)point.size(); ++d) {
    score += coords[d] * point[d];
  }
  return score;
}

double averageScoreWidthInSlabUnits(const SubRegion& region,
                                    const std::vector<User*>& users,
                                    const Facet& facet,
                                    int dim) {
  if (users.empty()) {
    return 0.0;
  }

  const bool hasNodeScoreBounds =
      region.nodeScoreMins.size() == users.size() &&
      region.nodeScoreMaxs.size() == users.size();
  const std::vector<std::vector<double>>& vertices = facet.getVertices();

  double totalWidth = 0.0;
  for (int u = 0; u < (int)users.size(); ++u) {
    const int slabRange = region.maxBorders[u] - region.minBorders[u];
    if (slabRange <= 0) {
      continue;
    }

    double scoreMin = 0.0;
    double scoreMax = 0.0;
    if (hasNodeScoreBounds) {
      scoreMin = region.nodeScoreMins[u];
      scoreMax = region.nodeScoreMaxs[u];
    } else {
      scoreMin = std::numeric_limits<double>::infinity();
      scoreMax = -std::numeric_limits<double>::infinity();
      for (const std::vector<double>& vertex : vertices) {
        const double score = userScoreAtPoint(users[u], vertex, dim);
        scoreMin = std::min(scoreMin, score);
        scoreMax = std::max(scoreMax, score);
      }
    }

    if (!std::isfinite(scoreMin) || !std::isfinite(scoreMax)) {
      totalWidth += slabRange;
      continue;
    }
    if (scoreMin > scoreMax) {
      std::swap(scoreMin, scoreMax);
    }

    const double scoreWidth = std::max(0.0, scoreMax - scoreMin);
    const double borderScoreWidth = std::fabs(
        users[u]->getSlab(region.maxBorders[u]).getLowerScore() -
        users[u]->getSlab(region.minBorders[u]).getLowerScore());

    if (borderScoreWidth <= 1e-12) {
      totalWidth += (scoreWidth <= 1e-12) ? 0.0 : slabRange;
      continue;
    }

    const double localSlabSpacing = borderScoreWidth / slabRange;
    totalWidth += scoreWidth / localSlabSpacing;
  }

  return totalWidth / users.size();
}

void solveRelaxedSubRegion(
    SubRegion& region,
    const std::vector<Facet>& facetList,
    const std::vector<User*>& users,
    int dim,
    uint64_t& relaxedLpSolveCount,
    double& relaxedLpSolveSeconds,
    std::vector<UserHullCache>* hullCache) {
  const auto lpStart = std::chrono::steady_clock::now();

  const std::vector<std::vector<double>>& facetVerts =
      facetList[region.facetIdx].getVertices();

  double slabLowerBound = 0.0;
  for (int slabIdx : region.minBorders) {
    slabLowerBound += slabIdx;
  }

  double lpRelaxedARF =
      checkARFFeasibility_LP(users, region.minBorders, region.maxBorders,
                             facetVerts, region.splitHalfspaces, dim, region.wLP,
                             region.relaxedRanks, &region.relaxedBasis, hullCache);

  // Store the strongest safe lower bound once. Every pruning, queue ordering,
  // and integrality check below reads region.relaxedARF instead of recomputing it.
  region.relaxedARF = std::max(lpRelaxedARF, slabLowerBound);
  relaxedLpSolveSeconds +=
      std::chrono::duration<double>(std::chrono::steady_clock::now() - lpStart)
          .count();
  ++relaxedLpSolveCount;
}

long long computeExactPointArfFromProductFile(
    const char* productFile,
    const std::vector<User*>& users,
    const std::vector<double>& point,
    int dim) {
  std::ifstream input(productFile);
  if (!input.is_open()) {
    throw std::runtime_error("Failed to open product file for exact ARF");
  }

  std::vector<double> pointScores(users.size(), 0.0);
  for (int u = 0; u < (int)users.size(); ++u) {
    const std::vector<double>& coords = users[u]->getCoordinates();
    for (int d = 0; d < dim; ++d) {
      pointScores[u] += coords[d] * point[d];
    }
  }

  long long arf = 0;
  int id = 0;
  std::vector<double> lower(dim, 0.0), upper(dim, 0.0), center(dim, 0.0);
  while (input >> id) {
    for (int d = 0; d < dim; ++d) {
      input >> lower[d];
    }
    for (int d = 0; d < dim; ++d) {
      input >> upper[d];
      center[d] = 0.5 * (lower[d] + upper[d]);
    }
    if (!input) {
      throw std::runtime_error("Malformed product file while computing exact ARF");
    }

    for (int u = 0; u < (int)users.size(); ++u) {
      const std::vector<double>& coords = users[u]->getCoordinates();
      double productScore = 0.0;
      for (int d = 0; d < dim; ++d) {
        productScore += coords[d] * center[d];
      }
      if (productScore > pointScores[u] + 1e-12) {
        ++arf;
      }
    }
  }

  return arf;
}
