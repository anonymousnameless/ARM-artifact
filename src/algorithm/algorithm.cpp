#include "algorithm.h"
#include "../data/user.h"
#include "../utils/utils.h"
#include "optimization.h"
#include "../data/cell.h"
#include "../data/halfspace.h"
#include "../data/facet.h"
#include <cmath>
#include <iostream>

#include <ilcplex/ilocplex.h>

namespace {
constexpr int kMaxMilpChoiceVars = 600;
constexpr int kMaxMilpCols = 650;
constexpr int kMaxMilpRows = 1300;
vector<int> collectEffectiveMilpSlabs(User* user,
                                      int minBorder,
                                      int maxBorder,
                                      double tol = 1e-9) {
    vector<int> effectiveSlabs;
    if (minBorder > maxBorder) {
        return effectiveSlabs;
    }

    for (int slab = minBorder; slab <= maxBorder; ++slab) {
        if (effectiveSlabs.empty() ||
            std::fabs(user->getSlab(slab).getLowerScore() -
                      user->getSlab(effectiveSlabs.back()).getLowerScore()) > tol) {
            effectiveSlabs.push_back(slab);
        }
    }
    return effectiveSlabs;
}

}


MILPModelEstimate estimateMILPModelInFacet(
    const vector<User*>& users,
    const vector<int>& minBorders,
    const vector<int>& maxBorders,
    int numFacetVertices,
    int numSplitHalfspaces,
    int dimen) {
    MILPModelEstimate estimate;
    const int numUsers = (int)users.size();
    if ((int)minBorders.size() != numUsers || (int)maxBorders.size() != numUsers) {
        return estimate;
    }

    bool allFixed = true;
    int nonFixedUsers = 0;
    for (int i = 0; i < numUsers; ++i) {
        const vector<int> effectiveSlabs =
            collectEffectiveMilpSlabs(users[i], minBorders[i], maxBorders[i]);
        const int rangeSize = (int)effectiveSlabs.size();
        if (rangeSize <= 0) {
            return estimate;
        }
        if (rangeSize > 1) {
            allFixed = false;
            ++nonFixedUsers;
            estimate.choiceVars += rangeSize;
        }
    }

    if (allFixed) {
        estimate.choiceVars = 0;
        estimate.cols = dimen + numFacetVertices;
        estimate.rows = 1 + dimen + numSplitHalfspaces + 2 * numUsers;
        return estimate;
    }

    estimate.cols = dimen + numFacetVertices + estimate.choiceVars;
    const int fixedUsers = numUsers - nonFixedUsers;
    estimate.rows =
        2 * fixedUsers + nonFixedUsers + 2 * estimate.choiceVars + 1 +
        dimen + numSplitHalfspaces;
    return estimate;
}

MILPResult findMinimumARF_InFacet(
    vector<User*>& users,
    const vector<int>& minBorders,
    const vector<int>& maxBorders,
    const vector<vector<double>>& facetVertices,
    const vector<Halfspace>& splitHalfspaces,
    int dimen,
    Cell& resultCell,
    vector<double>& solution,
    int& minARF,
    const MILPModelEstimate& sizeEstimate,
    const vector<double>* nodeScoreMins,
    const vector<double>* nodeScoreMaxs
) {
    int num_users = users.size();
    if (num_users == 0 || facetVertices.empty()) return MILP_INFEASIBLE;
    if ((int)minBorders.size() != num_users ||
        (int)maxBorders.size() != num_users ||
        !hasValidBorders(minBorders, maxBorders, users)) {
        return MILP_INFEASIBLE;
    }

    const bool innerTooLarge =
        !(sizeEstimate.cols > 0 && sizeEstimate.rows > 0 &&
          sizeEstimate.choiceVars <= kMaxMilpChoiceVars &&
          sizeEstimate.cols <= kMaxMilpCols &&
          sizeEstimate.rows <= kMaxMilpRows);
    if (innerTooLarge) {
        return MILP_TOO_LARGE;
    }

    IloEnv env;
    try {
        vector<vector<int>> effectiveSlabs(num_users);
        bool allFixed = true;
        int fixedArf = 0;
        vector<pair<Slab, User*>> fixedCombination;
        fixedCombination.reserve(num_users);
        for (int i = 0; i < num_users; ++i) {
            effectiveSlabs[i] =
                collectEffectiveMilpSlabs(users[i], minBorders[i], maxBorders[i]);
            if (effectiveSlabs[i].empty()) {
                env.end();
                return MILP_INFEASIBLE;
            }
            if (effectiveSlabs[i].size() != 1) {
                allFixed = false;
                continue;
            }
            const int fixedSlab = effectiveSlabs[i][0];
            fixedArf += fixedSlab;
            fixedCombination.push_back({users[i]->getSlab(fixedSlab), users[i]});
        }

        if (allFixed) {
            IloModel model(env);
            IloNumVarArray x(env, dimen, 0.0, 1.0, ILOFLOAT);

            const vector<vector<double>>& verts = facetVertices;
            int nVerts = (int)verts.size();
            IloNumVarArray lambda(env, nVerts, 0.0, 1.0, ILOFLOAT);
            IloExpr sumLambda(env);
            for (int i = 0; i < nVerts; ++i) sumLambda += lambda[i];
            model.add(sumLambda == 1.0);
            sumLambda.end();

            for (int d = 0; d < dimen; ++d) {
                IloExpr val(env);
                for (int i = 0; i < nVerts; ++i) val += lambda[i] * verts[i][d];
                model.add(x[d] == val);
                val.end();
            }

            for (const auto& hs : splitHalfspaces) {
                IloExpr lhs(env);
                const auto& c = const_cast<Halfspace&>(hs).getCoords();
                for (int d = 0; d < dimen; ++d) {
                    lhs += c[d] * x[d];
                }
                model.add(lhs <= const_cast<Halfspace&>(hs).getScore());
                lhs.end();
            }

            for (int i = 0; i < num_users; ++i) {
                IloExpr dot(env);
                for (int d = 0; d < dimen; ++d) {
                    dot += users[i]->getCoordinates()[d] * x[d];
                }

                const int slabIdx = minBorders[i];
                const double lower = users[i]->getSlab(slabIdx).getLowerScore();
                const double upper =
                    (slabIdx > 0) ? users[i]->getSlab(slabIdx - 1).getLowerScore()
                                  : 1.0;
                model.add(dot >= lower);
                model.add(dot <= upper);
                dot.end();
            }

            IloExpr zeroObjective(env);
            model.add(IloMinimize(env, zeroObjective));
            zeroObjective.end();

            IloCplex cplex(model);
            cplex.setOut(env.getNullStream());
            cplex.setWarning(env.getNullStream());
            cplex.setParam(IloCplex::Param::Emphasis::Numerical, 1);
            cplex.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-9);
            cplex.setParam(IloCplex::Param::Threads, 1);

            if (cplex.solve()) {
                solution.resize(dimen);
                for (int d = 0; d < dimen; ++d) {
                    solution[d] = cplex.getValue(x[d]);
                }
                resultCell.clean();
                resultCell = Cell(fixedCombination);
                minARF = fixedArf;
                env.end();
                return MILP_CERTIFIED;
            }

            auto status = cplex.getStatus();
            if (status == IloAlgorithm::Infeasible) {
                env.end();
                return MILP_INFEASIBLE;
            }

            throw std::runtime_error("CPLEX failed unexpectedly in all-fixed findMinimumARF_InFacet");
        }

        IloModel model(env);
        
        // --- Variables ---
        IloNumVarArray x(env, dimen, 0.0, 1.0, ILOFLOAT);
        vector<IloBoolVarArray> choice(num_users);
        vector<int> rangeSize(num_users);
        vector<double> scoreMin(num_users);
        vector<double> scoreMax(num_users);
        const bool hasStoredNodeScoreBounds =
            nodeScoreMins != nullptr && nodeScoreMaxs != nullptr &&
            (int)nodeScoreMins->size() == num_users &&
            (int)nodeScoreMaxs->size() == num_users;
        for (int i = 0; i < num_users; ++i) {
            rangeSize[i] = (int)effectiveSlabs[i].size();
            if (rangeSize[i] > 1) {
                choice[i] = IloBoolVarArray(env, rangeSize[i]);
            }

            double boundedMin = users[i]->getSlab(maxBorders[i]).getLowerScore();
            double boundedMax =
                (minBorders[i] > 0)
                    ? users[i]->getSlab(minBorders[i] - 1).getLowerScore()
                    : 1.0;

            if (hasStoredNodeScoreBounds) {
                boundedMin = std::max(boundedMin, (*nodeScoreMins)[i]);
                boundedMax = std::min(boundedMax, (*nodeScoreMaxs)[i]);
                if (boundedMin > boundedMax + 1e-9) {
                    env.end();
                    return MILP_INFEASIBLE;
                }
            }
            scoreMin[i] = boundedMin;
            scoreMax[i] = boundedMax;
        }
        
        // --- Objective: Minimize ARF ---
        IloExpr targetSum(env);
        for (int i = 0; i < num_users; ++i) {
            if (rangeSize[i] == 1) {
                continue;
            }
            for (int k = 0; k < rangeSize[i]; ++k) {
                int actualSlab = effectiveSlabs[i][k];
                targetSum += actualSlab * choice[i][k];
            }
        }
        model.add(IloMinimize(env, targetSum));
        targetSum.end();
        
        // --- Constraints ---
        // 1. One slab per user
        for (int i = 0; i < num_users; ++i) {
            if (rangeSize[i] == 1) {
                continue;
            }
            IloExpr oneSlab(env);
            for (int k = 0; k < rangeSize[i]; ++k) oneSlab += choice[i][k];
            model.add(oneSlab == 1);
            oneSlab.end();
        }
        
        // 2. Slab geometry with safe big-M constraints.
        for (int i = 0; i < num_users; ++i) {
            IloExpr dot(env);
            for (int d = 0; d < dimen; ++d) dot += users[i]->getCoordinates()[d] * x[d];

            if (rangeSize[i] == 1) {
                const int fixedSlab = effectiveSlabs[i][0];
                const double lower = users[i]->getSlab(fixedSlab).getLowerScore();
                const double upper =
                    (fixedSlab > 0) ? users[i]->getSlab(fixedSlab - 1).getLowerScore()
                                    : 1.0;
                model.add(dot >= lower);
                model.add(dot <= upper);
                dot.end();
                continue;
            }
            
            for (int k = 0; k < rangeSize[i]; ++k) {
                int actualSlab = effectiveSlabs[i][k];
                const double lower = users[i]->getSlab(actualSlab).getLowerScore();
                const double upper =
                    (actualSlab > 0) ? users[i]->getSlab(actualSlab - 1).getLowerScore()
                                     : 1.0;

                model.add(dot >= lower - (lower - scoreMin[i]) * (1 - choice[i][k]));
                model.add(dot <= upper + (scoreMax[i] - upper) * (1 - choice[i][k]));
            }
            dot.end();
        }
        
        // 3. Sub-region Geometry: hybrid V-rep base + split halfspace rows
        const vector<vector<double>>& verts = facetVertices;
        int nVerts = (int)verts.size();
        IloNumVarArray lambda(env, nVerts, 0.0, 1.0, ILOFLOAT);
        IloExpr sumLambda(env);
        for(int i=0; i<nVerts; ++i) sumLambda += lambda[i];
        model.add(sumLambda == 1.0);
        sumLambda.end();

        for(int d=0; d<dimen; ++d) {
            IloExpr val(env);
            for(int i=0; i<nVerts; ++i) val += lambda[i] * verts[i][d];
            model.add(x[d] == val);
            val.end();
        }

        // Add accumulated split halfspace constraints
        for (const auto& hs : splitHalfspaces) {
            IloExpr lhs(env);
            const auto& c = const_cast<Halfspace&>(hs).getCoords();
            for (int d = 0; d < dimen; ++d) {
                lhs += c[d] * x[d];
            }
            model.add(lhs <= const_cast<Halfspace&>(hs).getScore());
            lhs.end();
        }

        // --- Solve ---
        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        cplex.setParam(IloCplex::Param::Emphasis::MIP, 1);
        cplex.setParam(IloCplex::Param::Emphasis::Numerical, 1);
        cplex.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-9);
        cplex.setParam(IloCplex::Param::MIP::Tolerances::Integrality, 1e-9);
        cplex.setParam(IloCplex::Param::Threads, 1);

        const bool solved = cplex.solve();
        if (solved) {
            // --- Extract solution ---
            solution.resize(dimen);
            for(int k=0; k<dimen; ++k) solution[k] = cplex.getValue(x[k]);

            vector<pair<Slab, User*>> chosenCombination;
            int chosenArf = 0;
            for (int i = 0; i < num_users; ++i) {
                if (rangeSize[i] == 1) {
                    const int fixedSlab = effectiveSlabs[i][0];
                    chosenArf += fixedSlab;
                    chosenCombination.push_back({users[i]->getSlab(fixedSlab), users[i]});
                    continue;
                }

                int chosenK = -1;
                for (int k = 0; k < rangeSize[i]; ++k) {
                    if (cplex.getValue(choice[i][k]) > 0.5) {
                        if (chosenK >= 0) {
                            env.end();
                            return MILP_TOO_LARGE;
                        }
                        chosenK = k;
                    }
                }

                if (chosenK < 0) {
                    double bestValue = -1.0;
                    for (int k = 0; k < rangeSize[i]; ++k) {
                        const double value = cplex.getValue(choice[i][k]);
                        if (value > bestValue) {
                            bestValue = value;
                            chosenK = k;
                        }
                    }
                }

                if (chosenK < 0) {
                    env.end();
                    return MILP_TOO_LARGE;
                }

                const int chosenSlab = effectiveSlabs[i][chosenK];
                chosenArf += chosenSlab;
                chosenCombination.push_back({users[i]->getSlab(chosenSlab), users[i]});
            }

            resultCell.clean();
            resultCell = Cell(chosenCombination);
            minARF = chosenArf;
            env.end();
            return MILP_CERTIFIED;
        }
        
        auto status = cplex.getStatus();
        if (status == IloAlgorithm::Infeasible) {
            env.end();
            return MILP_INFEASIBLE;
        }

        throw std::runtime_error("CPLEX failed unexpectedly in findMinimumARF_InFacet");
        
    } catch (IloException& e) {
        cerr << "[MILP] CPLEX Error: " << e.getMessage() << endl;
        env.end();
        throw std::runtime_error("CPLEX Exception in findMinimumARF_InFacet: " + string(e.getMessage()));
    } catch (...) {
        cerr << "[MILP] Unknown error" << endl;
        env.end();
        throw;
    }
}





