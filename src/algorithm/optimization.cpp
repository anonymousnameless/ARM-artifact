#include "optimization.h"
#include "../utils/utils.h"
#include "../data/user.h"
#include "../data/facet.h"
#include "../spatial/rnode.h"
#include <lp_lib.h>
#include "../spatial/rentry.h"
#include "../data/cell.h"
#include "../data/slab.h"
#include <ilcplex/ilocplex.h>
#include <chrono>
#include <cstring>
#include <gmpxx.h>
#include <algorithm>

using namespace std::chrono;

namespace {

using ExactInt = mpz_class;
using ExactRat = mpq_class;

struct ExactHalfspaceConstraint {
    vector<ExactRat> coeffs;
    ExactRat rhs;
};

ExactRat exactRatFromDouble(double value)
{
    if (value == 0.0) {
        return ExactRat(0);
    }

    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));

    const bool negative = (bits >> 63) != 0;
    const int rawExponent = (int)((bits >> 52) & 0x7ffU);
    const std::uint64_t rawMantissa = bits & ((1ULL << 52) - 1ULL);

    ExactInt mantissa = 0;
    int exponent = 0;

    if (rawExponent == 0) {
        mantissa = rawMantissa;
        exponent = 1 - 1023 - 52;
    } else {
        mantissa = (ExactInt(1) << 52) + rawMantissa;
        exponent = rawExponent - 1023 - 52;
    }

    ExactRat result(mantissa);
    if (exponent > 0) {
        result *= ExactRat(ExactInt(1) << exponent);
    } else if (exponent < 0) {
        result /= ExactRat(ExactInt(1) << (-exponent));
    }

    if (negative) {
        result = -result;
    }
    result.canonicalize();
    return result;
}

double exactRatToDouble(const ExactRat& value)
{
    return value.get_d();
}

bool solveExactSystem(vector<vector<ExactRat>> matrix,
                      vector<ExactRat> rhs,
                      vector<ExactRat>& solution)
{
    const int n = (int)matrix.size();
    solution.assign(n, ExactRat(0));
    vector<int> pivotColumn(n, -1);

    int row = 0;
    for (int col = 0; col < n && row < n; ++col) {
        int pivot = -1;
        for (int r = row; r < n; ++r) {
            if (matrix[r][col] != 0) {
                pivot = r;
                break;
            }
        }
        if (pivot == -1) {
            continue;
        }

        if (pivot != row) {
            std::swap(matrix[pivot], matrix[row]);
            std::swap(rhs[pivot], rhs[row]);
        }

        ExactRat pivotValue = matrix[row][col];
        for (int c = col; c < n; ++c) {
            matrix[row][c] /= pivotValue;
        }
        rhs[row] /= pivotValue;

        for (int r = 0; r < n; ++r) {
            if (r == row || matrix[r][col] == 0) {
                continue;
            }
            ExactRat factor = matrix[r][col];
            for (int c = col; c < n; ++c) {
                matrix[r][c] -= factor * matrix[row][c];
            }
            rhs[r] -= factor * rhs[row];
        }

        pivotColumn[row] = col;
        ++row;
    }

    if (row != n) {
        return false;
    }

    for (int r = 0; r < n; ++r) {
        if (pivotColumn[r] == -1) {
            return false;
        }
        solution[pivotColumn[r]] = rhs[r];
    }
    return true;
}

bool satisfiesExactConstraints(const vector<ExactHalfspaceConstraint>& constraints,
                               const vector<ExactRat>& point)
{
    for (const auto& constraint : constraints) {
        ExactRat lhs = 0;
        for (int j = 0; j < (int)point.size(); ++j) {
            lhs += constraint.coeffs[j] * point[j];
        }
        if (lhs > constraint.rhs) {
            return false;
        }
    }
    return true;
}

void buildExactConstraints(const vector<pair<int, User*>>& combination,
                           vector<Halfspace>& inequalities,
                           int dimen,
                           vector<ExactHalfspaceConstraint>& constraints)
{
    constraints.clear();

    auto addConstraint = [&](const vector<double>& coeffs, double rhs) {
        ExactHalfspaceConstraint c;
        c.coeffs.reserve(dimen);
        for (int j = 0; j < dimen; ++j) {
            c.coeffs.push_back(exactRatFromDouble(coeffs[j]));
        }
        c.rhs = exactRatFromDouble(rhs);
        constraints.push_back(std::move(c));
    };

    for (const auto& entry : combination) {
        int slabIndex = entry.first;
        User* user = entry.second;
        const vector<double>& coords = user->getCoordinates();

        vector<double> negCoords(dimen);
        for (int j = 0; j < dimen; ++j) {
            negCoords[j] = -coords[j];
        }
        addConstraint(negCoords, -user->getSlab(slabIndex).getLowerScore());

        if (slabIndex != 0) {
            addConstraint(coords, user->getSlab(slabIndex - 1).getLowerScore());
        }
    }

    for (const auto& ineq : inequalities) {
        addConstraint(ineq.getCoords(), ineq.getScore());
    }

    for (int j = 0; j < dimen; ++j) {
        vector<double> upper(dimen, 0.0);
        upper[j] = 1.0;
        addConstraint(upper, 1.0);

        vector<double> lower(dimen, 0.0);
        lower[j] = -1.0;
        addConstraint(lower, 0.0);
    }
}

bool enumerateExactVertices(const vector<ExactHalfspaceConstraint>& constraints,
                            int dimen,
                            int startIndex,
                            vector<int>& active,
                            vector<double>& solution)
{
    if ((int)active.size() == dimen) {
        vector<vector<ExactRat>> matrix(dimen, vector<ExactRat>(dimen, ExactRat(0)));
        vector<ExactRat> rhs(dimen, ExactRat(0));
        for (int r = 0; r < dimen; ++r) {
            const auto& c = constraints[active[r]];
            matrix[r] = c.coeffs;
            rhs[r] = c.rhs;
        }

        vector<ExactRat> exactPoint;
        if (!solveExactSystem(matrix, rhs, exactPoint)) {
            return false;
        }

        if (!satisfiesExactConstraints(constraints, exactPoint)) {
            return false;
        }

        solution.resize(dimen);
        for (int j = 0; j < dimen; ++j) {
            solution[j] = exactRatToDouble(exactPoint[j]);
        }
        return true;
    }

    for (int idx = startIndex;
         idx <= (int)constraints.size() - (dimen - (int)active.size()); ++idx) {
        active.push_back(idx);
        if (enumerateExactVertices(constraints, dimen, idx + 1, active, solution)) {
            return true;
        }
        active.pop_back();
    }
    return false;
}

bool hasMatchingBasis(const LpBasisWarmStart* basisWarmStart,
                      int colCount,
                      int rowCount)
{
    return basisWarmStart != nullptr &&
           (int)basisWarmStart->colStatus.size() == colCount &&
           (int)basisWarmStart->rowStatus.size() == rowCount;
}

void applyBasisWarmStart(IloCplex& cplex,
                         IloEnv env,
                         const IloNumVarArray& vars,
                         const IloConstraintArray& rows,
                         const LpBasisWarmStart* basisWarmStart)
{
    if (!hasMatchingBasis(basisWarmStart, vars.getSize(), rows.getSize())) {
        return;
    }

    IloCplex::BasisStatusArray colStatus(env, vars.getSize());
    IloCplex::BasisStatusArray rowStatus(env, rows.getSize());

    for (int i = 0; i < vars.getSize(); ++i) {
        colStatus[i] =
            static_cast<IloCplex::BasisStatus>(basisWarmStart->colStatus[i]);
    }
    for (int i = 0; i < rows.getSize(); ++i) {
        rowStatus[i] =
            static_cast<IloCplex::BasisStatus>(basisWarmStart->rowStatus[i]);
    }

    cplex.setBasisStatuses(colStatus, vars, rowStatus, rows);
}

void captureBasisWarmStart(IloCplex& cplex,
                           IloEnv env,
                           const IloNumVarArray& vars,
                           const IloConstraintArray& rows,
                           LpBasisWarmStart* basisWarmStart)
{
    if (basisWarmStart == nullptr) {
        return;
    }

    IloCplex::BasisStatusArray colStatus(env, vars.getSize());
    IloCplex::BasisStatusArray rowStatus(env, rows.getSize());
    cplex.getBasisStatuses(colStatus, vars, rowStatus, rows);

    basisWarmStart->colStatus.assign(vars.getSize(), 0);
    basisWarmStart->rowStatus.assign(rows.getSize(), 0);

    for (int i = 0; i < vars.getSize(); ++i) {
        basisWarmStart->colStatus[i] = (int)colStatus[i];
    }
    for (int i = 0; i < rows.getSize(); ++i) {
        basisWarmStart->rowStatus[i] = (int)rowStatus[i];
    }
}

bool optimizeUserScoreOnNodeLP(
    User* user,
    const vector<vector<double>>& facetVertices,
    const vector<Halfspace>& splitHalfspaces,
    int dimen,
    bool maximize,
    double& scoreValue,
    LpBasisWarmStart* basisWarmStart = nullptr)
{
    if (facetVertices.empty()) {
        return false;
    }

    IloEnv env;
    try {
        IloModel model(env);

        IloNumVarArray x(env, dimen, 0.0, 1.0, ILOFLOAT);

        const vector<vector<double>>& vertices = facetVertices;
        int numVerts = (int)vertices.size();
        IloNumVarArray lambda(env, numVerts, 0.0, 1.0, ILOFLOAT);
        IloNumVarArray basisVars(env);
        for (int j = 0; j < dimen; ++j) {
            basisVars.add(x[j]);
        }
        for (int i = 0; i < numVerts; ++i) {
            basisVars.add(lambda[i]);
        }
        IloConstraintArray basisRows(env);

        IloExpr sumLambda(env);
        for (int i = 0; i < numVerts; ++i) {
            sumLambda += lambda[i];
        }
        IloConstraint sumLambdaConstraint = (sumLambda == 1.0);
        model.add(sumLambdaConstraint);
        sumLambda.end();

        for (int j = 0; j < dimen; ++j) {
            IloExpr conv(env);
            for (int i = 0; i < numVerts; ++i) {
                conv += lambda[i] * vertices[i][j];
            }
            IloConstraint convConstraint = (x[j] == conv);
            model.add(convConstraint);
            basisRows.add(convConstraint);
            conv.end();
        }

        for (const auto& hs : splitHalfspaces) {
            IloExpr lhs(env);
            const auto& c = const_cast<Halfspace&>(hs).getCoords();
            for (int j = 0; j < dimen; ++j) {
                lhs += c[j] * x[j];
            }
            IloConstraint splitConstraint =
                (lhs <= const_cast<Halfspace&>(hs).getScore());
            model.add(splitConstraint);
            basisRows.add(splitConstraint);
            lhs.end();
        }

        IloExpr score(env);
        const auto& coords = user->getCoordinates();
        for (int j = 0; j < dimen; ++j) {
            score += coords[j] * x[j];
        }

        IloObjective objective =
            maximize ? IloMaximize(env, score) : IloMinimize(env, score);
        model.add(objective);
        score.end();

        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        cplex.setParam(IloCplex::Param::Threads, 1);
        cplex.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-7);
        cplex.setParam(IloCplex::Param::RootAlgorithm, CPX_ALG_DUAL);
        applyBasisWarmStart(cplex, env, basisVars, basisRows, basisWarmStart);

        bool feasible = cplex.solve();
        if (feasible) {
            scoreValue = cplex.getObjValue();
            captureBasisWarmStart(cplex, env, basisVars, basisRows, basisWarmStart);
            env.end();
            return true;
        }

        auto status = cplex.getStatus();
        env.end();
        if (status == IloAlgorithm::Infeasible) {
            return false;
        }

        throw std::runtime_error(
            "LP solver failed unexpectedly in optimizeUserScoreOnNodeLP");
    } catch (IloException& e) {
        string message = e.getMessage();
        env.end();
        throw std::runtime_error(
            "CPLEX Exception in optimizeUserScoreOnNodeLP: " + message);
    } catch (...) {
        env.end();
        throw;
    }
}

} // namespace

bool checkFeasibleExact(const vector<pair<int,User*>> &combination,
                        vector<User*> users,
                        vector<Halfspace>& inequalities,
                        int dimen,
                        vector<double>& solution)
{
    (void)users;
    vector<ExactHalfspaceConstraint> constraints;
    buildExactConstraints(combination, inequalities, dimen, constraints);

    vector<int> active;
    return enumerateExactVertices(constraints, dimen, 0, active, solution);
}

bool checkFeasible(const vector<pair<int,User*>>  &combination, vector<User*> users, vector<User*> hullUsers, vector<Halfspace>& inequalities, int dimen, vector<double> &solution, bool hullOptimization) {
    
    solution.resize(dimen, 0.0); // prepare output

    try {
        IloEnv env;
        IloModel model(env);
        
        
        // --- 1. Decision variables: continuous x[0..dimen-1] in [0,1] ---
        IloNumVarArray x(env, dimen, 0.0, 1.0, ILOFLOAT);
        
        // --- 2. Add slab constraints based on combination ---
        if (hullOptimization == false) {
            // Add constraints for all users
            for (int i = 0; i < users.size(); ++i) {
                double score = users[i]->getSlab(combination[i].first).getLowerScore();
                const vector<double>& coordinates = users[i]->getCoordinates();
                
                // LOWER bound: dot(coordinates, x) >= score
                IloExpr lowerBound(env);
                for (int j = 0; j < dimen; ++j) {
                    lowerBound += coordinates[j] * x[j];
                }
                model.add(lowerBound >= score);
                lowerBound.end();
                
                // UPPER bound (if not the first slab)
                if (combination[i].first != 0) {
                    score = users[i]->getSlab(combination[i].first - 1).getLowerScore();
                    
                    IloExpr upperBound(env);
                    for (int j = 0; j < dimen; ++j) {
                        upperBound += coordinates[j] * x[j];
                    }
                    model.add(upperBound <= score);
                    upperBound.end();
                }
            }
        }
        else {
            // Hull optimization mode: only add constraints for users in hullUsers
            for (int i = 0; i < combination.size(); ++i) {
                auto it = find(hullUsers.begin(), hullUsers.end(), combination[i].second);
                if (it != hullUsers.end()) {
                    int index = distance(hullUsers.begin(), it);
                    double score = hullUsers[index]->getSlab(combination[i].first).getLowerScore();
                    const vector<double>& coordinates = hullUsers[index]->getCoordinates();
                    
                    // LOWER bound: dot(coordinates, x) >= score
                    IloExpr lowerBound(env);
                    for (int j = 0; j < dimen; ++j) {
                        lowerBound += coordinates[j] * x[j];
                    }
                    model.add(lowerBound >= score);
                    lowerBound.end();
                    
                    // UPPER bound (if not the first slab)
                    if (combination[i].first != 0) {
                        score = hullUsers[index]->getSlab(combination[i].first - 1).getLowerScore();
                        
                        IloExpr upperBound(env);
                        for (int j = 0; j < dimen; ++j) {
                            upperBound += coordinates[j] * x[j];
                        }
                        model.add(upperBound <= score);
                        upperBound.end();
                    }
                }
            }
        }
        
        // --- 3. Add global WCA inequalities ---
        for (auto& ineq : inequalities) {
            IloExpr wcaExpr(env);
            const auto& coeffs = ineq.getCoords();
            for (int j = 0; j < dimen; ++j) {
                wcaExpr += coeffs[j] * x[j];
            }
            model.add(wcaExpr <= ineq.getScore());
            wcaExpr.end();
        }
        
        // --- 4. Dummy objective (feasibility problem) ---
        model.add(IloMinimize(env, 0));
        
        // --- 5. Create solver and configure ---
        IloCplex cplex(model);
        cplex.setOut(env.getNullStream()); // Silent mode
        cplex.setParam(IloCplex::Param::Threads, 1);
        // Adjust tolerances to match lpsolve settings
        cplex.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-7);
        cplex.setParam(IloCplex::Param::Simplex::Tolerances::Optimality, 1e-7);
        
        // Quick solve - no need for aggressive parameters since this is LP
        cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 1);
        
        // --- 6. Solve ---
        bool feasible = cplex.solve();
        
        if (feasible) {
            // Extract solution
            for (int k = 0; k < dimen; ++k) {
                solution[k] = cplex.getValue(x[k]);
            }
        }
        
        // Clean up
        env.end();
        return feasible;
        
    } catch (IloException& e) {
        cerr << "CPLEX Error: " << e.getMessage() << endl;
        return false;
    }
}

bool checkFeasibleInFacet(vector<pair<int,User*>>& combination,vector<User*> users, vector<User*> hullUsers,Facet& facet,int dimen, vector<double>& solution, bool hullOptimization) 
{
    solution.resize(dimen, 0.0);

    try {
        IloEnv env;
        IloModel model(env);
        
        // --- 1. Decision variables: continuous x[0..dimen-1] in [0,1] ---
        IloNumVarArray x(env, dimen, 0.0, 1.0, ILOFLOAT);
        
        // --- 2. Add slab constraints based on combination ---
        if (hullOptimization == false) {
            for (int i = 0; i < users.size(); ++i) {
                double score = users[i]->getSlab(combination[i].first).getLowerScore();
                const vector<double>& coordinates = users[i]->getCoordinates();
                
                // LOWER bound
                IloExpr lowerBound(env);
                for (int j = 0; j < dimen; ++j) {
                    lowerBound += coordinates[j] * x[j];
                }
                model.add(lowerBound >= score);
                lowerBound.end();
                
                // UPPER bound (if not the first slab)
                if (combination[i].first != 0) {
                    score = users[i]->getSlab(combination[i].first - 1).getLowerScore();
                    
                    IloExpr upperBound(env);
                    for (int j = 0; j < dimen; ++j) {
                        upperBound += coordinates[j] * x[j];
                    }
                    model.add(upperBound <= score);
                    upperBound.end();
                }
            }
        }
        else {
            for (int i = 0; i < combination.size(); ++i) {
                auto it = find(hullUsers.begin(), hullUsers.end(), combination[i].second);
                if (it != hullUsers.end()) {
                    int index = distance(hullUsers.begin(), it);
                    double score = hullUsers[index]->getSlab(combination[i].first).getLowerScore();
                    const vector<double>& coordinates = hullUsers[index]->getCoordinates();
                    
                    IloExpr lowerBound(env);
                    for (int j = 0; j < dimen; ++j) {
                        lowerBound += coordinates[j] * x[j];
                    }
                    model.add(lowerBound >= score);
                    lowerBound.end();
                    
                    if (combination[i].first != 0) {
                        score = hullUsers[index]->getSlab(combination[i].first - 1).getLowerScore();
                        
                        IloExpr upperBound(env);
                        for (int j = 0; j < dimen; ++j) {
                            upperBound += coordinates[j] * x[j];
                        }
                        model.add(upperBound <= score);
                        upperBound.end();
                    }
                }
            }
        }
        
        // --- 3. Add convex combination constraints for facet vertices ---
        // Get vertices from your Facet class
        const vector<vector<double>>& facet_vertices = facet.getVertices();
        int num_vertices = facet_vertices.size();
        
        IloNumVarArray lambda(env, num_vertices, 0.0, 1.0, ILOFLOAT);
        
        // Constraint: Σλᵢ = 1
        IloExpr sumLambda(env);
        for (int i = 0; i < num_vertices; ++i) {
            sumLambda += lambda[i];
        }
        IloConstraint sumLambdaConstraint = (sumLambda == 1.0);
        model.add(sumLambdaConstraint);
        sumLambda.end();
        
        // Constraints: x[j] = Σ λᵢ * vertex[i][j] for each dimension j
        for (int j = 0; j < dimen; ++j) {
            IloExpr convexComb(env);
            for (int i = 0; i < num_vertices; ++i) {
                convexComb += lambda[i] * facet_vertices[i][j];
            }
            model.add(x[j] == convexComb);
            convexComb.end();
        }
        
        // --- 4. Dummy objective (feasibility problem) ---
        model.add(IloMinimize(env, 0));
        
        // --- 5. Create solver and configure ---
        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());
        
        cplex.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-7);
        cplex.setParam(IloCplex::Param::Simplex::Tolerances::Optimality, 1e-7);
        cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 1);
        
        // --- 6. Solve ---
        bool feasible = cplex.solve();
        
        if (feasible) {
            for (int k = 0; k < dimen; ++k) {
                solution[k] = cplex.getValue(x[k]);
            }
        }
        
        env.end();
        return feasible;
        
    } catch (IloException& e) {
        cerr << "CPLEX Error: " << e.getMessage() << endl;
        return false;
    }
}

static long long subtreeProductCount(RtreeNode* node) {
    if (node->cachedProductCount >= 0) return node->cachedProductCount;
    long long count = 0;
    for (int i = 0; i < node->m_count; ++i) {
        auto* entry = node->m_entries[i];
        count += entry->isLeaf() ? 1 : subtreeProductCount(entry->m_child);
    }
    return node->cachedProductCount = count;
}

void rtreePruning(User *u, RtreeNode* node, int dim, double minscore, double maxscore) {
    if (!node) return;

    // Cache user coordinates once — avoids repeated getCoordinates() calls
    const vector<double>& coords = u->getCoordinates();

    for (int i = 0; i < node->m_count; i++) {
        RtreeNodeEntry* entry = node->m_entries[i];

        const Point& lower = entry->m_mbr.getLower();
        const Point& upper = entry->m_mbr.getUpper();

        // Compute dot products directly from raw arrays — no vector allocation
        double lowerScore = 0.0, upperScore = 0.0;
        for (int j = 0; j < dim; j++) {
            lowerScore += coords[j] * lower.m_coor[j];
            upperScore += coords[j] * upper.m_coor[j];
        }

        if (lowerScore > maxscore) {
            u->addRankOffset(entry->isLeaf() ? 1 : subtreeProductCount(entry->m_child));
            continue;
        }
        if (upperScore < minscore) continue;

        if (entry->isLeaf()) {
            // Compute center score in-place: dot(coords, (lower + upper) / 2)
            // Avoids allocating a temporary vector<double> for every leaf
            double centerScore = 0.0;
            for (int j = 0; j < dim; j++) {
                centerScore += coords[j] * (lower.m_coor[j] + upper.m_coor[j]) * 0.5;
            }
            u->addProductScore(centerScore);
        } else {
            rtreePruning(u, entry->m_child, dim, minscore, maxscore);
        }
    }
}


void TrimWcaEdges(vector<std::vector<double>> intersections,vector<User*> hullUsers, vector<std::vector<double>> &TrimIntersections){

    struct ScoredIntersection {
        std::vector<double> score;
        int index;
    };

    std::vector<ScoredIntersection> scoredIntersections;
    scoredIntersections.reserve(intersections.size());

    for (int i = 0; i < (int)intersections.size(); ++i) {
        std::vector<double> currentScore(hullUsers.size(), 0.0);
        for (int j = 0; j < (int)hullUsers.size(); ++j) {
            currentScore[j] = hullUsers[j]->CalculateScore(intersections[i]);
        }

        scoredIntersections.push_back({currentScore, i});
    }

    std::sort(scoredIntersections.begin(), scoredIntersections.end(),
              [](const ScoredIntersection& a, const ScoredIntersection& b) {
                  if (a.score != b.score) {
                      return a.score > b.score;
                  }
                  return a.index < b.index;
              });

    std::vector<std::vector<double>> keptScores;

    for (const ScoredIntersection& current : scoredIntersections) {
        bool isDominated = false;

        for (int j = 0; j < (int)keptScores.size(); ++j) {
            int result = CheckCoordinateDomination(keptScores[j], current.score);
            if (result == 1) { // existing dominates current
                isDominated = true;
                break;
            }
        }

        if (!isDominated) {
            keptScores.push_back(current.score);
            TrimIntersections.push_back(intersections[current.index]);
        }
    }
}

bool CheckCellDomination(vector<Cell>& infeasibleCells, Cell& element) {
    for ( auto& c : infeasibleCells) {

        if(element.weaklyDominates(c))
            return true; // found a dominating Cell
    }

    return false; // no dominating Cell found
}

vector<Cell> CalculateOuterUBCCells(vector<vector<double>>& intersections,vector<User*>& users){
    vector<Cell> OuterUBCCells;
    vector<pair<Slab,User*>> OuterSlabCombination;
    for(int i=0;i<intersections.size();i++){
        for(int j=0;j<users.size();j++){
            double score = users[j]->CalculateScore(intersections[i]);
            int index = -1;
            users[j]->findClosestSlab(score, &index);
            if(index == 0)
                OuterSlabCombination.push_back({users[j]->getSlab(0),users[j]});
            else
                OuterSlabCombination.push_back({users[j]->getSlab(index-1),users[j]});
        }
        OuterUBCCells.push_back({OuterSlabCombination});
        OuterSlabCombination.clear();
    }
    return OuterUBCCells;
}

namespace {

// Helper function to compute lower convex hull of 2D points (x, y)
// Assumes points are sorted by x-coordinate
// Returns hull vertices in left-to-right order
vector<pair<double, double>> computeLowerConvexHull(vector<pair<double, double>>& points) {
    if (points.size() <= 2) return points;
    
    // Sort by x-coordinate
    sort(points.begin(), points.end());
    
    vector<pair<double, double>> hull;
    
    for (const auto& p : points) {
        // Remove points that make a right turn (not on lower hull)
        while (hull.size() >= 2) {
            auto& p1 = hull[hull.size() - 2];
            auto& p2 = hull[hull.size() - 1];
            
            // Compute cross product to check if we make a left turn
            double cross = (p2.first - p1.first) * (p.second - p1.second) - 
                          (p2.second - p1.second) * (p.first - p1.first);
            
            // If right turn or collinear, remove p2
            if (cross >= 0) {
                hull.pop_back();
            } else {
                break;
            }
        }
        hull.push_back(p);
    }
    
    return hull;
}

// Helper function to compute upper convex hull of 2D points (x, y)
// Assumes points are sorted by x-coordinate
// Returns hull vertices in left-to-right order
vector<pair<double, double>> computeUpperConvexHull(vector<pair<double, double>>& points) {
    if (points.size() <= 2) return points;
    
    // Sort by x-coordinate
    sort(points.begin(), points.end());
    
    vector<pair<double, double>> hull;
    
    for (const auto& p : points) {
        // Remove points that make a left turn (not on upper hull)
        while (hull.size() >= 2) {
            auto& p1 = hull[hull.size() - 2];
            auto& p2 = hull[hull.size() - 1];
            
            // Compute cross product to check if we make a right turn
            double cross = (p2.first - p1.first) * (p.second - p1.second) - 
                          (p2.second - p1.second) * (p.first - p1.first);
            
            // If left turn or collinear, remove p2
            if (cross <= 0) {
                hull.pop_back();
            } else {
                break;
            }
        }
        hull.push_back(p);
    }
    
    return hull;
}

uint64_t makeHullCacheKey(int minBorder, int maxBorder) {
    return (uint64_t(uint32_t(minBorder)) << 32) | uint32_t(maxBorder);
}

bool sameHullThreshold(double lhs, double rhs, double tol = 1e-12) {
    return std::fabs(lhs - rhs) <= tol;
}

bool sameLowerThreshold(User* user, int lhsSlab, int rhsSlab,
                        double tol = 1e-12) {
    return std::fabs(user->getSlab(lhsSlab).getLowerScore() -
                     user->getSlab(rhsSlab).getLowerScore()) <= tol;
}

} // namespace

int canonicalizeSplitSlab(User* user,
                          int splitSlab,
                          int minBorder,
                          int maxBorder) {
    int canonical = std::max(minBorder, std::min(splitSlab, maxBorder));
    while (canonical + 1 <= maxBorder &&
           sameLowerThreshold(user, canonical, canonical + 1)) {
        ++canonical;
    }
    return canonical;
}

vector<int> collectCanonicalSplitSlabs(User* user, int minBorder, int maxBorder) {
    vector<int> candidates;
    if (minBorder >= maxBorder) {
        return candidates;
    }

    int slab = minBorder;
    while (slab <= maxBorder) {
        int canonical = canonicalizeSplitSlab(user, slab, minBorder, maxBorder);
        if (canonical < maxBorder) {
            candidates.push_back(canonical);
        }
        slab = canonical + 1;
    }

    return candidates;
}

bool selectMidpointFallbackSplit(
    const vector<User*>& users,
    const vector<int>& minBorders,
    const vector<int>& maxBorders,
    int& userIndex,
    int& splitSlab) {
    struct MidpointSplit {
        int userIndex;
        int range;
        int splitSlab;
    };

    MidpointSplit bestMidpoint{-1, -1, -1};

    for (int u = 0; u < (int)users.size(); ++u) {
        const int range = maxBorders[u] - minBorders[u];
        if (range <= 0) {
            continue;
        }

        vector<int> splitCandidates =
            collectCanonicalSplitSlabs(users[u], minBorders[u], maxBorders[u]);
        if (splitCandidates.empty()) {
            continue;
        }

        const int midSlab = splitCandidates[splitCandidates.size() / 2];
        if (range > bestMidpoint.range) {
            bestMidpoint = {u, range, midSlab};
        }
    }

    if (bestMidpoint.userIndex == -1) {
        return false;
    }

    userIndex = bestMidpoint.userIndex;
    splitSlab = bestMidpoint.splitSlab;
    return true;
}

namespace {

vector<pair<double, double>> buildNormalizedHull(User* user,
                                                 int minBorder,
                                                 int maxBorder) {
    vector<pair<double, double>> points;
    points.reserve((size_t)maxBorder - (size_t)minBorder + 1);

    const double rankMin = (double)minBorder;
    const double rankRange = (double)(maxBorder - minBorder);
    for (int k = minBorder; k <= maxBorder; ++k) {
        double slabScore = user->getSlab(k).getLowerScore();
        double normalizedRank = ((double)k - rankMin) / rankRange;
        if (!points.empty() &&
            sameHullThreshold(points.back().first, slabScore)) {
            continue;
        }
        points.push_back({slabScore, normalizedRank});
    }

    return computeUpperConvexHull(points);
}

} // namespace

double checkARFFeasibility_LP(
    const vector<User*>& users,
    const vector<int>& minBorders,
    const vector<int>& maxBorders,
    const vector<vector<double>>& facetVertices,
    const vector<Halfspace>& splitHalfspaces,
    int dimen,
    vector<double>& solution,
    vector<double>& relaxedRanks,
    LpBasisWarmStart* basisWarmStart,
    vector<UserHullCache>* hullCache)
{
    solution.resize(dimen, 0.0);  // Prepare output for preference vector w
    int numUsers = users.size();
    relaxedRanks.resize(numUsers, 0.0);  // Prepare output for per-user relaxed ranks

    // Guard: empty facet vertex set is trivially infeasible
    if (facetVertices.empty()) return 1e9;
    if ((int)minBorders.size() != numUsers ||
        (int)maxBorders.size() != numUsers ||
        !hasValidBorders(minBorders, maxBorders, users)) {
        return 1e9;
    }

    try {
        IloEnv env;
        IloModel model(env);
        
        // Decision variables: x[0..dimen-1] in [0,1] (preference vector w)
        IloNumVarArray x(env, dimen, 0.0, 1.0, ILOFLOAT);
        IloNumVarArray basisVars(env);
        for (int j = 0; j < dimen; ++j) {
            basisVars.add(x[j]);
        }
        IloConstraintArray basisRows(env);
        
        // ========== FACET V-REP CONSTRAINT (base region) ==========
        // x must be a convex combination of the ROOT FACET's vertices.
        // This correctly constrains x to the WCA facet geometry.
        const vector<vector<double>>& vertices = facetVertices;
        int numVerts = (int)vertices.size();

        IloNumVarArray lambda(env, numVerts, 0.0, 1.0, ILOFLOAT);

        // Σλᵢ = 1
        IloExpr sumLambda(env);
        for (int i = 0; i < numVerts; ++i) {
            sumLambda += lambda[i];
        }
        model.add(sumLambda == 1.0);
        sumLambda.end();

        // x[j] = Σλᵢ * vertex[i][j]
        for (int j = 0; j < dimen; ++j) {
            IloExpr conv(env);
            for (int i = 0; i < numVerts; ++i) {
                conv += lambda[i] * vertices[i][j];
            }
            IloConstraint convConstraint = (x[j] == conv);
            model.add(convConstraint);
            basisRows.add(convConstraint);
            conv.end();
        }

        // ========== SPLIT HALFSPACE CONSTRAINTS ==========
        // Each split adds one halfspace: coords·x <= threshold (or -coords·x <= -threshold).
        // These incrementally refine the sub-region without touching the λ variables.
        for (const auto& hs : splitHalfspaces) {
            IloExpr lhs(env);
            const auto& c = const_cast<Halfspace&>(hs).getCoords();
            for (int j = 0; j < dimen; ++j)
                lhs += c[j] * x[j];
            IloConstraint splitConstraint =
                (lhs <= const_cast<Halfspace&>(hs).getScore());
            model.add(splitConstraint);
            basisRows.add(splitConstraint);
            lhs.end();
        }
        
        // ========== ARF CONSTRAINT (LP RELAXATION WITH NORMALIZED CONVEX HULL) ==========
        // For users with non-zero range, use a normalized rank variable in [0,1]:
        //   rank_i = minBorder_i + range_i * sNorm_i
        // This keeps the objective in original slab-index scale while shrinking the
        // hull coefficients substantially compared with raw rank variables.
        vector<IloNumVar> normalizedRanks(numUsers);
        vector<double> rankMins(numUsers, 0.0);
        vector<double> rankRanges(numUsers, 0.0);
        vector<bool> hasRankVar(numUsers, false);

        for (int i = 0; i < numUsers; ++i) {
            rankMins[i] = (double)minBorders[i];
            rankRanges[i] = (double)(maxBorders[i] - minBorders[i]);

            if (rankRanges[i] <= 0.0) {
                // User is already pinned to a single slab; no LP rank variable is needed.
                continue;
            }

            normalizedRanks[i] = IloNumVar(env, 0.0, 1.0, ILOFLOAT);
            hasRankVar[i] = true;
            basisVars.add(normalizedRanks[i]);

            User* u = users[i];
            const auto& coords = u->getCoordinates();

            // score = dot(coords, x)
            IloExpr score(env);
            for (int j = 0; j < dimen; ++j) {
                score += coords[j] * x[j];
            }

            const vector<pair<double, double>>* hullPtr = nullptr;
            vector<pair<double, double>> localHull;
            if (hullCache != nullptr && i < (int)hullCache->size()) {
                const uint64_t cacheKey =
                    makeHullCacheKey(minBorders[i], maxBorders[i]);
                UserHullCache& userCache = (*hullCache)[i];
                auto cachedHull = userCache.find(cacheKey);
                if (cachedHull == userCache.end()) {
                    cachedHull = userCache
                                     .emplace(cacheKey,
                                              buildNormalizedHull(
                                                  u, minBorders[i], maxBorders[i]))
                                     .first;
                }
                hullPtr = &cachedHull->second;
            } else {
                localHull = buildNormalizedHull(u, minBorders[i], maxBorders[i]);
                hullPtr = &localHull;
            }

            const vector<pair<double, double>>& hull = *hullPtr;

            for (size_t j = 0; j + 1 < hull.size(); ++j) {
                double x1 = hull[j].first;
                double y1 = hull[j].second;
                double x2 = hull[j + 1].first;
                double y2 = hull[j + 1].second;

                double m = (y2 - y1) / (x2 - x1 + 1e-10);
                double c = y1 - m * x1;

                // sNorm[i] >= m * score + c
                IloConstraint hullConstraint = (normalizedRanks[i] >= m * score + c);
                model.add(hullConstraint);
                basisRows.add(hullConstraint);
            }

            score.end();
        }

        // Objective stays in original slab-index scale.
        IloExpr arfSum(env);
        for (int i = 0; i < numUsers; ++i) {
            arfSum += rankMins[i];
            if (hasRankVar[i]) {
                arfSum += rankRanges[i] * normalizedRanks[i];
            }
        }
        
        // OBJECTIVE: Minimize ARF Sum
        model.add(IloMinimize(env, arfSum));
        arfSum.end();
        
        // Solve LP
        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        cplex.setParam(IloCplex::Param::Threads, 1);
        cplex.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-7);
        cplex.setParam(IloCplex::Param::RootAlgorithm, CPX_ALG_DUAL);
        applyBasisWarmStart(cplex, env, basisVars, basisRows, basisWarmStart);
        bool feasible = cplex.solve();

        if (feasible) {
            double objValue = cplex.getObjValue();
            
            // Extract preference vector w
            for (int j = 0; j < dimen; ++j) {
                solution[j] = cplex.getValue(x[j]);
            }
            
            // Extract per-user relaxed ranks in the original slab-index scale.
            for (int i = 0; i < numUsers; ++i) {
                if (hasRankVar[i]) {
                    const double normalizedValue = cplex.getValue(normalizedRanks[i]);
                    relaxedRanks[i] =
                        rankMins[i] + rankRanges[i] * normalizedValue;
                } else {
                    relaxedRanks[i] = rankMins[i];
                }
            }

            captureBasisWarmStart(cplex, env, basisVars, basisRows, basisWarmStart);

            cplex.end();
            env.end();
            return objValue;
        }

        auto status = cplex.getStatus();
        if (status == IloAlgorithm::Infeasible) {
            cplex.end();
            env.end();
            return 1e9;
        }

        throw std::runtime_error("LP solver failed unexpectedly in checkARFFeasibility_LP");
        
    } catch (IloException& e) {
        cerr << "CPLEX Error: " << e.getMessage() << endl;
        throw std::runtime_error("CPLEX Exception in checkARFFeasibility_LP: " + string(e.getMessage()));
    }
}

bool tightenBordersByNodeScoreLP(
    const vector<User*>& users,
    vector<int>& minBorders,
    vector<int>& maxBorders,
    const vector<vector<double>>& facetVertices,
    const vector<Halfspace>& splitHalfspaces,
    int dimen,
    uint64_t& lpSolveCount,
    double& lpSolveSeconds,
    double slabSnapSafetyMargin,
    vector<double>* nodeScoreMins,
    vector<double>* nodeScoreMaxs)
{
    if ((int)minBorders.size() != (int)users.size() ||
        (int)maxBorders.size() != (int)users.size() ||
        !hasValidBorders(minBorders, maxBorders, users)) {
        return false;
    }

    if (nodeScoreMins != nullptr) {
        nodeScoreMins->assign(users.size(), 0.0);
    }
    if (nodeScoreMaxs != nullptr) {
        nodeScoreMaxs->assign(users.size(), 0.0);
    }

    LpBasisWarmStart basisWarmStart;

    for (int u = 0; u < (int)users.size(); ++u) {
        double minScore = 0.0;
        double maxScore = 0.0;

        auto minSolveStart = steady_clock::now();
        bool minFeasible = optimizeUserScoreOnNodeLP(
            users[u], facetVertices, splitHalfspaces, dimen, false, minScore,
            &basisWarmStart);
        lpSolveSeconds += duration<double>(steady_clock::now() - minSolveStart).count();
        ++lpSolveCount;
        if (!minFeasible) {
            return false;
        }

        auto maxSolveStart = steady_clock::now();
        bool maxFeasible = optimizeUserScoreOnNodeLP(
            users[u], facetVertices, splitHalfspaces, dimen, true, maxScore,
            &basisWarmStart);
        lpSolveSeconds += duration<double>(steady_clock::now() - maxSolveStart).count();
        ++lpSolveCount;
        if (!maxFeasible) {
            return false;
        }

        if (nodeScoreMins != nullptr) {
            (*nodeScoreMins)[u] = minScore;
        }
        if (nodeScoreMaxs != nullptr) {
            (*nodeScoreMaxs)[u] = maxScore;
        }

        // Be conservative when translating floating-point LP score extrema into
        // discrete slab bounds. Expanding the reachable score interval avoids
        // cutting off near-threshold exact-feasible cells due to LP tolerance.
        int newMinRank =
            users[u]->findClosestSlab(maxScore + slabSnapSafetyMargin).getIndex();
        int newMaxRank =
            users[u]->findClosestSlab(minScore - slabSnapSafetyMargin).getIndex();

        minBorders[u] = max(minBorders[u], newMinRank);
        maxBorders[u] = min(maxBorders[u], newMaxRank);
        if (minBorders[u] > maxBorders[u]) {
            return false;
        }
    }

    return true;
}

