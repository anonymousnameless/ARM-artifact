#pragma once 

using namespace std;

#include <iostream>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <lp_lib.h>

struct RtreeNode;
struct User;
struct Inequality;
struct Cell;
struct Halfspace;
struct Facet;


void rtreePruning(User *u, RtreeNode* node, int dim, double minscore, double maxscore);

bool checkFeasible(const vector<pair<int,User*>> &combination, vector<User*> users, vector<User*> hullUsers, vector<Halfspace>& inequalities, int dimen, vector<double> &solution,bool hullOptimization);

bool checkFeasibleExact(const vector<pair<int,User*>> &combination,
                        vector<User*> users,
                        vector<Halfspace>& inequalities,
                        int dimen,
                        vector<double>& solution);

bool checkFeasibleInFacet(vector<pair<int,User*>>& combination,vector<User*> users, vector<User*> hullUsers,Facet& facet,int dimen, vector<double>& solution, bool hullOptimization);

void TrimWcaEdges(vector<std::vector<double>> intersections,vector<User*> hullUsers, vector<std::vector<double>> &TrimIntersections);

bool CheckCellDomination(vector<Cell>& infeasibleCells, Cell& element);

vector<Cell> CalculateOuterUBCCells(vector<vector<double>>& intersections,vector<User*>& users);

enum class ARFCheckMode {
    LESS_OR_EQUAL,    // Check if ARF <= target is feasible
    GREATER_OR_EQUAL  // Check if ARF >= target is feasible
};

struct LpBasisWarmStart {
    vector<int> colStatus;
    vector<int> rowStatus;

    bool empty() const {
        return colStatus.empty() || rowStatus.empty();
    }

    void clear() {
        colStatus.clear();
        rowStatus.clear();
    }
};

using UserHullCache = unordered_map<uint64_t, vector<pair<double, double>>>;

int canonicalizeSplitSlab(User* user,
                          int splitSlab,
                          int minBorder,
                          int maxBorder);

vector<int> collectCanonicalSplitSlabs(User* user,
                                       int minBorder,
                                       int maxBorder);

bool selectMidpointFallbackSplit(
    const vector<User*>& users,
    const vector<int>& minBorders,
    const vector<int>& maxBorders,
    int& userIndex,
    int& splitSlab);

// LP minimization: returns the minimum Relaxed ARF achievable within a sub-region.
// Hybrid approach:
//   - facetVertices: V-rep of the ROOT facet (convex combination, fixed per facet)
//   - splitHalfspaces: accumulated split halfspaces added one per split level
//     Each halfspace encodes one side of a user score split: coords·x <= threshold
//     (or -coords·x <= -threshold for the "keep above" side)
// Splitting is encoded via minBorders/maxBorders per user:
//   Child A (rank <= splitSlab): maxBorders[u] = splitSlab
//   Child B (rank >  splitSlab): minBorders[u] = splitSlab + 1
double checkARFFeasibility_LP(
    const vector<User*>& users,
    const vector<int>& minBorders,          // Min slab index per user
    const vector<int>& maxBorders,          // Max slab index per user
    const vector<vector<double>>& facetVertices,  // Root facet V-rep (λ combination)
    const vector<Halfspace>& splitHalfspaces,     // Accumulated split halfspaces
    int dimen,
    vector<double>& solution,               // Output: preference vector w
    vector<double>& relaxedRanks,           // Output: per-user relaxed rank s[i]
    LpBasisWarmStart* basisWarmStart = nullptr,
    vector<UserHullCache>* hullCache = nullptr
);

// Exact per-user border tightening:
// For each user i, solve two LPs on the node geometry to find the minimum and
// maximum achievable score, then convert those score extrema into tighter slab
// bounds and intersect them with [minBorders[i], maxBorders[i]].
// Returns false if the node becomes empty.
bool tightenBordersByNodeScoreLP(
    const vector<User*>& users,
    vector<int>& minBorders,
    vector<int>& maxBorders,
    const vector<vector<double>>& facetVertices,
    const vector<Halfspace>& splitHalfspaces,
    int dimen,
    uint64_t& lpSolveCount,
    double& lpSolveSeconds,
    double slabSnapSafetyMargin = 0.0,
    vector<double>* nodeScoreMins = nullptr,
    vector<double>* nodeScoreMaxs = nullptr);
