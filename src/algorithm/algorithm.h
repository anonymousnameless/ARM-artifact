#pragma once 

using namespace std;

#include <vector>

struct User;
struct Halfspace;
struct Cell;
struct Facet;

enum MILPResult {
    MILP_CERTIFIED,
    MILP_INFEASIBLE,
    MILP_TOO_LARGE
};

struct MILPModelEstimate {
    int choiceVars = 0;
    int cols = 0;
    int rows = 0;
};

MILPModelEstimate estimateMILPModelInFacet(
    const vector<User*>& users,
    const vector<int>& minBorders,
    const vector<int>& maxBorders,
    int numFacetVertices,
    int numSplitHalfspaces,
    int dimen);

MILPResult findMinimumARF_InFacet(
    vector<User*>& users,
    const vector<int>& minBorders,
    const vector<int>& maxBorders,
    const vector<vector<double>>& facetVertices,  // Root facet V-rep (λ combination)
    const vector<Halfspace>& splitHalfspaces,     // Accumulated split halfspaces
    int dimen,
    Cell& resultCell,
    vector<double>& solution,
    int& minARF,
    const MILPModelEstimate& sizeEstimate,
    const vector<double>* nodeScoreMins = nullptr,
    const vector<double>* nodeScoreMaxs = nullptr
);
