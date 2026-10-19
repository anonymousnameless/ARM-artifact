#pragma once 

using namespace std;

#include <cstdint>
#include <climits>
#include <vector>
#include <list>

struct Rtree; 
struct User;
struct RtreeNode;
struct Halfspace;
struct Facet;

struct Param {
    static const char* read(const int argc, const char** argv, const char* param, const char* def);
};

double computeARF(vector<int>& indices);

void appendFixedWCA(vector<Halfspace>& inequalities, int dim);

double dotProduct(const vector<double>& a, const vector<double>& b);

Rtree* loadProducts(const char* productFile,int dim);

void loadUsers(vector<User*> &usersLists, list<User> &tmpusersLists, const char* userFile,int dim);

void CreateWCA(vector<Halfspace>& inequalities,int d, int N, double r);

void CreateWCA(vector<Halfspace>& inequalities,
               int d,
               int N,
               double r,
               uint64_t seed);

int CheckCoordinateDomination(vector<double> a, vector<double> b);

static int computeMaxChild(int pageSize, int dimen);

// Vertex-tracking border tightening helpers
vector<vector<double>> clipVertices(
    const vector<vector<double>>& vertices,
    const vector<double>& userCoords,
    double threshold,
    bool keepAbove,
    int dim,
    bool& verticesComplete);

void deduplicateVertices(vector<vector<double>>& vertices, double eps = 1e-10);

bool recomputeBorders(
    const vector<vector<double>>& subVertices,
    vector<User*>& users,
    vector<int>& minBorders,
    vector<int>& maxBorders);

bool hasValidBorders(const vector<int>& minBorders,
                     const vector<int>& maxBorders,
                     const vector<User*>& users);

uint64_t getProcessPeakRssBytes();

uint64_t getAvailableMemoryBytes();

int getCpuSolveCap();

uint64_t addSolveSafetyMargin(uint64_t bytes);

uint64_t getFallbackLpJobBytes();

int getMemorySolveCap(uint64_t availableMemoryBytes,
                      uint64_t estimatedLpJobBytes);

int findContainingSlabIndex(User* user, double score, double tol = 1e-9);

int computePointArf(const vector<User*>& users,
                    const vector<double>& point,
                    vector<int>& slabIndices);

void eraseUsersInPlace(vector<User*>& users, const vector<User*>& toErase);
