#include "utils.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <climits> 
#include <string>
#include <stdexcept>
#include <thread>

#include "../spatial/rnode.h"
#include "../spatial/rtree.h"
#include "../spatial/rentry.h"
#include "../data/user.h"
#include "../data/facet.h"
#include "../data/halfspace.h"
#include "../spatial/hypercube.h"
#include "../data/slab.h"
#include <list>

namespace {
uint64_t readStatusValueBytes(const char *key) {
#if defined(__linux__)
    std::ifstream status("/proc/self/status");
    std::string label;
    uint64_t value = 0;
    std::string unit;
    while (status >> label >> value >> unit) {
        if (label == key) {
            return value * 1024ULL;
        }
    }
#else
    (void)key;
#endif
    return 0;
}

}

const char* Param::read(const int a_argc, const char** a_argv, const char* a_param, const char* a_def)
{
	for (int i = 0; i < a_argc; i++)
	{
		if (strcmp(a_argv[i], a_param) == 0)
		{
			if (i + 1 == a_argc)
				return "";
			else
				return a_argv[i + 1];
		}
	}
    
    return a_def;
}

uint64_t getProcessPeakRssBytes() {
    uint64_t highWaterMark = readStatusValueBytes("VmHWM:");
    if (highWaterMark != 0) {
        return highWaterMark;
    }
    return readStatusValueBytes("VmRSS:");
}

uint64_t getAvailableMemoryBytes() {
#if defined(__linux__)
    std::ifstream meminfo("/proc/meminfo");
    std::string label;
    uint64_t value = 0;
    std::string unit;
    while (meminfo >> label >> value >> unit) {
        if (label == "MemAvailable:") {
            return value * 1024ULL;
        }
    }
#endif
    return 0;
}

int getCpuSolveCap() {
    unsigned int hwThreads = std::max(1u, std::thread::hardware_concurrency());
    int reserveThreads = std::max(1, (int)hwThreads / 4);
    return std::max(1, (int)hwThreads - reserveThreads);
}

uint64_t addSolveSafetyMargin(uint64_t bytes) {
    return bytes + bytes / 2;
}

uint64_t getFallbackLpJobBytes() {
    return 256ULL * 1024ULL * 1024ULL;
}

int getMemorySolveCap(uint64_t availableMemoryBytes,
                      uint64_t estimatedLpJobBytes) {
    if (availableMemoryBytes == 0 || estimatedLpJobBytes == 0) {
        return INT_MAX;
    }

    uint64_t usableBytes = availableMemoryBytes * 60ULL / 100ULL;
    if (usableBytes < estimatedLpJobBytes) {
        return 1;
    }

    uint64_t cap = usableBytes / estimatedLpJobBytes;
    if (cap > (uint64_t)INT_MAX) {
        cap = (uint64_t)INT_MAX;
    }
    return std::max(1, (int)cap);
}

void appendFixedWCA(vector<Halfspace>& inequalities, int dim) {
    const double rhs = 0.16;
    vector<vector<double>> coeffs;

    switch (dim) {
    case 3:
        coeffs = {
            {0.29711254, 0.22283441, 0.14855627},
            {0.12977714, 0.32444284, 0.19466571},
            {0.19466571, 0.12977714, 0.32444284},
        };
        break;
    case 4:
        coeffs = {
            {0.28571429, 0.22857143, 0.11428571, 0.11428571},
            {0.11428607, 0.11428607, 0.34285821, 0.12777182},
            {0.17142420, 0.28570700, 0.17142420, 0.13999643},
        };
        break;
    case 6:
        coeffs = {
            {0.11428776, 0.17143164, 0.22857552, 0.05714388, 0.17143164,
             0.18068895},
            {0.28571605, 0.05714321, 0.05714321, 0.22857284, 0.11428642,
             0.08080050},
            {0.17142857, 0.17142857, 0.17142857, 0.17142857, 0.17142857,
             0.11428571},
        };
        break;
    default:
        throw std::invalid_argument(
            "Fixed WCA is only defined for d=3, d=4, and d=6");
    }

    for (const auto& coeff : coeffs) {
        inequalities.push_back({coeff, rhs});
    }
}

double dotProduct(const vector<double>& a, const vector<double>& b) {
    double value = 0.0;
    for (int i = 0; i < (int)a.size(); ++i) {
        value += a[i] * b[i];
    }
    return value;
}

int findContainingSlabIndex(User* user, double score, double tol) {
    const int slabCount = user->getSlabSize();
    const int baseIdx = user->findClosestSlab(score).getIndex();

    for (int delta = -1; delta <= 1; ++delta) {
        const int idx = baseIdx + delta;
        if (idx < 0 || idx >= slabCount) {
            continue;
        }

        const Slab& slab = user->getSlab(idx);
        if (score >= slab.getLowerScore() - tol &&
            score <= slab.getUpperScore() + tol) {
            return idx;
        }
    }

    return baseIdx;
}

int computePointArf(const vector<User*>& users,
                    const vector<double>& point,
                    vector<int>& slabIndices) {
    slabIndices.resize(users.size());
    int arf = 0;
    for (int i = 0; i < (int)users.size(); ++i) {
        const double score = dotProduct(users[i]->getCoordinates(), point);
        const int slabIdx = findContainingSlabIndex(users[i], score);
        slabIndices[i] = slabIdx;
        arf += slabIdx;
    }
    return arf;
}

void eraseUsersInPlace(vector<User*>& users, const vector<User*>& toErase) {
    if (toErase.empty()) {
        return;
    }

    users.erase(std::remove_if(users.begin(), users.end(),
                               [&](User* user) {
                                   return std::find(toErase.begin(), toErase.end(), user) !=
                                          toErase.end();
                               }),
                users.end());
}

Rtree* loadProducts(const char* productFile, int dim)
{
    std::fstream fpdata(productFile, std::ios::in);
    if (!fpdata.is_open()) {
        std::cerr << "Failed to open file: " << productFile << "\n";
        return nullptr;
    }

    int objCnt = 0;
    std::vector<RtreeNodeEntry*> entries;
    entries.reserve(131072);

    try {
        while (true) {
            int id;
            fpdata >> id;
            if (!fpdata) break;

            std::vector<float> cl(dim), cu(dim);
            for (int d = 0; d < dim; ++d) fpdata >> cl[d];
            for (int d = 0; d < dim; ++d) fpdata >> cu[d];

            Hypercube hc(dim, cl.data(), cu.data());
            entries.push_back(new RtreeNodeEntry(id, hc));
            ++objCnt;

            if (objCnt >= MAXPTS) {
                std::cerr << "Warning: Reached MAXPTS limit (" << MAXPTS << ")\n";
                break;
            }
        }
        fpdata.close();

        if (objCnt == 0) {
            std::cerr << "No entries loaded from file\n";
            return nullptr;
        }

        int maxChild = computeMaxChild(4096, dim);
        Rtree* rtree = new Rtree(dim, maxChild, maxChild, 
                                (int)(maxChild*0.3f), (int)(maxChild*0.3f));
        rtree->bulkLoad(entries.data(), objCnt, maxChild);
        return rtree;
        
    } catch (...) {
        // Clean up on exception
        for (RtreeNodeEntry* entry : entries) {
            delete entry;
        }
        throw;
    }
}

static int computeMaxChild(int pageSize, int dimen) {
    int nodeHeaderSize = sizeof(bool) + 3 * sizeof(int) + sizeof(RtreeNode*);
    int entrySize = Hypercube::size(dimen) + sizeof(int) + sizeof(RtreeNode*);
    return (pageSize - nodeHeaderSize) / entrySize;
}


void loadUsers(vector<User*> &users,list<User> &usersList, const char* userFile,int dim){
    fstream fpdata;
    int id;
    vector<double> coordinates;
    string tmpcord;

    fpdata.open(userFile, ios::in);
    
    while (true)
    {   
        fpdata >> id;
        if (!fpdata) break; 
        
        for (int d = 0; d < dim; d++){
            fpdata >> tmpcord;
            coordinates.push_back(floor(stod(tmpcord) * 10000.0) / 10000.0);
        }

       User user(coordinates);
    
       usersList.push_back(user);
       users.push_back(&usersList.back());
       coordinates.clear();
    }

}

void CreateWCA(vector<Halfspace>& inequalities, int d, int N, double r){
    CreateWCA(inequalities, d, N, r, 1);
}

void CreateWCA(vector<Halfspace>& inequalities,
               int d,
               int N,
               double r,
               uint64_t seed) {
    std::mt19937_64 gen(seed);
    std::normal_distribution<> dist(0.0, 1.0);
    
    inequalities.reserve(inequalities.size() + N);

    for(int i = 0; i < N; ++i) {
        std::vector<double> z(d);
        double norm = 0.0;

        /// Generate ONLY positive values from the start
        for(int j = 0; j < d; ++j) {
            z[j] = std::abs(dist(gen));  // Take abs of normal sample
            norm += z[j] * z[j];
        }

        norm = std::sqrt(norm);

        // Normalize to radius r (already positive)
        for(int j = 0; j < d; ++j) {
            z[j] = z[j] * r / norm;
        }

        // Tangent hyperplane: z·x ≤ r²
        inequalities.push_back({z, r*r});
    }

}


double computeARF(vector<int>& indices) {
    double arf = 0;
    for (int idx : indices) {
        arf += idx;
    }
    return arf;
}

int CheckCoordinateDomination(vector<double> a, vector<double> b){
    int betterA = 0;
    int betterB = 0;

    for (int i = 0; i < a.size(); ++i) {
        if (a[i] >  b[i]) {
            betterA ++; 
        }
        if (a[i] < b[i]) {
            betterB ++;
        }
    }
    if (betterA > 0 && betterB == 0) return 1; // a dominates b
    if (betterB > 0 && betterA == 0) return -1; // b dominates a
    return 0; // no domination
}

// ============================================================================
// Vertex-tracking border tightening helpers
// ============================================================================

vector<vector<double>> clipVertices(
    const vector<vector<double>>& vertices,
    const vector<double>& userCoords,
    double threshold,
    bool keepAbove,
    int dim,
    bool& verticesComplete)
{
    vector<vector<double>> result;
    vector<int> sideKeep, sideDiscard;
    const double eps = 1e-9;

    // Classify vertices
    for (int v = 0; v < (int)vertices.size(); ++v) {
        double score = 0.0;
        for (int d = 0; d < dim; ++d)
            score += userCoords[d] * vertices[v][d];
        bool onKeepSide = keepAbove ? (score >= threshold - eps) : (score <= threshold + eps);
        if (onKeepSide)
            sideKeep.push_back(v);
        else
            sideDiscard.push_back(v);
    }

    // Add kept vertices
    for (int k : sideKeep)
        result.push_back(vertices[k]);

    // Guard: if too many pairs, skip intersection computation to prevent O(V^2) blowup.
    // Mark vertices as INCOMPLETE — caller must NOT use them for border tightening.
    static const size_t MAX_PAIRS = 200;
    if (sideKeep.size() * sideDiscard.size() > MAX_PAIRS) {
        verticesComplete = false;
        return result;
    }

    // Compute intersection points for all pairs crossing the boundary
    for (int k : sideKeep) {
        double scoreK = 0.0;
        for (int d = 0; d < dim; ++d)
            scoreK += userCoords[d] * vertices[k][d];

        for (int di : sideDiscard) {
            double scoreD = 0.0;
            for (int d = 0; d < dim; ++d)
                scoreD += userCoords[d] * vertices[di][d];

            double denom = scoreD - scoreK;
            if (fabs(denom) < 1e-15) continue;

            double t = (threshold - scoreK) / denom;
            if (t < -eps || t > 1.0 + eps) continue;
            t = max(0.0, min(1.0, t));

            vector<double> newPt(dim);
            for (int d = 0; d < dim; ++d)
                newPt[d] = vertices[k][d] + t * (vertices[di][d] - vertices[k][d]);
            result.push_back(newPt);
        }
    }

    return result;
}

void deduplicateVertices(vector<vector<double>>& vertices, double eps) {
    vector<vector<double>> unique;
    for (const auto& v : vertices) {
        bool isDup = false;
        for (const auto& u : unique) {
            double dist2 = 0.0;
            for (int d = 0; d < (int)v.size(); ++d) {
                double diff = v[d] - u[d];
                dist2 += diff * diff;
            }
            if (dist2 < eps * eps) {
                isDup = true;
                break;
            }
        }
        if (!isDup)
            unique.push_back(v);
    }
    vertices = unique;
}

bool recomputeBorders(
    const vector<vector<double>>& subVertices,
    vector<User*>& users,
    vector<int>& minBorders,
    vector<int>& maxBorders)
{
    int numUsers = users.size();
    if (subVertices.empty()) return false;

    for (int u = 0; u < numUsers; ++u) {
        int newMin = INT_MAX, newMax = -1;
        for (const auto& v : subVertices) {
            double score = users[u]->CalculateScore(v);
            Slab& slab = users[u]->findClosestSlab(score);
            int rank = slab.getIndex();
            newMin = min(newMin, rank);
            newMax = max(newMax, rank);
        }
        minBorders[u] = max(minBorders[u], newMin);
        maxBorders[u] = min(maxBorders[u], newMax);
        if (minBorders[u] > maxBorders[u])
            return false;
    }
    return true;
}

bool hasValidBorders(const vector<int>& minBorders,
                     const vector<int>& maxBorders,
                     const vector<User*>& users)
{
    if (minBorders.size() != maxBorders.size() || minBorders.size() != users.size()) {
        return false;
    }

    for (size_t i = 0; i < minBorders.size(); ++i) {
        if (minBorders[i] > maxBorders[i]) {
            return false;
        }
        if (minBorders[i] < 0 || maxBorders[i] >= users[i]->getSlabSize()) {
            return false;
        }
    }

    return true;
}
