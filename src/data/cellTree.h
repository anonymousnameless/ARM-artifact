#ifndef _CELL_TREE_H_
#define _CELL_TREE_H_

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <map>
#include <set>
#include "lp_lib.h"

// Forward declarations
struct User;
struct Facet;
struct Halfspace;
class Cell;

typedef struct cell {
    // History of splits: ID -> Side (true=Left/Satisfied, false=Right/Unsatisfied)
    std::unordered_map<long int, bool> IntersectHP;
    
    // Fixed constraints
    std::unordered_set<long int> AboveHP; // Fixed Left
    std::unordered_set<long int> BelowHP; // Fixed Right
    
    int arf;           // Incremental ARF (relative to minBorders)
    int realarf;       // Absolute ARF (sum of actual slab indices)
    bool isPruned;
    cell *left, *right;
    
    cell() : left(nullptr), right(nullptr), isPruned(false), arf(0), realarf(0) {}
    cell(cell* copy) : left(nullptr), right(nullptr), arf(copy->arf), realarf(copy->realarf), isPruned(copy->isPruned) {}
    
    // CORRECTED: Merges history into 'node' instead of overwriting
    void copyleaf(cell* node) {
        for (auto& kv : IntersectHP) {
            node->IntersectHP[kv.first] = kv.second;
        }
        for (auto id : BelowHP) node->BelowHP.insert(id);
        for (auto id : AboveHP) node->AboveHP.insert(id);
        node->arf = arf;
        node->realarf = realarf;
        node->isPruned = isPruned;
    }
    
    void release() {
        IntersectHP.clear(); AboveHP.clear(); BelowHP.clear();
        std::unordered_map<long int, bool>().swap(IntersectHP);
        std::unordered_set<long int>().swap(AboveHP);
        std::unordered_set<long int>().swap(BelowHP);
    }
} cell;

struct HPInfo {
    User* user;
    int userIndex, slabIndex;
};

class cellTree {
private:
    cell* root;
    int dimen;
    Facet* rootFacet;
    
    // Mathematical definitions
    std::vector<std::vector<float>> HalfSpaces;
    std::vector<Halfspace> globalConstraints;
    
    // Metadata registry
    std::vector<HPInfo> hpRegistry;
    
    // Core Logic
    void inserthp(long int hpid, int userIndex, int targetARF, cell* node, std::unordered_map<long int, bool> touchhs);
    bool isFeasible(std::unordered_map<long int, bool>& touchhs, long int hpid, bool side);
    void updateRank(cell* node, const int mink);
    void pruneRecursive(cell* node, int maxARF);
    
    // LP Helpers
    void addGlobalBoundaries(lprec* lp);
    void addHP(lprec* model, long int hpid, bool sideindicator);
    void lpModel(lprec* lp, int& dimen, std::unordered_map<long int, bool>& touchhs);

public:
    cellTree(Facet& facet, std::vector<Halfspace>& inequalities, int dimen);
    ~cellTree();
    
    void releaseCell(cell* node);
    void insertSlab(User* user, int userIndex, int slabIndex, int targetARF);
    void pruneByARF(int maxARF) { if (root) pruneRecursive(root, maxARF); }
    
    // Utilities
    Cell convertToCell(cell* node, std::vector<User*>& users);
    
    // Traversal
    void collectLeaves(std::vector<cell*>& leaves);
    void dfsTraversal(cell* node, cell* path, std::vector<cell*>& leaves);
    
    // Checks
    int getMinLeafARF();
    bool isLeafFullyHMinus(cell* leaf, int totalUsers);
    bool isUserHMinus(cell* leaf, int userIndex);
    int countHMinus(cell* leaf);
};

#endif
