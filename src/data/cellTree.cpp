#include "cellTree.h"
#include "../data/user.h"
#include "../data/cell.h"
#include "../data/facet.h"
#include "../data/slab.h"
#include "../data/halfspace.h"
#include <climits>

cellTree::cellTree(Facet& facet, std::vector<Halfspace>& inequalities, int dimen)
    : dimen(dimen), rootFacet(&facet), globalConstraints(inequalities) {
    root = new cell();
}

cellTree::~cellTree() {
    releaseCell(root);
}

void cellTree::releaseCell(cell* node) {
    if (!node) return;
    releaseCell(node->left);
    releaseCell(node->right);
    node->release();
    delete node;
}

// ---------------------------------------------------------
// INSERTION
// ---------------------------------------------------------

void cellTree::insertSlab(User* user, int userIndex, int slabIndex, int targetARF) {
    Halfspace hp = user->getSlab(slabIndex).getLower();
    std::vector<float> coeffs(dimen + 2);
    for (int d = 0; d < dimen; ++d) coeffs[d] = hp.getCoords()[d];
    coeffs[dimen] = hp.getScore();
    coeffs[dimen + 1] = 0.0f;
    
    HalfSpaces.push_back(coeffs);
    long int hpID = HalfSpaces.size();
    
    hpRegistry.push_back({user, userIndex, slabIndex});
    
    std::unordered_map<long int, bool> emptyContext;
    inserthp(hpID, userIndex, targetARF, root, emptyContext);
}

void cellTree::inserthp(long int hpid, int userIndex, int targetARF, cell* node, std::unordered_map<long int, bool> touchhs) {
    
    // 1. PRUNING CHECK
    if (node->isPruned || node->arf >= targetARF) {
        node->isPruned = true;
        return;
    }

    // 2. Build Context
    for (auto& kv : node->IntersectHP) {
        touchhs[kv.first] = kv.second;
    }

    // 3. OPTIMIZATION: Check History for "H- (Right)"
    // If we are already Right of a tighter slab for this user, we are automatically Right of this one.
    for (auto& kv : touchhs) {
        long int ctxID = kv.first;
        bool side = kv.second;
        if (side == false) { // Right/H-
            HPInfo& ctxInfo = hpRegistry[ctxID - 1];
            if (ctxInfo.userIndex == userIndex) {
                node->BelowHP.insert(hpid);
                return;
            }
        }
    }
    
    for (long int ctxID : node->BelowHP) {
        HPInfo& ctxInfo = hpRegistry[ctxID - 1];
        if (ctxInfo.userIndex == userIndex) {
            node->BelowHP.insert(hpid);
            return;
        }
    }

    // 4. LP CHECKS
    bool feasLeft  = isFeasible(touchhs, hpid, true);  // Can we go Left?
    bool feasRight = isFeasible(touchhs, hpid, false); // Can we go Right?

    if (!feasLeft && !feasRight) {
        node->isPruned = true; // Impossible geometry
        return;
    }

    if (!feasLeft) {
        node->BelowHP.insert(hpid); // Forced Right (H-)
    }
    else if (!feasRight) {
        node->AboveHP.insert(hpid); // Forced Left (H+)
        updateRank(node, targetARF); 
    }
    else {
        // Split
        if (node->left == nullptr) {
            // NEW LEFT (H+): ARF Increases (user didn't achieve this slab)
            node->left = new cell(node);
            node->left->IntersectHP[hpid] = true;
            node->left->arf = node->arf + 1; 
            if (node->left->arf >= targetARF) node->left->isPruned = true;

            // NEW RIGHT (H-): ARF Same (user achieved this slab)
            node->right = new cell(node);
            node->right->IntersectHP[hpid] = false;
            node->right->arf = node->arf;
        }
        else {
            // RECURSE DOWN EXISTING CHILDREN
            auto leftContext = touchhs;
            inserthp(hpid, userIndex, targetARF, node->left, leftContext);
            
            auto rightContext = touchhs;
            inserthp(hpid, userIndex, targetARF, node->right, rightContext);
            
            // Garbage Collection
            if ((!node->left || node->left->isPruned) && 
                (!node->right || node->right->isPruned)) {
                
                node->isPruned = true;
                if (node->left) { releaseCell(node->left); node->left = nullptr; }
                if (node->right) { releaseCell(node->right); node->right = nullptr; }
            }
        }
    }
}

void cellTree::updateRank(cell* node, const int mink) {
    node->arf++; 
    
    if (node->arf >= mink) {
        node->isPruned = true;
        releaseCell(node->left); node->left = nullptr;
        releaseCell(node->right); node->right = nullptr;
    } else {
        if (node->left) updateRank(node->left, mink);
        if (node->right) updateRank(node->right, mink);
    }
}

void cellTree::pruneRecursive(cell* node, int maxARF) {
    if (!node || node->isPruned) return;
    
    if (node->arf >= maxARF) {
        node->isPruned = true;
        releaseCell(node->left); node->left = nullptr;
        releaseCell(node->right); node->right = nullptr;
        return;
    }
    
    if (node->left) pruneRecursive(node->left, maxARF);
    if (node->right) pruneRecursive(node->right, maxARF);
    
    if ((!node->left || node->left->isPruned) && (!node->right || node->right->isPruned) && (node->left != nullptr)) {
        node->isPruned = true;
        releaseCell(node->left); node->left = nullptr;
        releaseCell(node->right); node->right = nullptr;
    }
}

// ---------------------------------------------------------
// LP / JUDGE LOGIC
// ---------------------------------------------------------

void cellTree::addGlobalBoundaries(lprec* lp) {
    for (auto& ineq : globalConstraints) {
        std::vector<double> row(dimen + 1, 0.0);
        for (int i = 0; i < dimen; ++i) row[i + 1] = ineq.getCoords()[i];
        add_constraint(lp, row.data(), LE, ineq.getScore());
    }

    FacetVRep& vrep = rootFacet->getFacetVRep();
    std::vector<double> fRow(dimen + 1, 0.0);
    for (int i = 0; i < dimen; ++i) fRow[i + 1] = vrep.hyperplane_normal[i];
    add_constraint(lp, fRow.data(), EQ, vrep.hyperplane_offset);
}

void cellTree::addHP(lprec* model, long int hpid, bool sideindicator) {
    int index = hpid - 1; 
    std::vector<double> row(dimen + 1, 0.0);
    for (int i = 1; i <= dimen; i++) row[i] = HalfSpaces[index][i - 1];
    double rhs = HalfSpaces[index][dimen];

    // true=Left(LE), false=Right(GE)
    add_constraint(model, row.data(), sideindicator ? LE : GE, rhs);
}

void cellTree::lpModel(lprec* lp, int& dimen, std::unordered_map<long int, bool>& touchhs) {
    addGlobalBoundaries(lp);
    for (auto& kv : touchhs) {
        addHP(lp, kv.first, kv.second);
    }
    set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
}

bool cellTree::isFeasible(std::unordered_map<long int, bool>& touchhs, long int hpid, bool sideindicator) {
    lprec *lp = make_lp(0, dimen);
    set_verbose(lp, 1); 

    lpModel(lp, dimen, touchhs);
    addHP(lp, hpid, sideindicator);

    std::vector<double> obj(dimen + 1, 0.0);
    set_obj_fn(lp, obj.data());

    int ret = solve(lp);
    delete_lp(lp);
    
    // CORRECTED: Return true only if a solution was FOUND.
    // 0 = OPTIMAL, 1 = SUBOPTIMAL. (2 = INFEASIBLE).
    return (ret == 0 || ret == 1);
}

// ---------------------------------------------------------
// UTILITIES
// ---------------------------------------------------------

Cell cellTree::convertToCell(cell* node, std::vector<User*>& all_users) {
    std::map<User*, int> userSlabMap;

    // 1. Initialize with the Facet's geometric minimums
    std::vector<int> minBorders = rootFacet->getMinBorders();
    for (size_t i = 0; i < all_users.size(); ++i) {
        if (i < minBorders.size()) {
            userSlabMap[all_users[i]] = minBorders[i];
        } else {
            userSlabMap[all_users[i]] = 0;
        }
    }

    // Helper lambda to update the map safely
    auto updateRankIfWorse = [&](long int hpID) {
        HPInfo& info = hpRegistry[hpID - 1];
        
        // If we are on the Left (H+) of Slab K, it means we FAILED Slab K.
        // Therefore, our rank is at least K + 1.
        int impliedRank = info.slabIndex + 1;
        
        if (impliedRank > userSlabMap[info.user]) {
            userSlabMap[info.user] = impliedRank;
        }
    };

    // 2. Check AboveHP (Fixed Left/H+)
    for (long int id : node->AboveHP) {
        updateRankIfWorse(id);
    }

    // 3. Check IntersectHP (Path Left/H+)
    for (auto& kv : node->IntersectHP) {
        // kv.second == true means Left/H+ (Satisfied <=, meaning Worse Score/Rank)
        if (kv.second == true) { 
            updateRankIfWorse(kv.first);
        }
    }

    // 4. Calculate and store realarf (sum of actual slab indices)
    int computedRealARF = 0;
    for (const auto& [user, slabIdx] : userSlabMap) {
        computedRealARF += slabIdx;
    }
    node->realarf = computedRealARF;

    // 5. Construct the Cell
    std::vector<std::pair<Slab, User*>> combination;
    for (User* u : all_users) {
        int bestSlabIdx = userSlabMap[u];
        
        // Safety cap (can't go beyond last slab)
        if (bestSlabIdx >= u->getSlabSize()) {
            bestSlabIdx = u->getSlabSize() - 1;
        }
        
        combination.push_back({u->getSlab(bestSlabIdx), u});
    }
    
    return Cell(combination);
}

void cellTree::collectLeaves(std::vector<cell*>& leaves) {
    leaves.clear();
    if (root == nullptr || root->isPruned) return;

    // Start with a new accumulator for the root.
    cell* initialPath = new cell();
    root->copyleaf(initialPath);

    dfsTraversal(root, initialPath, leaves);
    
    std::sort(leaves.begin(), leaves.end(), [](cell* a, cell* b) { return a->realarf < b->realarf; });
}

void cellTree::dfsTraversal(cell* node, cell* currentPath, std::vector<cell*>& leaves) {
    // Update path with current node's info (MERGE)
    node->copyleaf(currentPath);
    
    // Leaf Check
    if (node->left == nullptr && node->right == nullptr) {
        // Calculate realarf for the leaf path
        std::vector<int> currentRanks = rootFacet->getMinBorders();
        
        for (long int id : currentPath->AboveHP) {
            if (id - 1 >= 0 && id - 1 < hpRegistry.size()) {
                HPInfo& info = hpRegistry[id - 1];
                if (info.userIndex < currentRanks.size()) {
                    if (info.slabIndex + 1 > currentRanks[info.userIndex]) {
                        currentRanks[info.userIndex] = info.slabIndex + 1;
                    }
                }
            }
        }
        for (auto& kv : currentPath->IntersectHP) {
            // kv.second == true means Left/H+ (Failed/Worse)
            if (kv.second == true) { 
                if (kv.first - 1 >= 0 && kv.first - 1 < hpRegistry.size()) {
                    HPInfo& info = hpRegistry[kv.first - 1];
                    if (info.userIndex < currentRanks.size()) {
                        if (info.slabIndex + 1 > currentRanks[info.userIndex]) {
                            currentRanks[info.userIndex] = info.slabIndex + 1;
                        }
                    }
                }
            }
        }

        int calcARF = 0;
        for (int r : currentRanks) calcARF += r;
        currentPath->realarf = calcARF;

        leaves.push_back(currentPath);
        return;
    }

    bool wentLeft = false;
    
    if (node->left != nullptr && !node->left->isPruned) {
        cell* leftPath = new cell();
        currentPath->copyleaf(leftPath); // Clone current history
        dfsTraversal(node->left, leftPath, leaves);
        wentLeft = true;
    }
    
    if (node->right != nullptr && !node->right->isPruned) {
        dfsTraversal(node->right, currentPath, leaves);
    } 
    else {
        // Did not go right. If we went left, 'leftPath' is the clone.
        // 'currentPath' is now unused and must be deleted.
        delete currentPath;
    }
}

int cellTree::getMinLeafARF() {
    std::vector<cell*> leaves;
    collectLeaves(leaves);
    if (leaves.empty()) return INT_MAX;
    int minARF = leaves[0]->realarf;
    for (auto* l : leaves) releaseCell(l); // Clean up temp leaves
    return minARF;
}

bool cellTree::isLeafFullyHMinus(cell* leaf, int totalUsers) {
    std::set<int> hMinusUsers;
    
    for (size_t i = 0; i < hpRegistry.size(); ++i) {
        long int hpID = i + 1;
        HPInfo& info = hpRegistry[i];

        bool isHMinus = false;

        // Check Fixed Right
        if (leaf->BelowHP.count(hpID)) {
            isHMinus = true;
        } 
        // Check Path Right
        else if (leaf->IntersectHP.count(hpID)) {
            if (leaf->IntersectHP.at(hpID) == false) {
                isHMinus = true;
            }
        }

        if (isHMinus) {
            hMinusUsers.insert(info.userIndex);
        }
    }

    return (int)hMinusUsers.size() == totalUsers;
}

bool cellTree::isUserHMinus(cell* leaf, int userIndex) {
    for (size_t i = 0; i < hpRegistry.size(); ++i) {
        HPInfo& info = hpRegistry[i];
        if (info.userIndex != userIndex) continue;

        long int hpID = i + 1;
        // Check Fixed Right
        if (leaf->BelowHP.count(hpID)) return true;
        // Check Path Right
        if (leaf->IntersectHP.count(hpID)) {
             if (leaf->IntersectHP.at(hpID) == false) return true;
        }
    }
    return false; 
}

int cellTree::countHMinus(cell* leaf) {
    std::set<int> hMinusUsers;
    for (size_t i = 0; i < hpRegistry.size(); ++i) {
        long int hpID = i + 1;
        HPInfo& info = hpRegistry[i];

        bool isHMinus = false;
        if (leaf->BelowHP.count(hpID)) isHMinus = true;
        else if (leaf->IntersectHP.count(hpID)) {
            if (leaf->IntersectHP.at(hpID) == false) isHMinus = true;
        }

        if (isHMinus) hMinusUsers.insert(info.userIndex);
    }
    return (int)hMinusUsers.size();
}

