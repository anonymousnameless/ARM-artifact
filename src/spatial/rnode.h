#pragma once

#include "../header.h"
#include "hypercube.h"

struct Rtree;            // forward
class RtreeNodeEntry;    // forward

class RtreeNode {
public:
    bool             m_isLeaf;
    Rtree*           m_rtree;      // owner tree
    int              m_count;      // used entries
    long long        cachedProductCount = -1; // populated after bulk loading, before read-only pruning
    int              m_level;      // 0 = leaf
    RtreeNode*       m_parent;
    int              m_usedspace;  // currently unused
    RtreeNodeEntry** m_entries;    // pointer array (capacity by level)

public:
    RtreeNode(Rtree* a_rtree, int a_level, RtreeNode* a_parent);
    ~RtreeNode();

    int        appendEntry(RtreeNodeEntry* entry);   // add if capacity, else overflow=1
    int        insertEntry(RtreeNodeEntry* entry);   // add or split; returns 0 if placed, 1 if split happened
    Hypercube  computeMBR() const;

    Hypercube  combine(const Hypercube& a, const Hypercube& b);
    void       replaceChild(RtreeNode* oldChild, RtreeNode* newChild1, RtreeNode* newChild2);
    void       splitNode(RtreeNode*& node1, RtreeNode*& node2);

    static int size() { return (int)sizeof(bool) + (int)sizeof(int) + (int)sizeof(RtreeNodeEntry*) * MAX_CHILDREN; }
};
