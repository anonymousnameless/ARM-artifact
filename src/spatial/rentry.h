#pragma once

struct RtreeNode;      // forward
#include "hypercube.h"

class RtreeNodeEntry {
public:
    Hypercube  m_mbr;      // Minimum bounding rectangle
    int        m_id;       // object id (for leaf nodes); -1 for internal
    RtreeNode* m_child;    // child node (nullptr if leaf)

public:
    RtreeNodeEntry();                                  // default (leaf, id=-1, empty MBR)
    RtreeNodeEntry(int id, const Hypercube& hc);       // leaf constructor
    RtreeNodeEntry(RtreeNode* child);                  // internal constructor

    bool isLeaf() const;
    void setLeaf(int id, const Hypercube& hc);
    void setChild(RtreeNode* child);

    static int size(int dimen) { return Hypercube::size(dimen) + (int)sizeof(int) + (int)sizeof(RtreeNode*); }
};