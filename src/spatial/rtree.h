#pragma once
#include <vector>
#include <algorithm>
#include "rentry.h"

struct RtreeNode;
class RtreeNodeEntry;

struct CompareByDim {
    int dim;
    CompareByDim(int d) : dim(d) {}
    bool operator()(RtreeNodeEntry* a, RtreeNodeEntry* b) const {
        // use raw array; no Point::operator[]
        return a->m_mbr.getLower().m_coor[dim] < b->m_mbr.getLower().m_coor[dim];
    }
};

class Rtree {
public:
    int m_dimen;
    RtreeNode* m_root;
    const int m_maxNodeChild;
    const int m_minNodeChild;
    const int m_maxLeafChild;
    const int m_minLeafChild;

    Rtree(const int a_dimen, const int a_maxNodefill, const int a_maxLeafFill,
          const int a_minNodeFill, const int a_minLeafFill);
    virtual ~Rtree();

    void bulkLoad(RtreeNodeEntry** entries, int count, int maxLeafChild);

    RtreeNode* bulkLoadSTR(RtreeNodeEntry** entries, int count, int maxChild,
                           int dim, int curDim, int level = 0);

    void destroy(RtreeNode* n);

private:
    RtreeNode* buildParents(std::vector<RtreeNode*>& children); 
};
