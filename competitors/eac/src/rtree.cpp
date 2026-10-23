#include "rtree.h"
#include "rentry.h"
#include "rnode.h"
#include "hypercube.h"
#include "point.h"
#include <cmath>
#include <cassert>

Rtree::Rtree(const int a_dimen, const int a_maxNodefill, const int a_maxLeafFill,
             const int a_minNodeFill, const int a_minLeafFill)
: m_dimen(a_dimen), m_root(nullptr),
  m_maxNodeChild(a_maxNodefill), m_minNodeChild(a_minNodeFill),
  m_maxLeafChild(a_maxLeafFill), m_minLeafChild(a_minLeafFill) { }

Rtree::~Rtree() { destroy(m_root); }

void Rtree::destroy(RtreeNode* n) {
    if (!n) return;
    for (int i = 0; i < n->m_count; ++i) {
        RtreeNodeEntry* e = n->m_entries[i];
        if (e && e->m_child) destroy(e->m_child);
    }
    delete n; // node dtor deletes entry wrappers
}

void Rtree::bulkLoad(RtreeNodeEntry** entries, int count, int maxLeafChild) {
    if (count <= 0) { m_root = nullptr; return; }
    m_root = bulkLoadSTR(entries, count, maxLeafChild, m_dimen, 0, 0);
    std::cout << "STR Bulkload: tree built with root node\n";
}

// Group children into parents of capacity m_maxNodeChild (no insert/split here)
RtreeNode* Rtree::buildParents(std::vector<RtreeNode*>& children) {
    if (children.empty()) return nullptr;
    if ((int)children.size() <= m_maxNodeChild) {
        int myLevel = children[0]->m_level + 1;
        RtreeNode* node = new RtreeNode(this, myLevel, nullptr);
        for (RtreeNode* ch : children) {
            ch->m_parent = node;
            node->appendEntry(new RtreeNodeEntry(ch));
        }
        return node;
    } else {
        std::vector<RtreeNode*> parents;
        parents.reserve((children.size() + (size_t)m_maxNodeChild - 1) / (size_t)m_maxNodeChild);
        int myLevel = children[0]->m_level + 1;
        for (size_t i = 0; i < children.size(); i += (size_t)m_maxNodeChild) {
            size_t end = std::min(children.size(), i + (size_t)m_maxNodeChild);
            RtreeNode* node = new RtreeNode(this, myLevel, nullptr);
            for (size_t j = i; j < end; ++j) {
                children[j]->m_parent = node;
                node->appendEntry(new RtreeNodeEntry(children[j]));
            }
            parents.push_back(node);
        }
        return buildParents(parents);
    }
}

// STR recursive slicing (classic S = ceil(L^(1/d)))
RtreeNode* Rtree::bulkLoadSTR(RtreeNodeEntry** entries, int count, int maxChild,
                              int dim, int curDim, int /*level*/)
{
    assert(dim > 0 && dim <= MAXDIMEN);

    if (count <= maxChild) {
        // Leaf
        RtreeNode* leaf = new RtreeNode(this, 0, nullptr);
        for (int i = 0; i < count; ++i) leaf->appendEntry(entries[i]);
        return leaf;
    }

    // leaves needed if this were the bottom
    int L = (count + maxChild - 1) / maxChild;
    int S = (int)std::ceil(std::pow((double)L, 1.0 / dim));
    if (S < 1) S = 1;

    std::sort(entries, entries + count, CompareByDim(curDim));

    int sliceSize = (count + S - 1) / S;
    std::vector<RtreeNode*> children;
    children.reserve(S);

    for (int s = 0; s < S; ++s) {
        int b = s * sliceSize;
        int e = std::min(count, b + sliceSize);
        if (b >= e) break;

        // recurse in next dimension
        RtreeNode* child = bulkLoadSTR(entries + b, e - b, maxChild,
                                       dim, (curDim + 1) % dim, 0);
        children.push_back(child);
    }

    // Build parent levels without overflow
    return buildParents(children);
}
