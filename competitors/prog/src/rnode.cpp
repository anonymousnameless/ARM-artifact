#include "rnode.h"
#include "rentry.h"
#include "rtree.h"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <vector>
#include <iostream>

RtreeNode::RtreeNode(Rtree* a_rtree, const int a_level, RtreeNode* a_parent)
: m_isLeaf(a_level == 0),
  m_rtree(a_rtree),
  m_count(0),
  m_level(a_level),
  m_parent(a_parent),
  m_usedspace(0),
  m_entries(nullptr)
{
    int capacity = m_isLeaf ? m_rtree->m_maxLeafChild : m_rtree->m_maxNodeChild;
    m_entries = new RtreeNodeEntry*[capacity];
    std::memset(m_entries, 0, sizeof(RtreeNodeEntry*) * capacity);
}

RtreeNode::~RtreeNode() {
    for (int i = 0; i < m_count; ++i) delete m_entries[i];
    delete[] m_entries;
}

int RtreeNode::appendEntry(RtreeNodeEntry* entry) {
    int capacity = (m_level == 0 ? m_rtree->m_maxLeafChild : m_rtree->m_maxNodeChild);
    if (m_count >= capacity) return 1;
    m_entries[m_count++] = entry;
    if (entry->m_child) entry->m_child->m_parent = this;
    return 0;
}

int RtreeNode::insertEntry(RtreeNodeEntry* entry)
{
    const int maxEntries = (m_level == 0 ? m_rtree->m_maxLeafChild : m_rtree->m_maxNodeChild);
    if (m_count < maxEntries) {
        m_entries[m_count++] = entry;
        if (entry->m_child) entry->m_child->m_parent = this;
        return 0;
    }

    // Split (quadratic), with capacity enforcement
    RtreeNode* node1 = new RtreeNode(m_rtree, m_level, m_parent);
    RtreeNode* node2 = new RtreeNode(m_rtree, m_level, m_parent);

    std::vector<RtreeNodeEntry*> allEntries(m_entries, m_entries + m_count);
    allEntries.push_back(entry);

    int seed1 = 0, seed2 = 1;
    float maxD = -1.0f;
    for (int i = 0; i < (int)allEntries.size(); ++i) {
        for (int j = i + 1; j < (int)allEntries.size(); ++j) {
            Hypercube comb = combine(allEntries[i]->m_mbr, allEntries[j]->m_mbr);
            float d = comb.volume() - allEntries[i]->m_mbr.volume() - allEntries[j]->m_mbr.volume();
            if (d > maxD) { maxD = d; seed1 = i; seed2 = j; }
        }
    }

    node1->appendEntry(allEntries[seed1]);
    node2->appendEntry(allEntries[seed2]);

    for (int i = 0; i < (int)allEntries.size(); ++i) {
        if (i == seed1 || i == seed2) continue;

        // Capacity enforcement:
        if (node1->m_count >= maxEntries) { node2->appendEntry(allEntries[i]); continue; }
        if (node2->m_count >= maxEntries) { node1->appendEntry(allEntries[i]); continue; }

        Hypercube mbr1 = combine(node1->computeMBR(), allEntries[i]->m_mbr);
        Hypercube mbr2 = combine(node2->computeMBR(), allEntries[i]->m_mbr);
        float increase1 = mbr1.volume() - node1->computeMBR().volume();
        float increase2 = mbr2.volume() - node2->computeMBR().volume();
        if (increase1 < increase2) node1->appendEntry(allEntries[i]);
        else                       node2->appendEntry(allEntries[i]);
    }

    if (m_parent) {
        m_parent->replaceChild(this, node1, node2);
    } else {
        RtreeNode* newRoot = new RtreeNode(m_rtree, m_level + 1, nullptr);
        node1->m_parent = newRoot;
        node2->m_parent = newRoot;
        newRoot->appendEntry(new RtreeNodeEntry(node1));
        newRoot->appendEntry(new RtreeNodeEntry(node2));
        if (m_rtree) m_rtree->m_root = newRoot;
    }

    // Old node entries are moved; detach to avoid double delete
    m_count = 0;
    return 1;
}

void RtreeNode::replaceChild(RtreeNode* oldChild, RtreeNode* newChild1, RtreeNode* newChild2)
{
    int idx = -1;
    for (int i = 0; i < m_count; ++i) {
        if (m_entries[i]->m_child == oldChild) { idx = i; break; }
    }
    if (idx == -1) return;

    delete m_entries[idx];
    delete oldChild;

    for (int j = idx; j < m_count - 1; ++j) m_entries[j] = m_entries[j + 1];
    --m_count;

    newChild1->m_parent = this;
    newChild2->m_parent = this;
    appendEntry(new RtreeNodeEntry(newChild1));
    appendEntry(new RtreeNodeEntry(newChild2));

    const int maxEntries = (m_level == 0 ? m_rtree->m_maxLeafChild : m_rtree->m_maxNodeChild);
    if (m_count > maxEntries) {
        if (m_parent) {
            RtreeNode *n1 = nullptr, *n2 = nullptr;
            splitNode(n1, n2);
            m_parent->replaceChild(this, n1, n2);
        } else {
            RtreeNode *n1 = nullptr, *n2 = nullptr;
            splitNode(n1, n2);
            RtreeNode* newRoot = new RtreeNode(m_rtree, m_level + 1, nullptr);
            n1->m_parent = newRoot; n2->m_parent = newRoot;
            newRoot->appendEntry(new RtreeNodeEntry(n1));
            newRoot->appendEntry(new RtreeNodeEntry(n2));
            if (m_rtree) m_rtree->m_root = newRoot;
        }
    }
}

Hypercube RtreeNode::computeMBR() const
{
    const int dim = m_rtree->m_dimen;

    if (m_count == 0) {
        float lo[MAXDIMEN] = {0};
        float hi[MAXDIMEN] = {0};
        return Hypercube(dim, lo, hi);
    }

    float lo[MAXDIMEN];
    float hi[MAXDIMEN];

    const Point& L0 = m_entries[0]->m_mbr.getLower();
    const Point& U0 = m_entries[0]->m_mbr.getUpper();
    for (int d = 0; d < dim; ++d) { lo[d] = L0.m_coor[d]; hi[d] = U0.m_coor[d]; }

    for (int i = 1; i < m_count; ++i) {
        const Point& L = m_entries[i]->m_mbr.getLower();
        const Point& U = m_entries[i]->m_mbr.getUpper();
        for (int d = 0; d < dim; ++d) {
            if (L.m_coor[d] < lo[d]) lo[d] = L.m_coor[d];
            if (U.m_coor[d] > hi[d]) hi[d] = U.m_coor[d];
        }
    }
    return Hypercube(dim, lo, hi);
}

Hypercube RtreeNode::combine(const Hypercube& a, const Hypercube& b)
{
    const int dim = a.dimen();
    float lo[MAXDIMEN], hi[MAXDIMEN];

    const Point& La = a.getLower(); const Point& Ua = a.getUpper();
    const Point& Lb = b.getLower(); const Point& Ub = b.getUpper();

    for (int d = 0; d < dim; ++d) {
        lo[d] = (La.m_coor[d] < Lb.m_coor[d]) ? La.m_coor[d] : Lb.m_coor[d];
        hi[d] = (Ua.m_coor[d] > Ub.m_coor[d]) ? Ua.m_coor[d] : Ub.m_coor[d];
    }
    return Hypercube(dim, lo, hi);
}

void RtreeNode::splitNode(RtreeNode*& node1, RtreeNode*& node2)
{
    const int maxEntries = (m_level == 0 ? m_rtree->m_maxLeafChild : m_rtree->m_maxNodeChild);

    node1 = new RtreeNode(m_rtree, m_level, m_parent);
    node2 = new RtreeNode(m_rtree, m_level, m_parent);

    int seed1 = 0, seed2 = 1;
    float maxD = -1.0f;
    for (int i = 0; i < m_count; ++i) {
        for (int j = i + 1; j < m_count; ++j) {
            Hypercube comb = combine(m_entries[i]->m_mbr, m_entries[j]->m_mbr);
            float d = comb.volume() - m_entries[i]->m_mbr.volume() - m_entries[j]->m_mbr.volume();
            if (d > maxD) { maxD = d; seed1 = i; seed2 = j; }
        }
    }

    node1->appendEntry(m_entries[seed1]);
    node2->appendEntry(m_entries[seed2]);

    for (int i = 0; i < m_count; ++i) {
        if (i == seed1 || i == seed2) continue;

        // capacity enforcement
        if (node1->m_count >= maxEntries) { node2->appendEntry(m_entries[i]); continue; }
        if (node2->m_count >= maxEntries) { node1->appendEntry(m_entries[i]); continue; }

        Hypercube mbr1 = combine(node1->computeMBR(), m_entries[i]->m_mbr);
        Hypercube mbr2 = combine(node2->computeMBR(), m_entries[i]->m_mbr);
        float increase1 = mbr1.volume() - node1->computeMBR().volume();
        float increase2 = mbr2.volume() - node2->computeMBR().volume();

        if (increase1 < increase2) node1->appendEntry(m_entries[i]);
        else                       node2->appendEntry(m_entries[i]);
    }

    m_count = 0;
}
