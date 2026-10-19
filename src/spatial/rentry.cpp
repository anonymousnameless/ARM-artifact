#include "rentry.h"
#include "rnode.h"
#include <cassert>

RtreeNodeEntry::RtreeNodeEntry() : m_mbr(), m_id(-1), m_child(nullptr) {}

RtreeNodeEntry::RtreeNodeEntry(int id, const Hypercube& hc) : m_mbr(hc), m_id(id), m_child(nullptr) {}

RtreeNodeEntry::RtreeNodeEntry(RtreeNode* child) : m_mbr(), m_id(-1), m_child(child) {
    assert(child != nullptr);
    m_mbr = child->computeMBR();
}

bool RtreeNodeEntry::isLeaf() const { return m_child == nullptr; }

void RtreeNodeEntry::setLeaf(int id, const Hypercube& hc) {
    m_id = id; 
    m_mbr = hc; 
    m_child = nullptr;
}

void RtreeNodeEntry::setChild(RtreeNode* child) {
    assert(child != nullptr);
    m_child = child; m_id = -1; m_mbr = child->computeMBR();
}