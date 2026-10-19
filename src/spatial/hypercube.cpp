#include "hypercube.h"
#include <cassert>
#include <cstring>

// default: zero points (dim=0)
Hypercube::Hypercube() : m_lower(), m_upper() {}

Hypercube::Hypercube(const int a_dimen, const float* a_l, const float* a_u) {
    assert(a_dimen >= 0 && a_dimen <= MAXDIMEN);
    float lo[MAXDIMEN], hi[MAXDIMEN];
    for (int d = 0; d < a_dimen; ++d) {
        const float x = a_l ? a_l[d] : 0.0f;
        const float y = a_u ? a_u[d] : 0.0f;
        if (x <= y) { lo[d] = x; hi[d] = y; } else { lo[d] = y; hi[d] = x; }
    }
    for (int d = a_dimen; d < MAXDIMEN; ++d) { lo[d] = 0.0f; hi[d] = 0.0f; }
    m_lower = Point(a_dimen, lo);
    m_upper = Point(a_dimen, hi);
}

Hypercube::Hypercube(const Point& a_lower, const Point& a_upper) {
    assert(a_lower.m_dimen == a_upper.m_dimen);
    const int dim = a_lower.m_dimen;
    float lo[MAXDIMEN], hi[MAXDIMEN];
    for (int d = 0; d < dim; ++d) {
        const float x = a_lower.m_coor[d];
        const float y = a_upper.m_coor[d];
        if (x <= y) { lo[d] = x; hi[d] = y; } else { lo[d] = y; hi[d] = x; }
    }
    for (int d = dim; d < MAXDIMEN; ++d) { lo[d] = 0.0f; hi[d] = 0.0f; }
    m_lower = Point(dim, lo);
    m_upper = Point(dim, hi);
}

Hypercube::Hypercube(const Hypercube& a_hc) : m_lower(a_hc.m_lower), m_upper(a_hc.m_upper) {}

Hypercube::~Hypercube() {}

Hypercube& Hypercube::operator=(const Hypercube& a_hc) {
    if (this != &a_hc) { m_lower = a_hc.m_lower; m_upper = a_hc.m_upper; }
    return *this;
}

const Point& Hypercube::getLower() const { return m_lower; }
const Point& Hypercube::getUpper() const { return m_upper; }

Point Hypercube::getCenter() const {
    const int dim = dimen();
    float c[MAXDIMEN];
    for (int i = 0; i < dim; ++i) c[i] = (m_lower.m_coor[i] + m_upper.m_coor[i]) * 0.5f;
    for (int i = dim; i < MAXDIMEN; ++i) c[i] = 0.0f;
    return Point(dim, c);
}

int Hypercube::dimen() const {
    assert(m_lower.m_dimen == m_upper.m_dimen);
    return m_lower.m_dimen;
}

float Hypercube::volume() const {
    const int dim = dimen();
    if (dim == 0) return 0.0f;
    float vol = 1.0f;
    for (int d = 0; d < dim; ++d) vol *= (m_upper.m_coor[d] - m_lower.m_coor[d]);
    return vol;
}

bool Hypercube::operator==(const Hypercube& a_hc) const {
    return (m_lower == a_hc.m_lower) && (m_upper == a_hc.m_upper);
}

bool Hypercube::enclose(const Hypercube& a_hc) const {
    assert(dimen() == a_hc.dimen());
    const int dim = dimen();
    for (int d = 0; d < dim; ++d) {
        if (m_lower.m_coor[d] > a_hc.m_lower.m_coor[d]) return false;
        if (m_upper.m_coor[d] < a_hc.m_upper.m_coor[d]) return false;
    }
    return true;
}

bool Hypercube::enclose(const Point& a_pt) const {
    assert(dimen() == a_pt.m_dimen);
    const int dim = dimen();
    for (int d = 0; d < dim; ++d) {
        const float v = a_pt.m_coor[d];
        if (v < m_lower.m_coor[d] || v > m_upper.m_coor[d]) return false;
    }
    return true;
}

int Hypercube::size(const int a_dimen) {
    assert(a_dimen >= 0 && a_dimen <= MAXDIMEN);
    return (int)sizeof(float) * 2 * a_dimen;
}