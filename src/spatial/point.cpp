#include "point.h"
#include <cstring>
#include <cassert>

Point::Point() : m_dimen(0) {
    std::memset(m_coor, 0, sizeof(m_coor));
}

Point::Point(const int a_dimen, const float* a_coor) : m_dimen(a_dimen) {
    assert(a_dimen >= 0 && a_dimen <= MAXDIMEN);
    if (a_coor != 0) {
        std::memcpy(m_coor, a_coor, sizeof(float) * m_dimen);
    } else {
        std::memset(m_coor, 0, sizeof(float) * m_dimen);
    }
    if (m_dimen < MAXDIMEN)
        std::memset(m_coor + m_dimen, 0, sizeof(float) * (MAXDIMEN - m_dimen));
}

Point::Point(const Point& a_pt) : m_dimen(a_pt.m_dimen) {
    std::memcpy(m_coor, a_pt.m_coor, sizeof(float) * m_dimen);
    if (m_dimen < MAXDIMEN)
        std::memset(m_coor + m_dimen, 0, sizeof(float) * (MAXDIMEN - m_dimen));
}

Point::~Point() {}

Point& Point::operator=(const Point& a_pt) {
    if (this != &a_pt) {
        m_dimen = a_pt.m_dimen;
        std::memcpy(m_coor, a_pt.m_coor, sizeof(float) * m_dimen);
        if (m_dimen < MAXDIMEN)
            std::memset(m_coor + m_dimen, 0, sizeof(float) * (MAXDIMEN - m_dimen));
    }
    return *this;
}

bool Point::operator==(const Point& o) const {
    if (m_dimen != o.m_dimen) return false;
    for (int i = 0; i < m_dimen; ++i)
        if (m_coor[i] != o.m_coor[i]) return false;
    return true;
}

bool Point::operator!=(const Point& o) const {
    return !(*this == o);
}
