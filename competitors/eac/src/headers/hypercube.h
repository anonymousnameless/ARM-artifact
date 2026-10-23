#pragma once
#include <cmath>
#include "point.h"

class Hypercube {
protected:
    Point m_lower;   // lower corner
    Point m_upper;   // upper corner
public:
    Hypercube();
    Hypercube(const int a_dimen, const float* a_l, const float* a_u);
    Hypercube(const Point& a_lower, const Point& a_upper);
    Hypercube(const Hypercube& a_hc);
    virtual ~Hypercube();

    Hypercube& operator=(const Hypercube& a_hc);

    const Point& getLower() const;
    const Point& getUpper() const;
    Point        getCenter() const;
    int          dimen() const;
    float        volume() const;

    bool operator==(const Hypercube& a_hc) const;
    bool enclose(const Hypercube& a_hc) const;
    bool enclose(const Point& a_pt) const;

    static int size(const int a_dimen);
};
