#pragma once
#include "header.h"  // must define MAXDIMEN

class Point {
public:
    int   m_dimen;              // 0..MAXDIMEN
    float m_coor[MAXDIMEN];

public:
    Point();
    Point(const int a_dimen, const float* a_coor = 0);
    Point(const Point& a_pt);
    ~Point();
    Point& operator=(const Point& a_pt);

    bool operator==(const Point& o) const;
    bool operator!=(const Point& o) const;
};
