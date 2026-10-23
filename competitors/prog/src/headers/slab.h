#pragma once

using namespace std;

#include <vector>

#include "halfspace.h"

class Slab{

    private:
        Halfspace upper;
        Halfspace lower;
        int index;

    public:
        Slab() = default;
        ~Slab() = default;

        Slab( Halfspace upper, Halfspace lower, int index);

        friend bool operator==(const Slab& a, const Slab& b) {
            return a.upper == b.upper && a.upper == b.upper && a.index == b.index;
        }

        Slab& operator=(const Slab& other);

        Halfspace getLower();
        Halfspace getUpper();
        int getIndex();
};
