#pragma once 

using namespace std;

#include <vector>

#include "halfspace.h"

class Slab{

    private: 
        const vector<double>* coords = nullptr;
        double upperScore = 0.0;
        double lowerScore = 0.0;
        int index = 0;

        bool sameCoordinates(const Slab& other) const;
    
    public:
        Slab() = default;
        ~Slab() = default;

        Slab(const vector<double>& coords, double upperScore, double lowerScore, int index);

        friend bool operator==(const Slab& a, const Slab& b) {
            return a.sameCoordinates(b) &&
                   a.upperScore == b.upperScore &&
                   a.lowerScore == b.lowerScore &&
                   a.index == b.index;
        }

        Slab& operator=(const Slab& other) = default;

        void rebindCoordinates(const vector<double>& coords);

        Halfspace getLower() const;
        Halfspace getUpper() const;
        double getLowerScore() const;
        double getUpperScore() const;
        int getIndex() const;
};
