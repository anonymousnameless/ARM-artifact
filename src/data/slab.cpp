#include "slab.h"

bool Slab::sameCoordinates(const Slab& other) const {
    if (coords == other.coords) {
        return true;
    }
    if (coords == nullptr || other.coords == nullptr) {
        return coords == other.coords;
    }
    return *coords == *other.coords;
}

Slab::Slab(const vector<double>& coords, double upperScore, double lowerScore, int index)
    : coords(&coords), upperScore(upperScore), lowerScore(lowerScore), index(index) {}

void Slab::rebindCoordinates(const vector<double>& coords) {
    this->coords = &coords;
}

Halfspace Slab::getLower() const {
    static const vector<double> emptyCoords;
    return Halfspace(coords ? *coords : emptyCoords, lowerScore);
}

Halfspace Slab::getUpper() const {
    static const vector<double> emptyCoords;
    return Halfspace(coords ? *coords : emptyCoords, upperScore);
}

double Slab::getLowerScore() const {
    return lowerScore;
}

double Slab::getUpperScore() const {
    return upperScore;
}

int Slab::getIndex() const {
    return index;
}
