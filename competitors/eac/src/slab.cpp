#include "slab.h"

Slab::Slab(Halfspace upper, Halfspace lower, int index){
    this->lower = lower;
    this->upper = upper;
    this->index = index;
}

Halfspace Slab::getLower(){
    return lower;
}

Halfspace Slab::getUpper(){
    return upper;
}

int Slab::getIndex(){
    return index;
}

Slab& Slab::operator=(const Slab& other) {
    if (this != &other) { // check for self-assignment
        upper = other.upper;
        lower = other.lower;
        index = other.index;
    }
    return *this;
}
