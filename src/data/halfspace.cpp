#include "halfspace.h"
#include <cmath>

Halfspace::Halfspace(const vector<double>& coords, double score)
    : coords(coords), score(score) {}

 void Halfspace::setScore(double score){
    this->score = score;
 }

vector<double>& Halfspace::getCoords(){
    return coords;
}

const vector<double>& Halfspace::getCoords() const {
    return coords;
}

double Halfspace::getScore(){
    return score;
}

double Halfspace::getScore() const {
    return score;
}

Halfspace& Halfspace::operator=(const Halfspace& other) {
    if (this != &other) {
        coords = other.coords;
        score = other.score;
    }
    return *this;
}
