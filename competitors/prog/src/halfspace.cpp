#include "halfspace.h"
#include <cmath>

Halfspace::Halfspace(vector<double> coords, double score){
    this->coords = coords;
    this->score = score;
}

 void Halfspace::setScore(double score){
    this->score = score;
 }

vector<double>& Halfspace::getCoords(){
    return coords;
}

double Halfspace::getScore(){
    return score;
}

Halfspace& Halfspace::operator=(const Halfspace& other) {
    if (this != &other) {
        coords = other.coords;
        score = other.score;
    }
    return *this;
}
