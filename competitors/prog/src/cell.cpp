#include "cell.h"

Cell::Cell(const vector<int>& combination){
    updateCell(combination);
}

void Cell::updateCell(const vector<int>& combination){
    combinationIndex = combination;
    setCombinationIndex();
}

void Cell::setCombinationIndex(){
   arf = computeARF(combinationIndex);
}

bool Cell::weaklyDominates(const Cell& other) const {
    bool dominates = true;
    for(int i=0;i<this->combinationIndex.size();i++){
        if(other.getCombinationIndex()[i] < this->combinationIndex[i])
            dominates = false;
    }

    if (dominates) return true;
    return false;
}

double Cell::getARF(){
    return arf;
}

const vector<int>& Cell::getCombinationIndex() const {
    return combinationIndex;
}

const vector<int>& Cell::getCombinationIndexOnly() const {
    return combinationIndex;
}

void Cell::clean(){
    arf = -1;
    combinationIndex.clear();
}

double Cell::computeARF(const vector<int>& indices) const {
    double arf = 0;
    for (int idx : indices) {
        arf += idx;
    }
    return arf;
}
