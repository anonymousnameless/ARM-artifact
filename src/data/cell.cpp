#include "cell.h"
#include "slab.h"
#include "user.h"
#include <iostream>

Cell::Cell(vector<pair<Slab,User*>> combination){
    for(int i=0;i<combination.size();i++){
        this->combination.push_back(combination[i]);
    }
    setCombinationIndex();
}

Cell::Cell(vector<int>& combination,vector<User*>& users){
    updateCell(combination,users);
}

void Cell::updateCell(vector<int>& combination, vector<User*>& users){
    this->combination.clear();
    // Guard: combination must match users size to avoid out-of-bounds access
    int safeSize = min((int)combination.size(), (int)users.size());
    for(int i=0;i<safeSize;i++){
        this->combination.push_back({users[i]->getSlab(combination[i]),users[i]});
    }  
    setCombinationIndex();
}

void Cell::setCombinationIndex(){
    // Rebuild the indices rather than retaining entries from a previous cell.
    combinationIndex.clear();
    arf = 0;
    vector<int> arfIndex;
    for(int i=0;i<combination.size();i++){
        combinationIndex.push_back({combination[i].first.getIndex(),combination[i].second});
        arfIndex.push_back(combination[i].first.getIndex());
    }
   arf = computeARF(arfIndex);
}

bool Cell::weaklyDominates(Cell& other){ 
    bool dominates = true; 
    for(int i=0;i<this->combinationIndex.size();i++){
        if(other.getCombinationIndex()[i].first < this->combinationIndex[i].first)
            dominates = false;
    }

    if (dominates) return true;
    return false;
}

double Cell::getARF() const{
    return arf;
}

vector<pair<Slab,User*>>& Cell::getCombination(){
    return combination;
}

vector<pair<int,User*>>& Cell::getCombinationIndex(){
    return combinationIndex;
}

void Cell::clean(){
    combination.clear();
    arf = -1;
    combinationIndex.clear();
}

double Cell::computeARF(vector<int>& indices) {
    double arf = 0;
    for (int idx : indices) {
        arf += idx;
    }
    return arf;
}
