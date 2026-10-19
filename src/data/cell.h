#pragma once 
#include <vector>

using namespace std;

#include "slab.h"

struct User;

class Cell
{
     private:
        vector<pair<Slab,User*>> combination;
        double arf;
        vector<pair<int,User*>> combinationIndex;

    public:
        Cell() = default;
        ~Cell() = default;

        Cell(vector<pair<Slab,User*>> combination);
        Cell(vector<int>& combination, vector<User*>& users);

        void updateCell(vector<int>& combination, vector<User*>& users);

        void setCombinationIndex();

        bool weaklyDominates(Cell& other);

        double computeARF(vector<int>& indices);

        vector<pair<Slab,User*>>& getCombination();

        double getARF() const;

         vector<pair<int,User*>>& getCombinationIndex();

        void clean();
        
};
