#pragma once
#include <vector>

using namespace std;

class Cell
{
     private:
        vector<int> combinationIndex;
        double arf = -1;

    public:
        Cell() = default;
        ~Cell() = default;

        Cell(const vector<int>& combination);

        void updateCell(const vector<int>& combination);

        void setCombinationIndex();

        bool weaklyDominates(const Cell& other) const;

        double computeARF(const vector<int>& indices) const;

        double getARF();

        const vector<int>& getCombinationIndex() const;

        const vector<int>& getCombinationIndexOnly() const;

        void clean();

};
