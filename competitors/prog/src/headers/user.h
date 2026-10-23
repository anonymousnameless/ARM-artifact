#pragma once

using namespace std;

#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>

#include "slab.h"

class User
{

    private:
        vector<double> coordinates;
        vector<double> products;
        vector<Slab> slabs;
        double minScore;
        double maxScore;
        long long rankOffset = 0;

    public:
        User(vector<double> coordinates);
        ~User();
        void addRankOffset(long long count) { rankOffset += count; }
        long long getRankOffset() const { return rankOffset; }

        bool operator==(const User& other) const { return coordinates == other.coordinates; }

        void addProduct(vector<double> product);

        void setEdgeScores(double minscore, double maxscore);

        double CalculateScore(vector<double> element);

        void createSlabs();

        Slab& findClosestSlab(double score);

        void checkProduct(double maxscore, double minscore);

        double getProductsScore(int index);

        int getProductsSize();

        Slab& getSlab(int index);

        int getSlabIndex(Slab item);

        int getSlabSize();

        vector<double> getCoordinates();

        bool checkSlabEmpty();

        vector<double> getProducts();

};
