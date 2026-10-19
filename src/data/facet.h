#pragma once
#include <vector>
#include "../geometry/qhull.h"
using namespace std;

struct User;;


struct UserBorders{
        int minBorder,maxBorder;
        User* u;

         UserBorders(int min, int max, User* user) 
        : minBorder(min), maxBorder(max), u(user) {}
};

struct BestVector{
    vector<double> vertex;
    vector<User*> mostPreferenceUsers;
    vector<User*> leastPreferenceUsers;
};

class Facet
{

    private:
        FacetVRep facet;
        vector<vector<double>>  vertices;
        vector<UserBorders> borders;
        BestVector v;
        double lowerBound = 0, maxBound = INFINITY;
        vector<vector<User*>>  verticesPopularity;
        vector<vector<User*>> verticesUnpopularity;


    public:
        Facet() = default;
        ~Facet() = default;

        Facet(FacetVRep& facetVRep,vector<std::vector<double>> &all_vertices);

        inline bool compareFacets(const Facet& a, const Facet& b) {return a.getLowerBound() < b.getLowerBound();}

        void calculateRankFromUser(User* u);

        BestVector findBestVertex();

        double getLowerBound() const;

        double getMaxBound() const;

        FacetVRep& getFacetVRep();

        vector<int> getMinBorders();

        vector<int> getMaxBorders();

        const vector<vector<double>>& getVertices() const { return vertices; }

};