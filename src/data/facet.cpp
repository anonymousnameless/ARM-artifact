#include "facet.h"
#include "user.h"
#include <climits>

Facet::Facet(FacetVRep& facetVRep,vector<std::vector<double>> &all_vertices){
    facet = facetVRep;
    for(int i=0;i<facetVRep.vertex_indices.size();i++){
        vertices.push_back(facetVRep.get_bounded_vertices(all_vertices).at(i));
    }

    verticesPopularity.resize(vertices.size());
    verticesUnpopularity.resize(vertices.size());
}

void Facet::calculateRankFromUser(User* u){
    double minRank = INFINITY;
    double maxRank = -1;
    double tmpRank;
    double score;
    int bestVert, worstVert;
    for(int i=0;i<vertices.size();i++){
        int slabIndex = -1;
        u->findClosestSlab(u->CalculateScore(vertices[i]), &slabIndex);
        tmpRank = slabIndex;
        if(tmpRank<minRank){
            minRank = tmpRank;
            bestVert = i;
        }
        
        if(tmpRank>maxRank){
            maxRank = tmpRank;
            worstVert = i;
        }
    }
    lowerBound = lowerBound + minRank;
    maxBound = maxBound + maxRank;
    borders.push_back(UserBorders(minRank,maxRank,u));
    verticesPopularity[bestVert].push_back(u);
    verticesUnpopularity[worstVert].push_back(u);
}  


BestVector Facet::findBestVertex(){
    int maxIndex = -1;
    int maxSize = 0;


    for (int i = 0; i < verticesPopularity.size(); ++i) {
        if (verticesPopularity[i].size() > maxSize) {
            maxSize = verticesPopularity[i].size();
            maxIndex = i;
        }

    }

    v.vertex = vertices[maxIndex];
    for(int i=0;i<verticesPopularity[maxIndex].size();i++){
        v.mostPreferenceUsers.push_back(verticesPopularity[maxIndex][i]);
    }

    for (int i=0;i<verticesUnpopularity[maxIndex].size();i++){
        v.leastPreferenceUsers.push_back(verticesUnpopularity[maxIndex][i]);
    }

    return v;

}

FacetVRep& Facet::getFacetVRep(){
    return facet;
}


double Facet::getLowerBound() const {
    return lowerBound;
}

double Facet::getMaxBound() const {
    return maxBound;
}

vector<int> Facet::getMinBorders(){
    vector<int> minBorder;
    for(int i=0;i<borders.size();i++){
        minBorder.push_back(borders[i].minBorder);
    }

    return minBorder;

}

vector<int> Facet::getMaxBorders(){
    vector<int> maxBorder;
    for(int i=0;i<borders.size();i++){
        maxBorder.push_back(borders[i].maxBorder);
    }

    return maxBorder;
}
