#include "user.h"
#include "slab.h"
#include <unordered_set>

User::User(vector<double> coordinates) : coordinates(coordinates) {}

User::~User() {}

void User::addProduct(vector<double> product) {

    double score = CalculateScore(product);
    products.push_back(score);
    sort(products.begin(), products.end(),
     [&](double a, double b) { return a > b; });

}

double User::CalculateScore(vector<double> element){
    double score = 0;
    for (int i=0;i<element.size();i++) {
        score += coordinates[i] * element[i];
    }
    return score;
}

void User::createSlabs(){
    slabs.clear();
    products.push_back(minScore);
    products.insert(products.begin(), maxScore);
    for(int i=1;i < products.size();i++){
        slabs.push_back(Slab(Halfspace(coordinates,products[i-1]),Halfspace(coordinates,products[i]),i-1));

    }
}

void User::setEdgeScores(double minscore, double maxscore){
    minScore = minscore;
    maxScore = maxscore;
}

Slab& User::findClosestSlab(double score) {
    const double eps = 1e-12;
    int n = (int)slabs.size();
    int left = 0, right = n;

    // Binary search in DESCENDING order to find first index with v <= score
    while (left < right) {
        int mid = left + (right - left) / 2;
        double v = slabs[mid].getLower().getScore();

        if (std::fabs(v - score) <= eps) {
            // treat as exact match
            return slabs[mid];
        }
        if (v > score) {
            // v is greater than target -> move right to find <=
            left = mid + 1;
        } else {
            // v <= target -> keep this side
            right = mid;
        }
    }

    // left is the first index with v <= score (or n if all v > score)
    if (left < n) {
        return slabs[left];        // guaranteed v <= score
    } else {
        // Clamp scores below the represented range to the final slab.
        return slabs[n - 1];
    }
}

void User::checkProduct(double maxscore,double minscore){
    unordered_set<double> todelete;
    for(int i=0;i<products.size();i++){
        if(products[i]<minscore || products[i] > maxscore){
            todelete.insert(products[i]);
        }
    }

    if (!todelete.empty()) {
        products.erase(
            std::remove_if(products.begin(), products.end(),
                [&](double u) {
                    return todelete.find(u) != todelete.end();
                }),
            products.end());
    }
}

double User::getProductsScore(int index){
    return products[index];
}

Slab& User::getSlab(int index){
    return slabs[index];
}

int User::getSlabIndex(Slab item){
    auto it = find(slabs.begin(), slabs.end(), item);
    if (it != slabs.end()) {
        return distance(slabs.begin(), it);  // get index
    }
    return -1;
}

vector<double> User::getCoordinates(){
    return coordinates;
}

bool User::checkSlabEmpty(){
    return slabs.empty();
}

vector<double> User::getProducts(){
    return products;
}

int User::getProductsSize(){
    return products.size();
}

int User::getSlabSize(){
    return slabs.size();
}
