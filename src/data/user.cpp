#include "user.h"
#include "slab.h"
#include <unordered_set>
#include <stdexcept>
#include <string>

User::User(vector<double> coordinates) : coordinates(coordinates) {}

User::User(const User& other)
    : coordinates(other.coordinates),
      products(other.products),
      slabs(other.slabs),
      minScore(other.minScore),
      maxScore(other.maxScore), rankOffset(other.rankOffset) {
    rebindSlabCoordinates();
}

User& User::operator=(const User& other) {
    if (this != &other) {
        coordinates = other.coordinates;
        products = other.products;
        slabs = other.slabs;
        minScore = other.minScore;
        maxScore = other.maxScore;
        rankOffset = other.rankOffset;
        rebindSlabCoordinates();
    }
    return *this;
}

User::~User() {}

void User::addProduct(vector<double> product) {
    // Just accumulate — sorted once in createSlabs()
    double score = CalculateScore(product);
    products.push_back(score);
}

void User::addProductScore(double score) {
    // Push a pre-computed score directly — avoids dot product + vector allocation
    products.push_back(score);
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
    // Sort once here instead of after every addProduct call (was O(N^2 log N))
    sort(products.begin(), products.end(),
         [](double a, double b) { return a > b; });
    products.push_back(minScore);
    products.insert(products.begin(), maxScore);
    slabs.reserve(products.size() > 0 ? products.size() - 1 : 0);
    for(int i=1;i < (int)products.size();i++){
        slabs.push_back(Slab(coordinates, products[i-1], products[i], i-1));
    }
    vector<double>().swap(products);
}

void User::setEdgeScores(double minscore, double maxscore){
    minScore = minscore;
    maxScore = maxscore;
}

Slab& User::findClosestSlab(double score, int* indexOut) {
    const double eps = 1e-12;
    int n = (int)slabs.size();
    int left = 0, right = n;

    // Binary search in DESCENDING order to find first index with v <= score
    while (left < right) {
        int mid = left + (right - left) / 2;
        double v = slabs[mid].getLowerScore();

        if (std::fabs(v - score) <= eps) {
            // treat as exact match
            if (indexOut) *indexOut = mid;
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
        if (indexOut) *indexOut = left;
        return slabs[left];        // guaranteed v <= score
    } else {
        // Clamp scores below the represented range to the final slab.
        if (indexOut) *indexOut = n - 1;
        return slabs[n - 1];
    }
}

void User::checkProduct(double maxscore,double minscore){
    unordered_set<double> todelete;
    for(int i=0;i<products.size();i++){
        if (products[i] > maxscore) ++rankOffset;
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
    if (index < 0 || index >= (int)slabs.size()) {
        throw std::out_of_range("User::getSlab index out of bounds: " + std::to_string(index) + 
                                " for user with " + std::to_string(slabs.size()) + " slabs.\n" +
                                "If this occurs during B&B, it indicates a bug in border checking.");
    }
    return slabs[index];
}

int User::getSlabIndex(Slab item){
    auto it = find(slabs.begin(), slabs.end(), item);
    if (it != slabs.end()) {
        return distance(slabs.begin(), it);  // get index   
    }
    return -1;
}

const vector<double>& User::getCoordinates() const {
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

void User::rebindSlabCoordinates() {
    for (Slab& slab : slabs) {
        slab.rebindCoordinates(coordinates);
    }
}
