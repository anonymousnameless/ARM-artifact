#pragma once

using namespace std;

#include <vector>


class Halfspace {

    private:
        vector<double> coords;
        double score;

    public:
        Halfspace() = default;
        ~Halfspace() = default;

        Halfspace(vector<double> coords, double score);

        friend bool operator==(const Halfspace& a, const Halfspace& b) {
            return a.coords == b.coords && a.score == b.score;
        }

        Halfspace& operator=(const Halfspace& other);

        void setScore(double score);

        vector<double>& getCoords();
        double getScore();

};
