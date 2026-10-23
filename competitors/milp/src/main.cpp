#include "param.h"
#include <iostream>
#include <vector>
#include "user.h"
#include "slab.h"
#include "halfspace.h"
#include "cell.h"
#include <cmath>
#include <queue>
#include <unordered_map>
#include <random>
#include <fstream>
#include <cstdint>
#include <chrono>
#include <stdexcept>
#include <ilcplex/ilocplex.h>
#include "rentry.h"
#include "rnode.h"
#include "rtree.h"
#include "qhull.h"
#include <set>
#include <iomanip>
#include <unordered_set>

using namespace std;

Rtree* loadProducts(const char* productFile, int dim);
bool loadUsers(vector<User> &usersLists,const char* userFile,int dim);
void CreateWCA(vector<Halfspace>& inequalities, int d, int N, double r, uint64_t seed);
void rteePruning(User &u, RtreeNode* node, int dim, double minscore, double maxscore);
Cell findMinimumARF_MILP(vector<User>& userLists, vector<Halfspace>& inequalities, int dimen, vector<double>& solution);
int CheckCoordinateDomination(const vector<double>& a, const vector<double>& b);
void TrimWcaEdges(vector<std::vector<double>> intersections, vector<User>& users, vector<std::vector<double>> &TrimIntersections);
static int computeMaxChild(int pageSize, int dimen);

struct VectorHash {
    size_t operator()(const vector<int>& v) const {
        size_t h = 0;
        size_t base = 131;
        for (int x : v)
            h = h * base + x;
        return h;
    }
};

int main(const int argc, const char** argv)
{
    const auto started = std::chrono::steady_clock::now();



	int dim = atoi(Param::read(argc, argv, "-d", "")); //number of Dimensions
	int N = atoi(Param::read(argc, argv, "-N", "")); // Number of constrains WCA
    const char* productFile = Param::read(argc, argv, "-pf", "");
    const char* userFile = Param::read(argc, argv, "-uf", "");
    if (dim < 2 || dim > MAXDIMEN || !*productFile || !*userFile) {
        cerr << "Usage: -d DIM -pf PRODUCTS -uf USERS [-N PLANES -wr RADIUS -ws SEED]" << endl;
        return 1;
    }
    if (!std::ifstream(productFile).good() || !std::ifstream(userFile).good()) {
        cerr << "Cannot open input file" << endl;
        return 1;
    }
    if (N <= 0) {
        N = dim;
    }
	double r = atof(Param::read(argc, argv, "-wr", "1.0"));
    if (!std::isfinite(r) || r <= 0.0) {
        cerr << "Radius must be finite and positive" << endl;
        return 1;
    }
    uint64_t wcaSeed =
        strtoull(Param::read(argc, argv, "-ws", "31"), nullptr, 10);
	vector<User> usersLists;
	vector<Halfspace> inequalities;   // Store inequalities as pair: (coefficients vector, RHS)
    vector<vector<double>> intersections;
    vector<Halfspace> H;
    vector<double> solution;


    if (!loadUsers(usersLists,userFile,dim)) {
        cerr << "Malformed user input" << endl;
        return 1;
    }
    if (usersLists.empty()) {
        cerr << "No users loaded" << endl;
        return 1;
    }
    Rtree* rtree = loadProducts(productFile, dim);
    if (!rtree) return 1;
    CreateWCA(inequalities, dim, N, r, wcaSeed);


    H.resize(inequalities.size());

    for (int i = 0; i < inequalities.size(); i++) {
        for (int j = 0; j <inequalities[i].getCoords().size(); j++) {
            H[i].getCoords().push_back(inequalities[i].getCoords()[j]);
        }
        H[i].setScore(inequalities[i].getScore());

    }

    intersections = vrep_halfspace_qhull(dim, H, true);

    std::puts("V-representation (one vertex per line):");
    for (auto& p: intersections){
        for (int i=0;i<dim;i++) std::printf("%s%.12g", (i?" ":"  "), p[i]);
        std::puts("");
    }


    const double EPS = 1e-16; // smaller than any meaningful number
    intersections.erase(
        std::remove_if(intersections.begin(), intersections.end(),
            [EPS](const std::vector<double>& v) {
                return std::all_of(v.begin(), v.end(),
                                [EPS](double x){ return std::fabs(x) < EPS; });
            }),
        intersections.end()
    );

    vector<vector<double>> TrimIntersections;
    TrimWcaEdges(intersections, usersLists, TrimIntersections);
    if (TrimIntersections.empty()) {
        TrimIntersections = intersections;
    }
    cout << "Trimmed UBS vertices: " << TrimIntersections.size()
         << " / " << intersections.size() << endl;

    vector<User> nonSlabsUsers;

    for(int i=0;i<usersLists.size();i++){
        double maxscore = -INFINITY;
        double minscore = INFINITY;
        for (const auto& pt : TrimIntersections) {
            double score = usersLists[i].CalculateScore(pt);
            maxscore = std::max(maxscore, score);
        }
        for (const auto& pt : intersections) {
            double score = usersLists[i].CalculateScore(pt);
            minscore = std::min(minscore, score);
        }
        if (!std::isfinite(maxscore)) {
            maxscore = minscore;
        }
        rteePruning(usersLists[i], rtree->m_root, dim, minscore,maxscore);
        usersLists[i].setEdgeScores(minscore,maxscore);
        usersLists[i].createSlabs();
        if(usersLists[i].checkSlabEmpty())
            nonSlabsUsers.push_back(usersLists[i]);
    }

   usersLists.erase(
    std::remove_if(usersLists.begin(), usersLists.end(),
        [&](User& u) {
            return std::find(nonSlabsUsers.begin(), nonSlabsUsers.end(), u) != nonSlabsUsers.end();
        }),
    usersLists.end());

    Cell bestCombination = findMinimumARF_MILP(usersLists, inequalities, dim, solution);


	cout << "Best combination of products indices for users: ";
    for (int idx : bestCombination.getCombinationIndex()) {
        cout << idx << " ";
    }
    if ((int)solution.size() != dim) {
        cerr << "No feasible solution was returned" << endl;
        delete rtree;
        return 1;
    }
    cout << "\nFeasible point in the intersection: " << std::setprecision(17);
    for (double val : solution) {
        cout << val << " ";
    }
    cout << endl;

    long long rankOffset = 0;
    for (const User& user : usersLists) rankOffset += user.getRankOffset();
    cout << "\nARF: " << (long long)bestCombination.getARF() + rankOffset << endl;
    cout << "TIME_SEC: " << std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count() << endl;
    delete rtree;

	return 0;
}



Cell findMinimumARF_MILP(vector<User>& userLists, vector<Halfspace>& inequalities, int dimen, vector<double>& solution) {
    const int numUsers = (int)userLists.size();
    if (numUsers == 0) {
        solution.assign(dimen, 0.0);
        return Cell(vector<int>());
    }

    int totalSlabVars = 0;
    for (int i = 0; i < numUsers; ++i) {
        const int slabCount = userLists[i].getSlabSize();
        if (slabCount <= 0) {
            cerr << "Cannot build MILP: user " << i << " has no slabs" << endl;
            return {};
        }
        totalSlabVars += slabCount;
    }

    const int totalVars = dimen + totalSlabVars;
    IloEnv env;
    try {
        IloModel model(env);
        IloNumVarArray x(env, dimen, 0.0, 1.0, ILOFLOAT);
        vector<IloBoolVarArray> choice(numUsers);

        for (int i = 0; i < numUsers; ++i) {
            choice[i] = IloBoolVarArray(env, userLists[i].getSlabSize());
        }

        IloExpr objective(env);
        for (int i = 0; i < numUsers; ++i) {
            IloExpr oneSlab(env);
            for (int k = 0; k < userLists[i].getSlabSize(); ++k) {
                oneSlab += choice[i][k];
                objective += userLists[i].getSlab(k).getIndex() * choice[i][k];
            }
            model.add(oneSlab == 1);
            oneSlab.end();
        }
        model.add(IloMinimize(env, objective));
        objective.end();

        // If y_{i,k}=1, force the chosen point x to lie in user i's slab k:
        // lower_{i,k} <= w_i · x <= upper_{i,k}.
        for (int i = 0; i < numUsers; ++i) {
            vector<double> coordinates = userLists[i].getCoordinates();

            double scoreMinPossible = 0.0;
            double scoreMaxPossible = 0.0;
            for (double c : coordinates) {
                if (c >= 0.0) {
                    scoreMaxPossible += c;
                } else {
                    scoreMinPossible += c;
                }
            }
            const double bigM =
                std::max(1.0, std::fabs(scoreMaxPossible - scoreMinPossible) + 1.0);

            IloExpr dot(env);
            for (int j = 0; j < dimen; ++j) {
                dot += coordinates[j] * x[j];
            }

            for (int k = 0; k < userLists[i].getSlabSize(); ++k) {
                const double lower = userLists[i].getSlab(k).getLower().getScore();
                const double upper = userLists[i].getSlab(k).getUpper().getScore();
                model.add(dot >= lower - bigM * (1 - choice[i][k]));
                model.add(dot <= upper + bigM * (1 - choice[i][k]));
            }
            dot.end();
        }

        // The point x must remain inside the WCA.
        for (auto& ineq : inequalities) {
            IloExpr lhs(env);
            vector<double>& coeffs = ineq.getCoords();
            for (int j = 0; j < dimen; ++j) {
                lhs += coeffs[j] * x[j];
            }
            model.add(lhs <= ineq.getScore());
            lhs.end();
        }

        cout << "Solving whole-problem CPLEX MILP with " << totalVars
             << " variables (" << totalSlabVars
             << " binary slab variables)..." << endl;

        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        cplex.setParam(IloCplex::Param::Threads, 1);
        cplex.setParam(IloCplex::Param::Emphasis::MIP, 1);
        cplex.setParam(IloCplex::Param::Emphasis::Numerical, 1);
        cplex.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-9);
        cplex.setParam(IloCplex::Param::MIP::Tolerances::Integrality, 1e-9);
        cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.0);
        cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.0);

        if (!cplex.solve() || cplex.getStatus() != IloAlgorithm::Optimal) {
            cerr << "CPLEX MILP did not solve. Status: " << cplex.getStatus()
                 << endl;
            env.end();
            return {};
        }

        solution.assign(dimen, 0.0);
        for (int j = 0; j < dimen; ++j) {
            solution[j] = cplex.getValue(x[j]);
        }

        vector<int> combination(numUsers, 0);
        for (int i = 0; i < numUsers; ++i) {
            int chosenSlab = 0;
            double bestValue = -1.0;
            for (int k = 0; k < userLists[i].getSlabSize(); ++k) {
                const double value = cplex.getValue(choice[i][k]);
                if (value > bestValue) {
                    bestValue = value;
                    chosenSlab = k;
                }
            }
            combination[i] = chosenSlab;
        }

        cout << "CPLEX MILP objective: " << cplex.getObjValue() << endl;
        env.end();
        return Cell(combination);
    } catch (IloException& e) {
        cerr << "CPLEX exception in whole-problem MILP: " << e.getMessage()
             << endl;
        env.end();
        return {};
    } catch (...) {
        cerr << "Unknown exception in whole-problem CPLEX MILP" << endl;
        env.end();
        return {};
    }
}


int CheckCoordinateDomination(const vector<double>& a, const vector<double>& b) {
    bool aBetter = false;
    bool bBetter = false;

    for (int i = 0; i < (int)a.size(); ++i) {
        if (a[i] > b[i]) {
            aBetter = true;
        } else if (a[i] < b[i]) {
            bBetter = true;
        }
    }

    if (aBetter && !bBetter) return 1;
    if (bBetter && !aBetter) return -1;
    return 0;
}

void TrimWcaEdges(vector<std::vector<double>> intersections,
                  vector<User>& users,
                  vector<std::vector<double>>& TrimIntersections) {
    TrimIntersections.clear();
    if (intersections.empty() || users.empty()) {
        TrimIntersections = intersections;
        return;
    }

    struct ScoredIntersection {
        vector<double> score;
        int index;
    };

    vector<ScoredIntersection> scoredIntersections;
    scoredIntersections.reserve(intersections.size());

    for (int i = 0; i < (int)intersections.size(); ++i) {
        vector<double> currentScore(users.size(), 0.0);
        for (int j = 0; j < (int)users.size(); ++j) {
            currentScore[j] = users[j].CalculateScore(intersections[i]);
        }
        scoredIntersections.push_back({currentScore, i});
    }

    sort(scoredIntersections.begin(), scoredIntersections.end(),
         [](const ScoredIntersection& a, const ScoredIntersection& b) {
             if (a.score != b.score) {
                 return a.score > b.score;
             }
             return a.index < b.index;
         });

    vector<vector<double>> keptScores;
    for (const ScoredIntersection& current : scoredIntersections) {
        bool isDominated = false;
        for (const vector<double>& kept : keptScores) {
            if (CheckCoordinateDomination(kept, current.score) == 1) {
                isDominated = true;
                break;
            }
        }
        if (!isDominated) {
            keptScores.push_back(current.score);
            TrimIntersections.push_back(intersections[current.index]);
        }
    }
}

static long long countProducts(const RtreeNode* node) {
    long long count = 0;
    for (int i = 0; i < node->m_count; ++i) {
        const RtreeNodeEntry* entry = node->m_entries[i];
        count += entry->m_child ? countProducts(entry->m_child) : 1;
    }
    return count;
}

void rteePruning(User &u, RtreeNode* node, int dim, double minscore, double maxscore){
   if (!node) return;

    for (int i = 0; i < node->m_count; i++) {
        RtreeNodeEntry* entry = node->m_entries[i];

        // Prepare Lower and Upper points
        vector<double> Lowerpt(dim), Upperpt(dim);
        for (int j = 0; j < dim; j++) {
            Lowerpt[j] = entry->m_mbr.getLower().m_coor[j];
            Upperpt[j] = entry->m_mbr.getUpper().m_coor[j];
        }

        double LowerScore = u.CalculateScore(Lowerpt);
        if (LowerScore > maxscore) {
            // These products outrank every feasible point and contribute a constant.
            u.addRankOffset(entry->m_child ? countProducts(entry->m_child) : 1);
            continue;
        }

        double UpperScore = u.CalculateScore(Upperpt);
        if (LowerScore < minscore && UpperScore < minscore) continue;  // prune

        // Leaf node → add product
        if (entry->m_child == nullptr) {
            Point center = entry->m_mbr.getCenter();
            u.addProduct(vector<double>(center.m_coor, center.m_coor + center.m_dimen));
        }
        // Internal node → recurse
        else {
            rteePruning(u, entry->m_child, dim, minscore, maxscore);
        }
    }

    return;
}

Rtree* loadProducts(const char* productFile, int dim)
{
    std::fstream fpdata(productFile, std::ios::in);
    if (!fpdata.is_open()) {
        std::cerr << "Failed to open file: " << productFile << "\n";
        return nullptr;
    }

    int objCnt = 0;
    RtreeNodeEntry** p = new RtreeNodeEntry*[MAXPTS]; // ensure MAXPTS >= #lines

    while (true) {
        int id;
        fpdata >> id;
        if (!fpdata) break; // handles EOF and errors

        std::vector<float> cl(dim), cu(dim);
        for (int d = 0; d < dim; ++d) fpdata >> cl[d];
        for (int d = 0; d < dim; ++d) fpdata >> cu[d];
        bool validBounds = true;
        for (int d = 0; d < dim; ++d) {
            validBounds = validBounds && std::isfinite(cl[d]) &&
                std::isfinite(cu[d]) && cl[d] >= 0.0f &&
                cl[d] <= cu[d] && cu[d] <= 1.0f;
        }
        if (!fpdata || !validBounds || objCnt >= MAXPTS) {
            cerr << "Malformed product row or product capacity exceeded" << endl;
            for (int i = 0; i < objCnt; ++i) delete p[i];
            delete[] p;
            return nullptr;
        }

        Hypercube hc(dim, cl.data(), cu.data());
        p[objCnt++] = new RtreeNodeEntry(id, hc);

        if ((objCnt % 1000) == 0) std::cout << ".";
        if ((objCnt % 10000) == 0) std::cout << objCnt << " objects loaded\n";
    }
    fpdata.close();
    if (objCnt == 0) {
        cerr << "No products loaded" << endl;
        delete[] p;
        return nullptr;
    }

    std::cout << "Bulkloading R-tree...\n";
    int maxChild = computeMaxChild(4096, dim); // OK
    Rtree* rtree = new Rtree(dim, maxChild, maxChild, (int)(maxChild*0.3f), (int)(maxChild*0.3f));
    rtree->bulkLoad(p, objCnt, maxChild);
    std::cout << "[Rtree build done]\n";

    // entries are owned by the tree now; free the pointer array
    delete[] p;

    return rtree;
}

static int computeMaxChild(int pageSize, int dimen) {
    int nodeHeaderSize = sizeof(bool) + 3 * sizeof(int) + sizeof(RtreeNode*);
    int entrySize = Hypercube::size(dimen) + sizeof(int) + sizeof(RtreeNode*);
    return (pageSize - nodeHeaderSize) / entrySize;
}





bool loadUsers(vector<User>& usersLists, const char* userFile, int dim) {
    std::ifstream input(userFile);
    int id;
    while (input >> id) {
        vector<double> coordinates(dim);
        double sum = 0.0;
        for (double& weight : coordinates) {
            if (!(input >> weight) || !std::isfinite(weight) || weight < 0.0)
                return false;
            weight = std::floor(weight * 10000.0) / 10000.0;
            sum += weight;
        }
        if (sum <= 0.0) return false;
        usersLists.emplace_back(coordinates);
    }
    return input.eof();
}

void CreateWCA(vector<Halfspace>& inequalities, int d, int N, double r, uint64_t seed){
    std::mt19937_64 gen(seed);
    std::normal_distribution<> dist(0.0, 1.0);  // standard normal

    for(int i = 0; i < N; ++i) {
        std::vector<double> z(d);
        double norm = 0.0;

        // Step 3: generate standard normals
        for(int j = 0; j < d; ++j) {
            double val = std::abs(dist(gen)); // positive orthant
            z[j] = val;
            norm += val * val;
        }

        norm = std::sqrt(norm);

        // normalize to lie on sphere
        for(int j = 0; j < d; ++j) {
            z[j] = z[j] * r / norm;
        }

        // Step 4: store inequality (tangent hyperplane)
        inequalities.push_back({z, r*r});
    }

    // print
    for(int i = 0; i < N; ++i) {
        std::cout << "Inequality " << i+1 << ": ";
        for(int j = 0; j < d; ++j) {
            std::cout << inequalities[i].getCoords()[j] << "*x" << j+1;
            if(j != d-1) std::cout << " + ";
        }
        std::cout << " <= " << inequalities[i].getScore() << "\n";
    }

    return;
}
