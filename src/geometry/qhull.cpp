#include "qhull.h"
#include "../data/user.h"
#include "../utils/utils.h"

extern "C" {
#include <libqhull_r/qhull_ra.h>
}

// CPLEX for Phase I LP (find_feasible)
#include <ilcplex/cplex.h>

#include <iostream>
#include <sstream>
#include <list>
#include <map>
#include <algorithm>
#include <cmath>

// ===== Internal Helpers =====

namespace {

    static inline double snap0(double v, double eps) {
        return (std::fabs(v) < eps) ? 0.0 : v;
    }

    static inline double eval_leq(Halfspace& L, const std::vector<double>& x) {
        double s = 0.0;
        for (size_t i = 0; i < L.getCoords().size(); ++i)
            s += L.getCoords()[i] * x[i];
        return s - L.getScore();
    }

    // -------------------------------------------------------------------------
    // FIX 1: find_feasible using CPLEX Phase I LP
    //
    // Solves:  minimize s
    //          subject to: a_i·x - s <= b_i   for all halfspaces i
    //                      0 <= x_j <= 1       (unit box)
    //                      s unbounded
    //
    // If optimal s < -1e-8, then x* is a strictly interior feasible point.
    // This replaces the old heuristic anchor approach which silently failed
    // on tight polytopes, causing Qhull to produce wrong vertices.
    // -------------------------------------------------------------------------
    static bool find_feasible( std::vector<Halfspace>& I, int d, std::vector<double>& fp) {
        int status = 0;

        CPXENVptr env = CPXopenCPLEX(&status);
        if (!env) return false;

        CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_OFF);
        CPXsetintparam(env, CPXPARAM_Simplex_Display, 0);

        CPXLPptr lp = CPXcreateprob(env, &status, "feasible");
        if (!lp) {
            CPXcloseCPLEX(&env);
            return false;
        }

        // Variables: x_0 .. x_{d-1}, s (index d)
        // Objective: minimize s
        int numCols = d + 1;
        std::vector<double> obj(numCols, 0.0);
        obj[d] = 1.0; // only s in objective

        // Box: 0 <= x_i <= 1, s unbounded
        std::vector<double> lb(numCols, 0.0);
        std::vector<double> ub(numCols, 1.0);
        lb[d] = -CPX_INFBOUND;
        ub[d] =  CPX_INFBOUND;

        status = CPXnewcols(env, lp, numCols,
                            obj.data(), lb.data(), ub.data(),
                            nullptr, nullptr);
        if (status) {
            CPXfreeprob(env, &lp);
            CPXcloseCPLEX(&env);
            return false;
        }

        // Constraints: a_i·x - s <= b_i
        int nRows = (int)I.size();
        std::vector<int>    rmatbeg(nRows);
        std::vector<int>    rmatind;
        std::vector<double> rmatval;
        std::vector<double> rhs(nRows);
        std::vector<char>   sense(nRows, 'L');

        rmatind.reserve((size_t)nRows * (d + 1));
        rmatval.reserve((size_t)nRows * (d + 1));

        for (int i = 0; i < nRows; i++) {
            rmatbeg[i] = (int)rmatind.size();
            const auto& coords = I[i].getCoords();
            for (int j = 0; j < d; j++) {
                rmatind.push_back(j);
                rmatval.push_back(coords[j]);
            }
            // -s coefficient
            rmatind.push_back(d);
            rmatval.push_back(-1.0);
            rhs[i] = I[i].getScore();
        }

        status = CPXaddrows(env, lp, 0, nRows, (int)rmatind.size(),
                            rhs.data(), sense.data(),
                            rmatbeg.data(), rmatind.data(), rmatval.data(),
                            nullptr, nullptr);
        if (status) {
            CPXfreeprob(env, &lp);
            CPXcloseCPLEX(&env);
            return false;
        }

        status = CPXlpopt(env, lp);

        bool found = false;
        if (status == 0) {
            double objval = 0.0;
            CPXgetobjval(env, lp, &objval);

            if (objval < -1e-10) { // accept interior points with very small slack
                fp.resize(d);
                std::vector<double> sol(numCols);
                CPXgetx(env, lp, sol.data(), 0, numCols - 1);
                for (int i = 0; i < d; i++) fp[i] = sol[i];
                found = true;
            }
        }

        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return found;
    }

    static void add_box_constraints(std::vector<Halfspace>& ineq, int d, BoxConstraints box_type) {
        if (box_type == BoxConstraints::NONE) return;

        // Add lower bounds: x_i >= 0  →  -x_i <= 0
        for (int i = 0; i < d; i++) {
            Halfspace ge0;
            ge0.getCoords().assign(d, 0.0);
            ge0.getCoords()[i] = -1.0;
            ge0.setScore(0.0);
            ineq.push_back(ge0);
        }

        if (box_type == BoxConstraints::FULL_UNIT_BOX) {
            // Add upper bounds: x_i <= 1
            for (int i = 0; i < d; i++) {
                Halfspace le1;
                le1.getCoords().assign(d, 0.0);
                le1.getCoords()[i] = 1.0;
                le1.setScore(1.0);
                ineq.push_back(le1);
            }
        }
    }

    static void sort_facet_vertices_2d(FacetVRep& facet, const std::vector<std::vector<double>>& vertices) {
        if (vertices[0].size() != 2 || facet.vertex_indices.size() < 2) return;

        double cx = 0.0, cy = 0.0;
        for (int idx : facet.vertex_indices) {
            cx += vertices[idx][0];
            cy += vertices[idx][1];
        }
        cx /= facet.vertex_indices.size();
        cy /= facet.vertex_indices.size();

        std::sort(facet.vertex_indices.begin(), facet.vertex_indices.end(),
            [&](int a, int b) {
                double angle_a = std::atan2(vertices[a][1] - cy, vertices[a][0] - cx);
                double angle_b = std::atan2(vertices[b][1] - cy, vertices[b][0] - cx);
                return angle_a < angle_b;
            });
    }

    static void sort_facet_vertices_general(FacetVRep& facet, const std::vector<std::vector<double>>& vertices, int d) {
        if (facet.vertex_indices.size() < 3) return;

        std::vector<int> sorted;
        std::vector<bool> used(facet.vertex_indices.size(), false);

        sorted.push_back(facet.vertex_indices[0]);
        used[0] = true;

        for (size_t i = 1; i < facet.vertex_indices.size(); i++) {
            int last_idx = sorted.back();
            double min_dist = 1e100;
            int best_idx = -1;

            for (size_t j = 0; j < facet.vertex_indices.size(); j++) {
                if (used[j]) continue;
                double dist = 0.0;
                for (int k = 0; k < d; k++) {
                    double diff = vertices[facet.vertex_indices[j]][k] - vertices[last_idx][k];
                    dist += diff * diff;
                }
                if (dist < min_dist) {
                    min_dist = dist;
                    best_idx = (int)j;
                }
            }

            if (best_idx >= 0) {
                sorted.push_back(facet.vertex_indices[best_idx]);
                used[best_idx] = true;
            }
        }
        facet.vertex_indices = sorted;
    }

} // namespace

// ===== Main UBC Function =====

WCAResult get_ubc_vrep(int d, std::vector<Halfspace>& ineq_in, BoxConstraints box_type, double eps) {
    WCAResult result;

    // 1. Combine input constraints with box constraints
    std::vector<Halfspace> full_constraints = ineq_in;
    add_box_constraints(full_constraints, d, box_type);

    result.all_inequalities = full_constraints;

    // 2. Find feasible interior point via CPLEX Phase I LP (FIX 1)
    std::vector<double> feas;
    if (!find_feasible(full_constraints, d, feas))
        return result; // polytope has no interior — return empty

    // 3. Prepare Qhull halfspace data
    const int dimHS = d + 1;
    const int nHS   = (int)full_constraints.size();
    std::vector<coordT> hs;
    hs.reserve((size_t)nHS * dimHS);

    for (auto& L : full_constraints) {
        for (int j = 0; j < d; j++)
            hs.push_back((coordT)L.getCoords()[j]);
        hs.push_back((coordT)(-L.getScore()));
    }

    // 4. Run Qhull halfspace intersection
    std::ostringstream oss;
    oss << "qhull H";
    for (int i = 0; i < d; i++) {
        if (i) oss << ",";
        oss << feas[i];
    }
    oss << " Qt";
    if (d >= 5) oss << " Qx";

    qhT qh_qh;
    qhT* qh = &qh_qh;
    qh_zero(qh, stderr);
    int exitcode = qh_new_qhull(qh, dimHS, nHS, hs.data(), qh_False,
                                (char*)oss.str().c_str(), nullptr, stderr);

    if (exitcode != 0) {
        qh_freeqhull(qh, !qh_ALL);
        int a = 0, b = 0;
        qh_memfreeshort(qh, &a, &b);
        return result;
    }

    // 5. Collect raw vertices — FIX 2: proper O(log V) map lookup
    std::vector<std::vector<double>> raw_vertices;

    auto vec_less = [&](const std::vector<double>& a, const std::vector<double>& b) {
        for (size_t i = 0; i < a.size(); i++) {
            if (a[i] < b[i] - eps) return true;
            if (a[i] > b[i] + eps) return false;
        }
        return false;
    };
    std::map<std::vector<double>, int, decltype(vec_less)> vertex_map(vec_less);

    facetT* facet;
    FORALLfacets {
        double off = facet->offset;
        if (std::fabs(off) < 1e-18) continue;

        std::vector<double> vertex(d);
        for (int i = 0; i < d; i++)
            vertex[i] = snap0(-facet->normal[i] / off + feas[i], eps);

        // O(log V) lookup — no linear scan
        if (vertex_map.find(vertex) == vertex_map.end()) {
            vertex_map[vertex] = (int)raw_vertices.size();
            raw_vertices.push_back(vertex);
        }
    }

    qh_freeqhull(qh, !qh_ALL);
    int curl = 0, tot = 0;
    qh_memfreeshort(qh, &curl, &tot);

    const double hyperplane_eps = 1e-6;

    // 6. Filter vertices — keep only those touching a user inequality
    for (auto& v : raw_vertices) {
        bool touches_user_ineq = false;

        for (auto& ineq : ineq_in) {
            double score = 0.0;
            const auto& coords = ineq.getCoords();
            for (size_t i = 0; i < (size_t)d; i++)
                score += coords[i] * v[i];

            if (std::fabs(score - ineq.getScore()) < hyperplane_eps) {
                touches_user_ineq = true;
                break;
            }
        }

        if (touches_user_ineq)
            result.vertices.push_back(v);
    }

    // Exact fallback when Qhull produced raw vertices but the active-plane
    // filter lost every input halfspace vertex.
    if (result.vertices.empty() && !raw_vertices.empty()) {
        auto solve_active_set = [&](const std::vector<int>& active,
                                    std::vector<double>& x) {
            std::vector<std::vector<double>> a(d, std::vector<double>(d + 1, 0.0));
            for (int row = 0; row < d; ++row) {
                const auto& constraint = full_constraints[active[row]];
                const auto& coords = constraint.getCoords();
                for (int col = 0; col < d; ++col)
                    a[row][col] = coords[col];
                a[row][d] = constraint.getScore();
            }

            for (int col = 0; col < d; ++col) {
                int pivot = col;
                for (int row = col + 1; row < d; ++row) {
                    if (std::fabs(a[row][col]) > std::fabs(a[pivot][col]))
                        pivot = row;
                }
                if (std::fabs(a[pivot][col]) < 1e-10)
                    return false;
                if (pivot != col)
                    std::swap(a[pivot], a[col]);

                const double div = a[col][col];
                for (int j = col; j <= d; ++j)
                    a[col][j] /= div;
                for (int row = 0; row < d; ++row) {
                    if (row == col)
                        continue;
                    const double factor = a[row][col];
                    for (int j = col; j <= d; ++j)
                        a[row][j] -= factor * a[col][j];
                }
            }

            x.assign(d, 0.0);
            for (int row = 0; row < d; ++row)
                x[row] = a[row][d];
            return true;
        };

        auto satisfies_all_constraints = [&](const std::vector<double>& v) {
            for (auto& ineq : full_constraints) {
                double score = 0.0;
                const auto& coords = ineq.getCoords();
                for (int i = 0; i < d; ++i)
                    score += coords[i] * v[i];
                if (score > ineq.getScore() + 1e-8)
                    return false;
            }
            return true;
        };

        auto touches_input_halfspace = [&](const std::vector<double>& v) {
            for (auto& ineq : ineq_in) {
                double score = 0.0;
                const auto& coords = ineq.getCoords();
                for (int i = 0; i < d; ++i)
                    score += coords[i] * v[i];
                if (std::fabs(score - ineq.getScore()) < hyperplane_eps)
                    return true;
            }
            return false;
        };

        auto already_added = [&](const std::vector<double>& v) {
            for (const auto& existing : result.vertices) {
                double max_delta = 0.0;
                for (int i = 0; i < d; ++i)
                    max_delta = std::max(max_delta, std::fabs(existing[i] - v[i]));
                if (max_delta < 1e-8)
                    return true;
            }
            return false;
        };

        std::vector<int> active(d, 0);
        auto enumerate_active_sets = [&](auto&& self, int start, int depth) -> void {
            if (depth == d) {
                std::vector<double> vertex;
                if (solve_active_set(active, vertex) &&
                    satisfies_all_constraints(vertex) &&
                    touches_input_halfspace(vertex) &&
                    !already_added(vertex)) {
                    result.vertices.push_back(vertex);
                }
                return;
            }

            const int remaining = d - depth;
            for (int idx = start;
                 idx <= (int)full_constraints.size() - remaining;
                 ++idx) {
                active[depth] = idx;
                self(self, idx + 1, depth + 1);
            }
        };
        enumerate_active_sets(enumerate_active_sets, 0, 0);
    }

    // 7. Filter facets — build FacetVRep for each user inequality
    for (size_t k = 0; k < ineq_in.size(); ++k) {
        auto& inequality = ineq_in[k];
        std::vector<int> vertices_on_plane;

        for (size_t v_idx = 0; v_idx < result.vertices.size(); v_idx++) {
            const auto& vertex_coords = result.vertices[v_idx];
            double score = 0.0;
            const auto& coords = inequality.getCoords();
            for (size_t i = 0; i < coords.size(); ++i)
                score += coords[i] * vertex_coords[i];

            if (std::fabs(score - inequality.getScore()) < hyperplane_eps)
                vertices_on_plane.push_back((int)v_idx);
        }

        if (vertices_on_plane.size() >= (size_t)(d - 1)) {
            FacetVRep f;
            f.vertex_indices    = vertices_on_plane;
            f.defining_halfspace = inequality;
            f.hyperplane_normal  = inequality.getCoords();
            f.hyperplane_offset  = inequality.getScore();
            f.original_index     = (int)k;

            f.representative.assign(d, 0.0);
            for (int idx : f.vertex_indices)
                for (int j = 0; j < d; j++)
                    f.representative[j] += result.vertices[idx][j];
            for (int j = 0; j < d; j++)
                f.representative[j] /= f.vertex_indices.size();

            result.facets.push_back(f);
        }
    }

    // 8. Sort vertices within each facet
    for (auto& f : result.facets) {
        if (d == 2)
            sort_facet_vertices_2d(f, result.vertices);
        else
            sort_facet_vertices_general(f, result.vertices, d);
    }

    return result;
}

// ===== Convex Hull Functions =====

HullVRep convex_hull_vrep_nd(int d, const std::vector<double>& pts) {
    HullVRep R;
    R.d = d;
    const int n = (int)pts.size() / d;
    if (n == 0 || d < 2) return R;

    const double eps = 1e-12;
    double s0 = 0.0;
    for (int j = 0; j < d; ++j) s0 += pts[j];
    bool on_simplex = true;
    for (int i = 1; i < n && on_simplex; ++i) {
        double s = 0.0;
        for (int j = 0; j < d; ++j) s += pts[i * d + j];
        if (std::fabs(s - s0) > eps) on_simplex = false;
    }

    int d_eff = d;
    std::vector<double> data_vec;
    if (on_simplex) {
        d_eff = d - 1;
        if (d_eff >= 2) {
            data_vec.resize((size_t)n * d_eff);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < d_eff; ++j)
                    data_vec[i * d_eff + j] = pts[i * d + j];
        } else {
            int imin = 0, imax = 0;
            for (int i = 1; i < n; ++i) {
                if (pts[i * d] < pts[imin * d]) imin = i;
                if (pts[i * d] > pts[imax * d]) imax = i;
            }
            R.indices = {imin, imax};
            for (int j = 0; j < d; ++j) R.coords.push_back(pts[imin * d + j]);
            for (int j = 0; j < d; ++j) R.coords.push_back(pts[imax * d + j]);
            return R;
        }
    } else {
        data_vec = pts;
    }

    // FIX 3: safe coordT cast with static_assert
    static_assert(sizeof(coordT) == sizeof(double), "Qhull coordT must be double");
    coordT* data = reinterpret_cast<coordT*>(const_cast<double*>(data_vec.data()));

    std::string flags = (d_eff >= 5) ? "qhull Qt Qx" : "qhull Qt";

    qhT qh_qh;
    qhT* qh = &qh_qh;
    qh_zero(qh, stderr);
    int exitcode = qh_new_qhull(qh, d_eff, n, data, qh_False,
                                (char*)flags.c_str(), nullptr, stderr);

    if (exitcode) {
        qh_freeqhull(qh, !qh_ALL);
        int a = 0, b = 0;
        qh_memfreeshort(qh, &a, &b);
        qh_zero(qh, stderr);
        std::string flags2 = flags + " QJ";
        exitcode = qh_new_qhull(qh, d_eff, n, data, qh_False,
                                (char*)flags2.c_str(), nullptr, stderr);
        if (exitcode) {
            qh_freeqhull(qh, !qh_ALL);
            a = b = 0;
            qh_memfreeshort(qh, &a, &b);
            return R;
        }
    }

    vertexT* vertex;
    FORALLvertices {
        int id = qh_pointid(qh, vertex->point);
        R.indices.push_back(id);
        for (int j = 0; j < d; ++j)
            R.coords.push_back(pts[id * d + j]);
    }

    qh_freeqhull(qh, !qh_ALL);
    int curl = 0, tot = 0;
    qh_memfreeshort(qh, &curl, &tot);
    return R;
}

HullVRep convex_hull_vrep_from_users(std::vector<User*>& users) {
    HullVRep R;
    if (users.empty()) return R;

    const int d = (int)users[0]->getCoordinates().size();
    for (auto& u : users)
        if ((int)u->getCoordinates().size() != d) return R;

    std::vector<double> flat;
    flat.reserve(users.size() * (size_t)d);
    for (auto& u : users) {
        const auto& coords = u->getCoordinates();
        flat.insert(flat.end(), coords.begin(), coords.end());
    }

    R = convex_hull_vrep_nd(d, flat);
    return R;
}

void hull_vertices_users(const HullVRep& H, std::vector<User*>& users, std::vector<User*>& hullUsers) {
    if (H.d <= 0) return;
    hullUsers.reserve(H.indices.size());
    for (int idx : H.indices)
        hullUsers.push_back(users[idx]);
}
