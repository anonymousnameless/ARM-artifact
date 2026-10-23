#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdio>

extern "C" {
  #include <libqhull_r/qhull_ra.h>   // sudo apt install -y libqhull-dev
}


static inline double snap0(double v, double eps) {
    return (std::fabs(v) < eps) ? 0.0 : v;
}

static inline double eval_leq( Halfspace& L, const std::vector<double>& x) {
    double s = 0.0;
    for (size_t i = 0; i < L.getCoords().size(); ++i) s += L.getCoords()[i] * x[i];
    return s - L.getScore();
}

static bool find_feasible(const std::vector<Halfspace>& I, int d, std::vector<double>& fp) {
    const double anchors[] = {0.1, 0.2, 0.3};
    for (double a : anchors){
        fp.assign(d, a);
        bool ok = true;
        for (auto& L : I) { if (eval_leq((Halfspace&)L, fp) > -1e-12){ ok = false; break; } }
        if (ok) return true;
    }
    fp.assign(d, 0.1);
    for (int t = 0; t < 8; ++t){
        bool ok = true;
        for (auto& L : I) { if (eval_leq((Halfspace&)L, fp) > -1e-12){ ok = false; break; } }
        if (ok) return true;
        for (int i = 0; i < d; i++) fp[i] *= 0.5;
    }
    return false;
}

/*  ===== main API =====
 * Computes the V-representation (intersection vertices) of halfspaces:
 *    - d: dimension (>=2)
 *    - ineq_in: rows of (a, c) meaning a·x <= c
 *    - add_unit_box: if true, also intersect with [0,1]^d
 *    - eps: numeric snap-to-zero for clean output
 * Returns: vector of vertices (each has size d).  No dedup/sorting (keep minimal).
 */
std::vector<std::vector<double>> vrep_halfspace_qhull(int d, const std::vector<Halfspace>& ineq_in, bool add_unit_box, double eps = 1e-12){
    std::vector<Halfspace> ineq = ineq_in;

    if (add_unit_box){
        for (int i = 0; i < d; i++){
            Halfspace ge0; ge0.getCoords().assign(d, 0.0); ge0.getCoords()[i] = -1.0; ge0.setScore(0.0); ineq.push_back(ge0);
            Halfspace le1; le1.getCoords().assign(d, 0.0); le1.getCoords()[i] =  1.0; le1.setScore(1.0); ineq.push_back(le1);
        }
    }

    std::vector<double> feas;
    if (!find_feasible(ineq, d, feas))
        return {};

    const int dimHS = d + 1;
    const int nHS   = (int)ineq.size();
    std::vector<coordT> hs; hs.reserve(nHS * dimHS);
    for (auto& L : ineq){
        for (int j = 0; j < d; j++) hs.push_back((coordT)L.getCoords()[j]);
        hs.push_back((coordT)(-L.getScore()));
    }

    std::ostringstream oss;
    oss << "qhull H";
    for (int i = 0; i < d; i++){ if (i) oss << ","; oss << feas[i]; }
    oss << " Qt";
    if (d >= 5) oss << " Qx";

    qhT qh_qh; qhT* qh = &qh_qh;
    qh_zero(qh, stderr);
    int exitcode = qh_new_qhull(qh, dimHS, nHS, hs.data(), qh_False, (char*)oss.str().c_str(), nullptr, stderr);

    std::vector<std::vector<double>> out;
    if (exitcode != 0) { qh_freeqhull(qh, !qh_ALL); int a=0,b=0; qh_memfreeshort(qh,&a,&b); return out; }

    facetT* facet;
    FORALLfacets {
        double off = facet->offset;
        if (std::fabs(off) < 1e-18) continue;
        std::vector<double> x(d);
        for (int i = 0; i < d; i++) x[i] = snap0((-facet->normal[i]/off + feas[i]), eps);
        out.push_back(std::move(x));
    }

    qh_freeqhull(qh, !qh_ALL);
    int curl=0, tot=0;
    qh_memfreeshort(qh,&curl,&tot);
    return out;
}
