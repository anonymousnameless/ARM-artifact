#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "../data/halfspace.h"

// Forward declarations
struct User;

// Hull representation (standard convex hull)
struct HullVRep {
    int d;
    std::vector<int> indices;
    std::vector<double> coords;
};

// Facet representation (for WCA/UBC)
struct FacetVRep {
    std::vector<int> vertex_indices;           // Indices of vertices in the global list
    std::vector<double> representative;        // Centroid of the facet

    // The specific halfspace inequality that generated this facet.
    Halfspace defining_halfspace;

    // Index in the original input vector.
    // Returns -1 if this facet was created by a Box Constraint (e.g., x >= 0).
    int original_index;

    // Helpers (derived from defining_halfspace)
    std::vector<double> hyperplane_normal;
    double hyperplane_offset;

    // --- Helper Functions ---

    // Get the actual bounded vertices coordinates for this facet
    std::vector<std::vector<double>> get_bounded_vertices(const std::vector<std::vector<double>>& all_vertices) const {
        std::vector<std::vector<double>> facet_vertices;
        facet_vertices.reserve(vertex_indices.size());
        for (int idx : vertex_indices) {
            facet_vertices.push_back(all_vertices[idx]);
        }
        return facet_vertices;
    }

    // Get bounded segment length (for edges in 2D/3D)
    double get_segment_length(const std::vector<std::vector<double>>& all_vertices) const {
        if (vertex_indices.size() != 2) return 0.0;
        const auto& v1 = all_vertices[vertex_indices[0]];
        const auto& v2 = all_vertices[vertex_indices[1]];
        double dist_sq = 0.0;
        for (size_t i = 0; i < v1.size(); ++i) {
            double diff = v2[i] - v1[i];
            dist_sq += diff * diff;
        }
        return std::sqrt(dist_sq);
    }

    // Simple check if a point lies on the hyperplane of this facet
    bool is_point_on_plane(const std::vector<double>& point, double eps = 1e-9) const {
        double score = 0.0;
        for (size_t i = 0; i < hyperplane_normal.size(); ++i) {
            score += hyperplane_normal[i] * point[i];
        }
        return std::fabs(score - hyperplane_offset) < eps;
    }
};

// WCA Result with vertices, facets, and the full constraint set
struct WCAResult {
    std::vector<std::vector<double>> vertices;
    std::vector<FacetVRep> facets;

    // Contains ALL inequalities used to build the shape (Input + Box constraints).
    std::vector<Halfspace> all_inequalities;
};

// Box constraint options
enum class BoxConstraints {
    NONE,           // No box constraints
    LOWER_ONLY,     // Only x_i >= 0
    FULL_UNIT_BOX   // Both 0 <= x_i <= 1
};

// ===== API =====

WCAResult get_ubc_vrep(
    int d,
    std::vector<Halfspace>& ineq_in,
    BoxConstraints box_type = BoxConstraints::LOWER_ONLY,
    double eps = 1e-12
);

HullVRep convex_hull_vrep_nd(int d, const std::vector<double>& pts);
HullVRep convex_hull_vrep_from_users(std::vector<User*>& users);
void hull_vertices_users(const HullVRep& H, std::vector<User*>& users, std::vector<User*>& hullUsers);