#pragma once

#include "optimization.h"
#include "search_types.h"
#include <cstdint>
#include <vector>

class Facet;
class User;

double userScoreAtPoint(User* user, const std::vector<double>& point, int dim);

double averageScoreWidthInSlabUnits(const SubRegion& region,
                                    const std::vector<User*>& users,
                                    const Facet& facet,
                                    int dim);

void solveRelaxedSubRegion(
    SubRegion& region,
    const std::vector<Facet>& facetList,
    const std::vector<User*>& users,
    int dim,
    uint64_t& relaxedLpSolveCount,
    double& relaxedLpSolveSeconds,
    std::vector<UserHullCache>* hullCache = nullptr);

long long computeExactPointArfFromProductFile(
    const char* productFile,
    const std::vector<User*>& users,
    const std::vector<double>& point,
    int dim);
