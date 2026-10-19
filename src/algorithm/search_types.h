#pragma once

#include "optimization.h"
#include "halfspace.h"
#include <queue>
#include <vector>

struct UserGap {
  int userIndex;
  double gap;
  int splitSlab;
};

// One search node inside a fixed facet. The node is defined by per-user slab
// ranges plus any halfspaces introduced by branching decisions.
struct SubRegion {
  int facetIdx;
  double relaxedARF;
  std::vector<int> minBorders;
  std::vector<int> maxBorders;
  std::vector<double> nodeScoreMins;
  std::vector<double> nodeScoreMaxs;
  std::vector<double> relaxedRanks;
  std::vector<double> wLP;
  LpBasisWarmStart relaxedBasis;
  std::vector<Halfspace> splitHalfspaces;
};

struct RelaxedArfGreater {
  bool operator()(const SubRegion& a, const SubRegion& b) const {
    return a.relaxedARF > b.relaxedARF;
  }
};

using SubRegionQueue =
    std::priority_queue<SubRegion, std::vector<SubRegion>, RelaxedArfGreater>;
