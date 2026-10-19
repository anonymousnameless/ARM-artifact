#pragma once

#include "search_types.h"
#include <chrono>
#include <string>

struct AnytimeGapStatus {
  bool hasLowerBound = false;
  double lowerBound = 0.0;
  long long integerLowerBound = 0;
  long long certifiedGap = 0;
  double relativeGap = 0.0;
};

double bestAnytimeLowerBound(const SubRegion* current,
                             const SubRegionQueue* innerQueue,
                             const SubRegionQueue* globalQueue);

AnytimeGapStatus computeAnytimeGap(
    int incumbentArf,
    const SubRegion* current,
    const SubRegionQueue* innerQueue,
    const SubRegionQueue* globalQueue, long long rankOffset = 0);

bool anytimeGapTargetReached(
    int incumbentArf,
    double targetGapFraction,
    const SubRegion* current,
    const SubRegionQueue* innerQueue,
    const SubRegionQueue* globalQueue,
    AnytimeGapStatus* status = nullptr, long long rankOffset = 0);

void logAnytimeSnapshot(
    const std::string& event,
    int incumbentArf,
    const std::chrono::steady_clock::time_point& searchStart,
    const SubRegion* current,
    const SubRegionQueue* innerQueue,
    const SubRegionQueue* globalQueue, long long rankOffset = 0);
