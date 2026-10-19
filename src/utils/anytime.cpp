#include "anytime.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

double bestAnytimeLowerBound(const SubRegion* current,
                             const SubRegionQueue* innerQueue,
                             const SubRegionQueue* globalQueue) {
  double lowerBound = std::numeric_limits<double>::infinity();
  if (current != nullptr && std::isfinite(current->relaxedARF)) {
    lowerBound = std::min(lowerBound, current->relaxedARF);
  }
  if (innerQueue != nullptr && !innerQueue->empty() &&
      std::isfinite(innerQueue->top().relaxedARF)) {
    lowerBound = std::min(lowerBound, innerQueue->top().relaxedARF);
  }
  if (globalQueue != nullptr && !globalQueue->empty() &&
      std::isfinite(globalQueue->top().relaxedARF)) {
    lowerBound = std::min(lowerBound, globalQueue->top().relaxedARF);
  }
  return lowerBound;
}

AnytimeGapStatus computeAnytimeGap(
    int incumbentArf,
    const SubRegion* current,
    const SubRegionQueue* innerQueue,
    const SubRegionQueue* globalQueue, long long rankOffset) {
  AnytimeGapStatus status;
  if (incumbentArf == INT_MAX) {
    status.relativeGap = std::numeric_limits<double>::infinity();
    return status;
  }

  const double lowerBound =
      bestAnytimeLowerBound(current, innerQueue, globalQueue);
  if (!std::isfinite(lowerBound)) {
    status.relativeGap = std::numeric_limits<double>::infinity();
    return status;
  }

  status.hasLowerBound = true;
  status.lowerBound = lowerBound + rankOffset;
  status.integerLowerBound = std::min(
      (long long)incumbentArf + rankOffset,
      std::max(0LL, (long long)std::ceil(lowerBound - 1e-6)) + rankOffset);
  status.certifiedGap =
      std::max(0LL, (long long)incumbentArf + rankOffset - status.integerLowerBound);

  if (status.integerLowerBound > 0) {
    status.relativeGap =
        (double)status.certifiedGap / (double)status.integerLowerBound;
  } else {
    status.relativeGap =
        (status.certifiedGap == 0)
            ? 0.0
            : std::numeric_limits<double>::infinity();
  }
  return status;
}

bool anytimeGapTargetReached(
    int incumbentArf,
    double targetGapFraction,
    const SubRegion* current,
    const SubRegionQueue* innerQueue,
    const SubRegionQueue* globalQueue,
    AnytimeGapStatus* status, long long rankOffset) {
  AnytimeGapStatus computed =
      computeAnytimeGap(incumbentArf, current, innerQueue, globalQueue, rankOffset);
  if (status != nullptr) {
    *status = computed;
  }
  return computed.hasLowerBound &&
         std::isfinite(computed.relativeGap) &&
         std::isfinite(targetGapFraction) &&
         targetGapFraction >= 0.0 &&
         computed.relativeGap <= targetGapFraction + 1e-12;
}

void logAnytimeSnapshot(
    const std::string& event,
    int incumbentArf,
    const std::chrono::steady_clock::time_point& searchStart,
    const SubRegion* current,
    const SubRegionQueue* innerQueue,
    const SubRegionQueue* globalQueue, long long rankOffset) {
  const double elapsed =
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    searchStart)
          .count();
  const AnytimeGapStatus gapStatus =
      computeAnytimeGap(incumbentArf, current, innerQueue, globalQueue, rankOffset);

  std::ios::fmtflags oldFlags = std::cout.flags();
  std::streamsize oldPrecision = std::cout.precision();
  std::cout << std::fixed << std::setprecision(3)
            << "[ANYTIME] elapsed=" << elapsed
            << " incumbent=" << (long long)incumbentArf + rankOffset;
  if (gapStatus.hasLowerBound) {
    std::cout << " lp_bound=" << std::setprecision(6) << gapStatus.lowerBound
              << " int_bound=" << gapStatus.integerLowerBound
              << " certified_gap=" << gapStatus.certifiedGap;
  } else {
    std::cout << " lp_bound=NA int_bound=NA certified_gap=NA";
  }
  std::cout << " event=" << event << std::endl;
  std::cout.flags(oldFlags);
  std::cout.precision(oldPrecision);
}
