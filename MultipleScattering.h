#include "Track.h"
#include "Event.h"
#include "Geometry.h"

void applyMultipleScattering(const Geometry& geom, const TrackState& propState, const float z1, const float z2, const float phismear, SVector6& deltaPar);
