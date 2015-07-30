#ifndef _validation_
#define _validation_
#include "Track.h"
#include "Hit.h"

class Validation {
public:
  virtual void fillSimHists(const TrackVec&, const std::vector<HitVec>&) {}
  virtual void fillCandidateHists(const TrackVec&) {}
  virtual void fillBuildHists(unsigned int, unsigned int, unsigned int) {}
  //  virtual void fillAssociationHists(TrackVec& evt_track_candidates, TrackVec& evt_sim_tracks, TrackVec& evt_assoc_tracks_RD, TrackVec& evt_assoc_tracks_SD) {}
  virtual void fillAssociationHists(const TrackVec& evt_track_candidates, const TrackVec& evt_sim_tracks, const std::vector<HitVec>& layerHits) {}
  virtual void fillFitStateHists(const TrackState&, const TrackState&) {}
  virtual void fillFitHitHists(unsigned int, const HitVec&, const MeasurementState&, const TrackState&, const TrackState&) {}
  virtual void fillFitTrackHists(const TrackState&, const TrackState&) {}
  virtual void saveHists() {}
  virtual void deleteHists() {}
};

#endif
