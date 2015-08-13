#ifndef _validation_
#define _validation_
#include "Track.h"
#include "Hit.h"

class Event;

class Validation {
public:
  virtual void fillSimHists(const TrackVec&, const MCHitVec&) {}
  virtual void fillCandidateHists(const TrackVec&) {}
  virtual void fillBuildHists(unsigned int, unsigned int, unsigned int) {}
  //  virtual void fillAssociationHists(TrackVec& evt_track_candidates, TrackVec& evt_sim_tracks, TrackVec& evt_assoc_tracks_RD, TrackVec& evt_assoc_tracks_SD) {}
  virtual void fillAssociationHists(const Event& ev) {}
  virtual void fillFitStateHists(const TrackState&, const TrackState&) {}
  virtual void fillFitHitHists(const MeasurementState&, const MeasurementState&, const TrackState&, const TrackState&) {}
  virtual void fillFitTrackHists(const TrackState&, const TrackState&) {}
  virtual void saveHists() {}
  virtual void deleteHists() {}
};

#endif
