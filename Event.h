#ifndef _event_
#define _event_

#include "Track.h"
#include "Validation.h"
#include "Geometry.h"

typedef std::pair<unsigned int,unsigned int> BinInfo;

// phi-eta partitioning map: vector of vector of vectors of std::pairs. 
// vec[nLayers][nEtaBins][nPhiBins]
typedef std::vector<std::vector<std::vector<BinInfo> > > BinInfoMap;

namespace Config {
  static constexpr const unsigned int nlayers_per_seed = 3;
  static constexpr const unsigned int maxCand = 10;
  static constexpr const float chi2Cut = 15.;
  static constexpr const float nSigma = 3.;
  static constexpr const float minDPhi = 0.;
};

class Event {
public:
  Event(const Geometry& g, Validation& v, int threads = 1);
  void Simulate(unsigned int nTracks);
  void Segment();
  void Seed();
  void Find();
  void Fit();

  SimTkIDInfo SimTrackIDInfo(const Track& trk) const;

  // from mchitid to hitid
  const HitID HitIDFromMCID(uint32_t id) const
  {
    return mcHits_[id].hitID();
  }
  const HitID HitIDFromMCID(HitID id) const
  {
    assert(HitID::MCLayerID == id.layer_);
    return HitIDFromMCID(id.index_);
  }

  const Hit& HitFromID(HitID id) const
  {
    assert(layerHits_.size() > id.layer_);
    assert(layerHits_[id.layer_].size() > id.index_);
    return layerHits_[id.layer_][id.index_];
  }

  const Hit& HitFromMCID(uint32_t id) const
  {
    return HitFromID(HitIDFromMCID(id));
  }

  const Hit& HitFromMCID(HitID id) const
  {
    return HitFromID(HitIDFromMCID(id));
  }

  const MCHit& MCHitFromHitID(HitID id) const
  {
    assert(HitID::MCLayerID != id.layer_);
    return mcHits_[HitFromID(id).mcHitID()];
  }

  const MCHit& MCHitFromHit(const Hit& hit) const
  {
    return mcHits_[hit.mcHitID()];
  }

  const MCHit& MCHitFromMCID(HitID id) const
  {
    assert(HitID::MCLayerID == id.layer_);
    return mcHits_[id.index_];
  }

  const Geometry& geom_;
  Validation& validation_;
  std::vector<HitVec> layerHits_;
  MCHitVec mcHits_;
  std::vector<float> minR_, maxR_;
  TrackVec simTracks_, seedTracks_, candidateTracks_;
  int threads_;

  BinInfoMap segmentMap_;
};

typedef std::vector<Event> EventVec;

#endif
