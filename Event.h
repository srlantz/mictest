#ifndef _event_
#define _event_

#define TEST_CLONE_ENGINE
//#define CLONE_ENGINE_SINGLE_THREAD
//#define BEST_OF_TEN

#include "Track.h"
#include "Validation.h"
#include "Geometry.h"

typedef std::pair<unsigned int,unsigned int> BinInfo;

// phi-eta partitioning map: vector of vector of vectors of std::pairs. 
// vec[nLayers][nEtaBins][nPhiBins]
typedef std::vector<std::vector<std::vector<BinInfo> > > BinInfoMap;

namespace Config {
  constexpr unsigned int nlayers_per_seed = 3;
  constexpr unsigned int maxCand = 10;
  constexpr float chi2Cut = 15.;
  constexpr float nSigma = 3.;
  constexpr float minDPhi = 0.;
};

class Event {
public:
  Event(const Geometry& g, Validation& v, int threads = 1);
  void Simulate(unsigned int nTracks);
  void Segment();
  void Seed();
  void Find();
  void Fit();

  const Geometry& geom_;
  Validation& validation_;
  std::vector<HitVec> layerHits_;
  //simTracks_, simHits_, simHitIdxs_ and initialHits_ work on sync, same index
  TrackVec simTracks_;
  std::vector<HitVec> simHits_;
  std::vector<std::vector<int> > simHitIdxs_;
  std::vector<HitVec>       initialHits_;
  std::vector<MCHitInfoVec> initialMCHitsInfo_;
  TrackVec seedTracks_, candidateTracks_;
  int threads_;

  BinInfoMap segmentMap_;
};

typedef std::vector<Event> EventVec;

#endif
