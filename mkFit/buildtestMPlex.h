#ifndef _buildtest_mplex_
#define _buildtest_mplex_

#include "Event.h"
#include "Track.h"

#ifdef TEST_CLONE_ENGINE
#include "CandCloner.h"
#endif

double runBuildingTest(Event& ev);
double runBuildingTestBestHit(Event& ev);

void buildTestParallel(std::vector<Track>& evt_seeds,std::vector<Track>& evt_track_candidates,
		       std::vector<std::vector<Hit> >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,
		       const int& nhits_per_seed,const unsigned int& maxCand,const float& chi2Cut,const float& nSigma,const float& minDPhi);

void processCandidates(std::pair<Track, TrackState>& cand,std::vector<std::pair<Track, TrackState> >& tmp_candidates,
		       unsigned int ilay,std::vector<std::vector<Hit> >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx,
		       const int& nhits_per_seed,const unsigned int& maxCand,const float& chi2Cut,const float& nSigma,const float& minDPhi);

#ifdef TEST_CLONE_ENGINE
double runBuildingTestPlex(Event& ev, CandCloner& cloner);
#else
double runBuildingTestPlex(Event& ev);
#endif
double runBuildingTestPlexOld(Event& ev);
double runBuildingTestPlexBestHit(Event& ev);

#endif
