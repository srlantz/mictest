#include "buildtest.h"

#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"
#include "Event.h"
#include "Debug.h"

#include <cmath>
#include <iostream>
#include <mutex>

typedef std::pair<Track, TrackState> cand_t;
typedef TrackVec::const_iterator TrkIter;

#ifndef TBB
typedef std::vector<cand_t> candvec;
#else
#include "tbb/tbb.h"
typedef tbb::concurrent_vector<cand_t> candvec;
#endif
typedef candvec::const_iterator canditer;

void processCandidates(const Event& ev, const cand_t& cand, candvec& tmp_candidates, unsigned int ilay, bool debug);

inline float normalizedPhi(float phi) {
  return std::fmod(phi, M_PI);
}

static bool sortByHitsChi2(cand_t cand1, cand_t cand2)
{
  if (cand1.first.nHits()==cand2.first.nHits()) return cand1.first.chi2()<cand2.first.chi2();
  return cand1.first.nHits()>cand2.first.nHits();
}

static std::mutex evtlock;

void processSeed(Event& ev, const Track& seed, candvec& candidates, unsigned int ilay, const bool debug)
{
  auto& evt_track_candidates(ev.candidateTracks_);

  dprint("processing seed # " << seed.SimTrackID() << " par=" << seed.parameters() << " candidates=" << candidates.size());

  candvec tmp_candidates;
  tmp_candidates.reserve(3*candidates.size()/2);

  if (candidates.size() > 0) {
#ifdef TBB0 // turned off, loop is too small
    //loop over running candidates
    parallel_for( tbb::blocked_range<canditer>(candidates.begin(),candidates.end()), 
        [&](const tbb::blocked_range<canditer>& cands) {
      for (auto&& cand : cands) {
        processCandidates(ev, cand, tmp_candidates, ilay, debug);
      }
    }); //end of running candidates loop
#else
    //loop over running candidates
    for (auto&& cand : candidates) {
      processCandidates(ev, cand, tmp_candidates, ilay, debug);
    }
#endif

    ev.validation_.fillBuildHists(ilay, tmp_candidates.size(), candidates.size());

    if (tmp_candidates.size()>Config::maxCand) {
      dprint("huge size=" << tmp_candidates.size() << " keeping best "<< Config::maxCand << " only");
      std::partial_sort(tmp_candidates.begin(),tmp_candidates.begin()+(Config::maxCand-1),tmp_candidates.end(),sortByHitsChi2);
      tmp_candidates.resize(Config::maxCand); // thread local, so ok not thread safe
    } else if (tmp_candidates.size()==0) {
      dprint("no more candidates, saving best");
      // save the best candidate from the previous iteration and then swap in
      // the empty new candidate list; seed will be skipped on future iterations
      auto&& best = std::max_element(candidates.begin(),candidates.end(),sortByHitsChi2);
      std::lock_guard<std::mutex> evtguard(evtlock); // should be rare
      evt_track_candidates.push_back(best->first);
    }
    dprint("swapping with size=" << tmp_candidates.size());
    candidates.swap(tmp_candidates);
    tmp_candidates.clear();
  }
}

void buildTracks(Event& ev)
{
  auto& evt_track_candidates(ev.candidateTracks_);
  const auto& evt_lay_hits(ev.layerHits_);
  const auto& evt_seeds(ev.seedTracks_);
  bool debug(false);

  std::vector<candvec> track_candidates(evt_seeds.size());
  for (auto iseed = 0U; iseed < evt_seeds.size(); iseed++) {
    const auto& seed(evt_seeds[iseed]);
    track_candidates[iseed].push_back(cand_t(seed, seed.state()));
  }

  //loop over layers, starting from after the seed
  for (auto ilay = Config::nlayers_per_seed; ilay < evt_lay_hits.size(); ++ilay) {
    dprint("going to layer #" << ilay << " with N cands = " << track_candidates.size());

#ifdef TBB
    //process seeds
    parallel_for( tbb::blocked_range<size_t>(0, evt_seeds.size()), 
        [&](const tbb::blocked_range<size_t>& seediter) {
      for (auto iseed = seediter.begin(); iseed != seediter.end(); ++iseed) {
        const auto& seed(evt_seeds[iseed]);
        auto&& candidates(track_candidates[iseed]);
        processSeed(ev, seed, candidates, ilay, debug);
      }
    }); //end of process seeds loop
#else
    //process seeds
    for (auto iseed = 0U; iseed != evt_seeds.size(); ++iseed) {
      const auto& seed(evt_seeds[iseed]);
      auto&& candidates(track_candidates[iseed]);
      processSeed(ev, seed, candidates, ilay, debug);
    }
#endif
  } //end of layer loop

  //std::lock_guard<std::mutex> evtguard(evtlock);
  for (const auto& cand : track_candidates) {
    if (cand.size()>0) {
      // only save one track candidate per seed, one with lowest chi2
      //std::partial_sort(cand.begin(),cand.begin()+1,cand.end(),sortByHitsChi2);
      auto&& best = std::max_element(cand.begin(),cand.end(),sortByHitsChi2);
      evt_track_candidates.push_back(best->first);
    }
  }
}

void processCandidates(const Event& ev,
                       const cand_t& cand,
                       candvec& tmp_candidates,
                       unsigned int ilayer,
                       bool debug)
{
  const Track& tkcand = cand.first;
  const TrackState& updatedState = cand.second;
  const auto& evt_lay_hits(ev.layerHits_);
  const auto& evt_lay_phi_hit_idx(ev.lay_phi_hit_idx_);
    
  dprint("processing candidate with nHits=" << tkcand.nHits());
#ifdef LINEARINTERP
  TrackState propState = propagateHelixToR(updatedState,ev.geom_.Radius(ilayer));
#else
#ifdef TBB
#error "Invalid combination of options (thread safety)"
#endif
  TrackState propState = propagateHelixToLayer(updatedState,ilayer,&ev.geom_);
#endif // LINEARINTERP
#ifdef CHECKSTATEVALID
  if (!propState.valid) {
    return;
  }
#endif
  const float predx = propState.parameters.At(0);
  const float predy = propState.parameters.At(1);
#ifdef DEBUG
  if (debug) {
    const float predz = propState.parameters.At(2);
    std::cout << "propState at hit#" << ilayer << " r/phi/z : " << sqrt(pow(predx,2)+pow(predy,2)) << " "
              << std::atan2(predy,predx) << " " << predz << std::endl;
    dumpMatrix(propState.errors);
  }
#endif
  
  const float phi = std::atan2(predy,predx);

  const float px2py2 = predx*predx+predy*predy;
  const float dphidx = -predy/px2py2;
  const float dphidy =  predx/px2py2;
  const float dphi2  = dphidx*dphidx*(propState.errors.At(0,0)) +
                       dphidy*dphidy*(propState.errors.At(1,1)) +
                     2*dphidy*dphidx*(propState.errors.At(0,1));
  const float dphi   =  sqrt(std::abs(dphi2));//how come I get negative squared errors sometimes?
  
  const float nSigmaDphi = std::min(std::max(Config::nSigma*dphi, (float) Config::minDPhi), (float) M_PI);
  const float dphiMinus = normalizedPhi(phi-nSigmaDphi);
  const float dphiPlus  = normalizedPhi(phi+nSigmaDphi);
  
  const unsigned int binMinus = getPhiPartition(dphiMinus);
  const unsigned int binPlus  = getPhiPartition(dphiPlus);
  
  dprint("phi: " << phi << " binMinus: " << binMinus << " binPlus: " << binPlus << " dphi2: " << dphi2);
  
  const BinInfo& binInfoMinus = evt_lay_phi_hit_idx[ilayer][int(binMinus)];
  const BinInfo& binInfoPlus  = evt_lay_phi_hit_idx[ilayer][int(binPlus)];
 
  unsigned int firstIndex = binInfoMinus.first;
  unsigned int maxIndex   = binInfoPlus.first+binInfoPlus.second;
  unsigned int lastIndex  = -1;
  unsigned int totalSize  = evt_lay_hits[ilayer].size(); 

  // Branch here from wrapping
  if (binMinus<=binPlus){
    lastIndex = maxIndex;
  } else { // loop wrap around end of array for binMinus > binPlus, for dPhiMinus < 0 or dPhiPlus > 0 at initialization
    lastIndex = totalSize+maxIndex;
  }

  dprint("Count: " << lastIndex-firstIndex << "/" << totalSize << " firstIndex: " << firstIndex 
                   << " maxIndex: " << maxIndex << " lastIndex: " << lastIndex);

#ifdef LINEARINTERP
  const float minR = ev.geom_.Radius(ilayer);
  float maxR = minR;
  for (unsigned int ihit=firstIndex;ihit<lastIndex;++ihit) {//loop over hits on layer (consider only hits from partition)
    const float candR = evt_lay_hits[ilayer][ihit % totalSize].r();
    if (candR > maxR) maxR = candR;
  }
  const float deltaR = maxR - minR;

  dprint("min, max: " << minR << ", " << maxR);
  const TrackState propStateMin = propState;
  const TrackState propStateMax = propagateHelixToR(updatedState,maxR);
#ifdef CHECKSTATEVALID
  if (!propStateMax.valid) {
    return;
  }
#endif
#endif
    
#ifdef TBB0 // turned off, loop is too small
  //loop over hits on layer (consider only hits from partition)
  parallel_for( tbb::blocked_range<size_t>(firstIndex, lastIndex), 
      [&](const tbb::blocked_range<size_t>& ihits) {
    for (auto ihit = ihits.begin(); ihit != ihits.end(); ++ihit) {
#else
    for (auto ihit = firstIndex; ihit != lastIndex; ++ihit) {
#endif
      const auto& hitCand = evt_lay_hits[ilayer][ihit % totalSize];
    
      const MeasurementState& hitMeas = hitCand.measurementState();

#ifdef LINEARINTERP
      const float ratio = (hitCand.r() - minR)/deltaR;
      propState.parameters = (1.0-ratio)*propStateMin.parameters + ratio*propStateMax.parameters;
      dprint(std::endl << ratio << std::endl << propStateMin.parameters << std::endl << propState.parameters << std::endl
                       << propStateMax.parameters << std::endl << propStateMax.parameters - propStateMin.parameters
                       << std::endl << std::endl << hitMeas.parameters << std::endl);
#endif
      const float chi2 = computeChi2(propState,hitMeas);
    
#ifdef DEBUG
      if (debug) {
        const float hitx = hitCand.position()[0];
        const float hity = hitCand.position()[1];
        const float hitz = hitCand.position()[2];
        std::cout << "consider hit r/phi/z : " << sqrt(pow(hitx,2)+pow(hity,2)) << " "
                  << std::atan2(hity,hitx) << " " << hitz << " chi2=" << chi2 << std::endl;
      }
#endif
    
      if ((chi2<Config::chi2Cut)&&(chi2>0.)) {//fixme 
        dprint("found hit with index: " << ihit << " chi2=" << chi2);
        const TrackState tmpUpdatedState = updateParameters(propState, hitMeas);
        Track tmpCand = tkcand.clone();
        tmpCand.addHit(hitCand,chi2);
        tmp_candidates.push_back(cand_t(tmpCand,tmpUpdatedState));
      }
    }
#ifdef TBB0
  }); //end of consider hits on layer loop
#endif

  //add also the candidate for no hit found
  if (tkcand.nHits()==ilayer) { //only if this is the first missing hit
    dprint("adding candidate with no hit");
    tmp_candidates.push_back(cand_t(tkcand,propState));
  }
}
