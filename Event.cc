#include "Event.h"
#include "Simulation.h"
#include "KalmanUtils.h"
#include "buildtest.h"
#include "fittest.h"
#include "Debug.h"

#ifdef TBB
#include "tbb/tbb.h"
#endif

// find the simtrack that provided the most hits
SimTkIDInfo Event::SimTrackIDInfo(const Track& trk) const
{
  const auto& hitids = trk.hitIDs();

  if (hitids.size() == 0){
    return SimTkIDInfo(0,0);
  }

  std::vector<unsigned int> mctrack;
  for (auto&& ihit : hitids){
    auto&& mchit = MCHitFromHitID(ihit);
    mctrack.push_back(mchit.mcTrackID());
  }
  std::sort(mctrack.begin(), mctrack.end()); // ensures all elements are checked properly

  auto mtrk(mctrack[0]), mcount(0U), m(mctrack[0]), c(0U);

  for (auto i : mctrack) {
    if (i == m) { ++c; } else { c = 0; }
    if (c >= mcount) { mtrk = m; mcount = c; }
    m = i;
  }

  return SimTkIDInfo(mtrk,mcount);
}

static bool sortByPhi(const Hit& hit1, const Hit& hit2)
{
  return hit1.phi()<hit2.phi();
}

static bool tracksByPhi(const Track& t1, const Track& t2)
{
  return t1.posPhi()<t2.posPhi();
}

#ifdef ETASEG
const float etaDet = 2.0;
static bool sortByEta(const Hit& hit1, const Hit& hit2){
  return hit1.eta()<hit2.eta();
}

// within a layer with a "reasonable" geometry, ordering by Z is the same as eta
static bool sortByZ(const Hit& hit1, const Hit& hit2){
  return hit1.z()<hit2.z();
}
#endif

Event::Event(const Geometry& g, Validation& v, int threads) : geom_(g), validation_(v), threads_(threads)
{
  MCHit::mcHitIDCounter_ = 0;
  layerHits_.resize(geom_.CountLayers());
  minR_.resize(geom_.CountLayers());
  maxR_.resize(geom_.CountLayers());
  segmentMap_.resize(geom_.CountLayers());
  for (auto ilay = 0U; ilay < geom_.CountLayers(); ++ilay) {
    minR_[ilay] = maxR_[ilay] = geom_.Radius(ilay);
  }
}

void Event::Simulate(unsigned int nTracks)
{
  simTracks_.resize(nTracks);
  for (auto&& l : layerHits_) {
    l.resize(nTracks);
  }
  mcHits_.resize(geom_.CountLayers()*nTracks);

#if defined(TBB)
  parallel_for( tbb::blocked_range<size_t>(0, nTracks, 100), 
      [&](const tbb::blocked_range<size_t>& itracks) {

    const Geometry tmpgeom(geom_.clone()); // USolids isn't thread safe??
    for (auto itrack = itracks.begin(); itrack != itracks.end(); ++itrack) {
#else
    const Geometry& tmpgeom(geom_);
    for (unsigned int itrack=0; itrack<nTracks; ++itrack) {
#endif
      //create the simulated track
      SVector3 pos;
      SVector3 mom;
      SMatrixSym66 covtrk;
      HitVec hits;
      MCHitVec mchits;
      // unsigned int starting_layer  = 0; --> for displaced tracks, may want to consider running a separate Simulate() block with extra parameters

      int q=0;//set it in setup function
      float pt = 0.5+g_unif(g_gen)*9.5;//this input, 0.5<pt<10 GeV (below ~0.5 GeV does not make 10 layers)
      setupTrackByToyMC(itrack,pt,pos,mom,covtrk,mchits,hits,q,tmpgeom);
      for (auto ihit = 0U; ihit < hits.size(); ++ihit) {
        const auto& hit = hits[ihit];
        auto ilay = mchits[ihit].layer();
        layerHits_[ilay][itrack] = hit;
        auto r = hit.r();
        if (r < minR_[ilay]) minR_[ilay] = r;
        if (r > maxR_[ilay]) maxR_[ilay] = r;
      }
      HitIDVec hitids;
      for (auto&& hit : mchits) {
        mcHits_[hit.mcHitID()] = hit;
        hitids.push_back(hit.mcHitID());
      }
      simTracks_[itrack] = Track(q,pos,mom,covtrk,hitids,0.0);
    }
#if defined(TBB)
  });
#endif
  validation_.fillSimHists(simTracks_, mcHits_);
}

template<typename Container, typename Counter>
void partition_sort(Container& c, Counter& count, unsigned int maxKey) {
  std::vector<Container> buckets (maxKey);

  for (auto&& ele : c) {
    const auto bin = ele.binIndex();
    buckets[bin].push_back(ele);
    ++count[bin];
  }
  auto i(0U);
  for (auto&& b : buckets) {
    for (auto&& e : b) {
      c[i++] = e;
    }
  }
}

void Event::Segment()
{
  for (unsigned int ilayer=0; ilayer<layerHits_.size(); ++ilayer) {

    std::vector<int> hitsinseg(Config::nPhiPart*Config::nEtaPart);
    partition_sort(layerHits_[ilayer], hitsinseg, Config::nPhiPart*Config::nEtaPart);

    BinInfo cbin(0,0);
    segmentMap_[ilayer].resize(Config::nEtaPart);
    for (auto ieta = 0U; ieta < Config::nEtaPart; ++ieta) {
      segmentMap_[ilayer][ieta].resize(Config::nPhiPart);
      for (auto iphi = 0U; iphi < Config::nPhiPart; ++iphi) {
        cbin.second = hitsinseg[binNumber(ieta, iphi)];
        segmentMap_[ilayer][ieta][iphi] = cbin;
        cbin.first += cbin.second;
      }
    }
    for (auto ihit = 0U; ihit < layerHits_[ilayer].size(); ++ihit) {
      auto mcid = layerHits_[ilayer][ihit].mcHitID();
      mcHits_[mcid].setHitID(HitID(ilayer, ihit));
    }
  }
}

void Event::Seed()
{
#define SIMSEEDS
#ifdef SIMSEEDS
  //create seeds (from sim tracks for now)
  for (unsigned int itrack=0;itrack<simTracks_.size();++itrack) {
    const Track& trk = simTracks_[itrack];
    const auto& simhitids = trk.hitIDs();
    TrackState updatedState = trk.state();
    HitIDVec seedhits;

    for (auto ilayer=0U;ilayer<Config::nlayers_per_seed;++ilayer) {//seeds have first three layers as seeds
      const auto& seed_hit = HitFromMCID(simhitids[ilayer]);
      TrackState propState = propagateHelixToR(updatedState,seed_hit.r());
      MeasurementState measState = seed_hit.measurementState();
      updatedState = updateParameters(propState, measState);
      seedhits.push_back(HitIDFromMCID(simhitids[ilayer]));
    }
    Track seed(updatedState,seedhits,0.);//fixme chi2
    seedTracks_.push_back(seed);
  }
#else

  // follow CMSSW -> start with hits in 2nd layer, build seed by then projecting back to beamspot.
  // build seed pairs... then project to third layer.

  // p=qBR => R is rad of curvature (in meters), q = 0.2997 in natural units, B = 3.8T, p = min pt = 0.5 GeV. 
  // so R = (0.5/0.2997...*3.8) * 100 -> max rad. of curv. in cm 

  const float curvRad = 4.38900125;
  const float d0 = 0.1; // 1 mm x,y beamspot from simulation
  const float dZ = 2.0;

  for (unsigned int ihit=0;ihit<layerHits_[1].size();++ihit) { // 1 = second layer
    float outerrad  = layerHits_[1][ihit].r();
    float ouerphi   = layerHits_[1][ihit].phi();

    unsigned int mcID = layerHits_[1][ihit].mcTrackID();
    Hit innerHit = simTracks_[mcID].hitsVector()[0];
    float innerrad  = innerHit.r();
    innerrad = 4.0;
  
    std::cout << "Diff: " << innerrad - 4.0 <<std::endl;

    float innerPhiPlusCentral  = outerphi-acos(outerrad/(2.*curvRad))+acos(innerrad/(2.*curvRad));
    float innerPhiMinusCentral = outerphi-acos(outerrad/(-2.*curvRad))+acos(-innerrad/(2.*curvRad));

    // for d0 displacements

    float alphaPlus = acos(-((d0*d0)+(outerrad*outerrad)-2.*(d0*curvRad))/((2.*outerrad)*(d0-curvRad)));
    float betaPlus  = acos(-((d0*d0)+(innerrad*innerrad)-2.*(d0*curvRad))/((2.*innerrad)*(d0-curvRad)));

    float alphaMinus = acos(-((d0*d0)+(outerrad*outerrad)+2.*(d0*curvRad))/((2.*outerrad)*(d0+curvRad)));
    float betaMinus  = acos(-((d0*d0)+(innerrad*innerrad)+2.*(d0*curvRad))/((2.*innerrad)*(d0+curvRad)));

    float innerPhiPlus = outerphi-alphaPlus+betaPlus;
    float innerPhiMinus = outerphi-alphaMinus+betaMinus;

    float innerZPlus  = (innerrad/outerrad)*(outerhitz-dZ)+dZ;
    float innerZMinus = (innerrad/outerrad)*(outerhitz+dZ)-dZ;
    float centralZ    = (innerrad/outerrad)*outerhitz;

    printf("ihit: %1u \n   iphi: %5f ophi: %5f \n   iphiP: %5f iphiM: %5f \n   iphiPC: %5f iphiMC: %5f \n",
           ihit,
           innerHit.phi(), outerphi,
           innerPhiPlus,innerPhiMinus,
           innerPhiPlusCentral,innerPhiMinusCentral
           );
    printf("   innerZ: %5f iZM: %5f iZP: %5f cZ: %5f \n",
           innerhitz,innerZMinus,innerZPlus,centralZ
           );

    unsigned int phibinPlus  = getPhiPartition(innerPhiPlus);
    unsigned int phibinMinus = getPhiPartition(innerPhiMinus);

  }
#endif
  std::sort(seedTracks_.begin(), seedTracks_.end(), tracksByPhi);
}

void Event::Find()
{
  //buildTracksBySeeds(*this);
  buildTracksByLayers(*this);

  // From CHEP-2015
  // buildTestSerial(*this, Config::nlayers_per_seed, Config::maxCand, Config::chi2Cut, Config::nSigma, Config::minDPhi);

  validation_.fillAssociationHists(*this);
  validation_.fillCandidateHists(candidateTracks_);
}

void Event::Fit()
{
#ifdef ENDTOEND
  runFittingTest(*this, candidateTracks_);
#else
  runFittingTest(*this, simTracks_);
#endif
}

