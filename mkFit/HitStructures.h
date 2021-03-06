#ifndef HitStructures_H
#define HitStructures_H

#include "Config.h"
#include "Hit.h"
#include "Track.h"
//#define DEBUG
#include "Debug.h"

// for each layer
//   Config::nEtaBin vectors of hits, resized to large enough N
//   filled with corresponding hits
//   sorted in phi
//   and with a corresponding BinInfoSomething
//
// My gut feeling is that we could have them all in one place and just store
// indices into small enough eta-phi bins (oh, and eta is a horrible separation variable)
// 10*50000*64/1024/1024
//    30.51757812500000000000
// Or make them smaller ... if we could use short indices that might make a difference.
//
// Need to add eta/phi ... or something significantly better to the Hit struct
//
// At least ... fast eta/phi estimators.
//
// Oh, and we should use radix sort. Could one vectorize this?

// Need a good "array of pods" class with aligned alloc and automatic growth.
// For now just implement the no-resize / no-destroy basics in the BoH.

typedef std::pair<int, int> PhiBinInfo_t;

//==============================================================================

inline bool sortHitsByPhiMT(const Hit& h1, const Hit& h2)
{
  return std::atan2(h1.position()[1],h1.position()[0])<std::atan2(h2.position()[1],h2.position()[0]);
}

inline bool sortTrksByPhiMT(const Track& t1, const Track& t2)
{
  return t1.momPhi() < t2.momPhi();
}

//==============================================================================

class BunchOfHits
{
public:
  // This eventually becomes a vector of pointers / refs / multi indices
  //std::vector<Hit>          m_hits;
  Hit                      *m_hits;
  std::vector<PhiBinInfo_t> m_phi_bin_infos;

  int     m_real_size;
  int     m_fill_index;
  int     m_fill_index_old;

public:
  BunchOfHits() :
    //m_hits          (Config::maxHitsPerBunch),
    m_phi_bin_infos (Config::nPhiPart),
    m_real_size     (Config::maxHitsPerBunch),
    m_fill_index    (0),
    m_fill_index_old(0)
  {
    m_hits = (Hit*) _mm_malloc(sizeof(Hit)*Config::maxHitsPerBunch, 64);
    Reset();
  }

  ~BunchOfHits()
  {
    _mm_free(m_hits);
  }
  
  void Reset();

  void InsertHit(const Hit& hit)
  {
    assert (m_fill_index < m_real_size); // or something

    m_hits[m_fill_index] = hit;
    ++m_fill_index;
  }

  void SortByPhiBuildPhiBins();
};

//==============================================================================

class LayerOfHits
{
public:
  std::vector<BunchOfHits> m_bunches_of_hits;

public:
  LayerOfHits() :
    m_bunches_of_hits(Config::nEtaBin)
  {}

  void Reset()
  {
    for (auto &i : m_bunches_of_hits)
    {
      i.Reset();
    }
  }

  void InsertHit(const Hit& hit)
  {
    int b1, b2;
    int cnt = getBothEtaBins(hit.eta(), b1, b2);

    if (b1 != -1) m_bunches_of_hits[b1].InsertHit(hit);
    if (b2 != -1) m_bunches_of_hits[b2].InsertHit(hit);
  }

  void SortByPhiBuildPhiBins()
  {
    for (auto &i : m_bunches_of_hits)
    {
      i.SortByPhiBuildPhiBins();
    }
  }
};

//==============================================================================

class EventOfHits
{
public:
  std::vector<LayerOfHits> m_layers_of_hits;
  int                      m_n_layers;

public:
  EventOfHits(int n_layers) :
    m_layers_of_hits(n_layers),
    m_n_layers(n_layers)
  {}

  void Reset()
  {
    for (auto &i : m_layers_of_hits)
    {
      i.Reset();
    }
  }

  void InsertHit(const Hit& hit, int layer)
  {
    m_layers_of_hits[layer].InsertHit(hit);
  }

  void SortByPhiBuildPhiBins()
  {
    for (auto &i : m_layers_of_hits)
    {
      i.SortByPhiBuildPhiBins();
    }
  }
};


//==============================================================================
//==============================================================================

// This should actually be a BunchOfCandidates that share common hit vector.
// At the moment this is an EtaBin ...

class EtaBinOfCandidates
{
public:
  std::vector<Track> m_candidates;

  int     m_real_size;
  int     m_fill_index;

public:
  EtaBinOfCandidates() :
    m_candidates (Config::maxCandsPerEtaBin),
    m_real_size  (Config::maxCandsPerEtaBin),
    m_fill_index (0)
  {}

  void Reset()
  {
    m_fill_index = 0;
  }

  void InsertTrack(const Track& track)
  {
    assert (m_fill_index < m_real_size); // or something

    m_candidates[m_fill_index] = track;
    ++m_fill_index;
  }

  void SortByPhi()
  {
    std::sort(m_candidates.begin(), m_candidates.begin() + m_fill_index, sortTrksByPhiMT);
  }
};

class EventOfCandidates
{
public:
  std::vector<EtaBinOfCandidates> m_etabins_of_candidates;

public:
  EventOfCandidates() :
    m_etabins_of_candidates(Config::nEtaBin)
  {}

  void Reset()
  {
    for (auto &i : m_etabins_of_candidates)
    {
      i.Reset();
    }
  }

  void InsertCandidate(const Track& track)
  {
    int bin = getEtaBinExtendedEdge(track.posEta());
    // XXXX MT Had to add back this conditional for best-hit (bad seeds)
    // Note also the ExtendedEdge above, this practically removes bin = -1
    // occurence and improves efficiency.
    if (bin != -1)
    {
      m_etabins_of_candidates[bin].InsertTrack(track);
    }
  }

  void SortByPhi()
  {
    for (auto &i : m_etabins_of_candidates)
    {
      i.SortByPhi();
    }
  }
};



//-------------------------------------------------------
// for combinatorial version, switch to vector of vectors
//-------------------------------------------------------

class EtaBinOfCombCandidates
{
public:
  std::vector<std::vector<Track> > m_candidates;

  //these refer to seeds
  int     m_real_size;
  int     m_fill_index;

public:
  EtaBinOfCombCandidates() :
    m_candidates (Config::maxCandsPerEtaBin / Config::maxCandsPerSeed),
    m_real_size  (Config::maxCandsPerEtaBin / Config::maxCandsPerSeed),
    m_fill_index (0)
  {
    for (int s=0;s<m_real_size;++s)
    {
      m_candidates[s].reserve(Config::maxCandsPerSeed);//we should never exceed this
    }
  }

  void Reset()
  {
    m_fill_index = 0;
    for (int s=0;s<m_real_size;++s)
    {
      m_candidates[s].clear();
    }
  }

  void InsertSeed(const Track& seed)
  {
    assert (m_fill_index < m_real_size); // or something

    m_candidates[m_fill_index].push_back(seed);
    ++m_fill_index;
  }

  void InsertTrack(const Track& track, int seed_index)
  {
    assert (seed_index <= m_fill_index); // or something

    m_candidates[seed_index].push_back(track);
  }

  /* void SortByPhi() */
  /* { */
  /*   std::sort(m_candidates.begin(), m_candidates.begin() + m_fill_index, sortTrksByPhiMT); */
  /* } */
};

class EventOfCombCandidates
{
public:
  std::vector<EtaBinOfCombCandidates> m_etabins_of_comb_candidates;

public:
  EventOfCombCandidates() :
    m_etabins_of_comb_candidates(Config::nEtaBin)
  {}

  void Reset()
  {
    for (auto &i : m_etabins_of_comb_candidates)
    {
      i.Reset();
    }
  }

  void InsertSeed(const Track& seed)
  {
    int bin = getEtaBin(seed.posEta());
    if ( bin != -1 )
    {
      m_etabins_of_comb_candidates[bin].InsertSeed(seed);
    } 
#ifdef DEBUG
    else { dprint("excluding seed with r=" << seed.posR() << " etaBin=" << bin); };
#endif
  }

  void InsertCandidate(const Track& track, int seed_index)
  {
    int bin = getEtaBin(track.posEta());
    m_etabins_of_comb_candidates[bin].InsertTrack(track,seed_index);
  }

};

#endif
