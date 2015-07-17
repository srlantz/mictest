#ifndef CandCloner_h
#define CandCloner_h

#include "MkFitter.h"

#include <vector>


class CandCloner
{
public:
  CandCloner()
  {
    m_fitter = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(0);
  }

  ~CandCloner()
  {
     _mm_free(m_fitter);
  }

  void begin_eta_bin(EtaBinOfCombCandidates * eb_o_ccs, int start_seed, int n_seeds)
  {
    mp_etabin_of_comb_candidates = eb_o_ccs;
    m_start_seed = start_seed;
    m_n_seeds    = n_seeds;
    m_hits_to_add.resize(n_seeds);
    // XXX Should resize vectors in m_hits_to_add to whatever makes sense.
  }

  void begin_layer(BunchOfHits *b_o_hs, int lay)
  {
    mp_bunch_of_hits = b_o_hs;

    m_layer = lay;

    m_fitter->SetNhits(m_layer);//here again assuming one hit per layer

    m_idx_max = m_idx_max_prev = 0;
  }

  void begin_iteration()
  {
    m_idx_max_prev = m_idx_max;
  }

  void add_cand(int idx, const MkFitter::IdxChi2List& cand_info)
  {
    m_hits_to_add[idx].push_back(cand_info);

    m_idx_max = std::max(m_idx_max, idx);
  }

  void end_iteration()
  {
    int proc_n   = m_idx_max - m_idx_max_prev;

    printf("CandCloner::end_iteration process %d, max_prev=%d, max=%d\n", proc_n, m_idx_max_prev, m_idx_max);

    // At this point can forward proc_n seed vectors for further processing.
    // Need to:
    // - know where to put the good candidates
    // - have a thread that does it
    //
    // Hmmh, it seems we would prefer to wait until we have at least NN of them.
    // Or have this explicitly processed in the cloner thread ... will require more
    // communication though.
  }

  void CutAndPastaFromBuildTestMPlex();

  // ----------------------------------------------------------------

  // eventually, protected or private

  int  m_idx_max, m_idx_max_prev;
  std::vector<std::vector<MkFitter::IdxChi2List>> m_hits_to_add;

  BunchOfHits            *mp_bunch_of_hits;
  EtaBinOfCombCandidates *mp_etabin_of_comb_candidates;
  int m_start_seed, m_n_seeds;
  int m_layer;

  MkFitter *m_fitter;
};

#endif
