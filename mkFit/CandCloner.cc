#include "CandCloner.h"

namespace
{
bool sortCandListByHitsChi2(MkFitter::IdxChi2List cand1,MkFitter::IdxChi2List cand2)
{
  if (cand1.nhits==cand2.nhits) return cand1.chi2<cand2.chi2;
  return cand1.nhits>cand2.nhits;
}
}

void CandCloner::CutAndPastaFromBuildTestMPlex()
{
  //----------------------- BEGIN CLONE ENGINE -----------------------//
  //now we should figure out which maxCand candidates per seed to create
  //this is vectorized, along the lines of the loop above - but now we know the hit to add
  //will have to try to ship to a separate thread	   
  //1) sort the candidates per each seed
  for (int is = 0; is < m_n_seeds; ++is)
  {
    std::vector<MkFitter::IdxChi2List>& hitsToAddForThisSeed = m_hits_to_add[is];
#ifdef DEBUG
    std::cout << "dump seed n " << is << " with input candidates=" << hitsToAddForThisSeed.size() << std::endl;
    for (int ih=0;ih<hitsToAddForThisSeed.size();ih++)
    {
      std::cout << "trkIdx=" << hitsToAddForThisSeed[ih].trkIdx << " hitIdx=" << hitsToAddForThisSeed[ih].hitIdx << " chi2=" <<  hitsToAddForThisSeed[ih].chi2 << std::endl;
      std::cout << "original pt=" << etabin_of_comb_candidates.m_candidates[th_start_seed+is][hitsToAddForThisSeed[ih].trkIdx].pT() << " " 
                << "nHits=" << etabin_of_comb_candidates.m_candidates[th_start_seed+is][hitsToAddForThisSeed[ih].trkIdx].nHits() << " " 
                << "nHitIdx=" << etabin_of_comb_candidates.m_candidates[th_start_seed+is][hitsToAddForThisSeed[ih].trkIdx].nHitIdx() << " " 
                << "chi2=" << etabin_of_comb_candidates.m_candidates[th_start_seed+is][hitsToAddForThisSeed[ih].trkIdx].chi2() << " " 
                << std::endl;
    }
#endif	     
    //sort the damn thing
    std::sort(hitsToAddForThisSeed.begin(), hitsToAddForThisSeed.end(), sortCandListByHitsChi2);
  }
  //2) now create the candidates for the best maxCand, we'll do it vectorized with MkFitter
  //2a) first we need to unroll the vectors
  std::vector<std::pair<int,MkFitter::IdxChi2List> > seed_newcand_idx;
  for (int is = 0; is < m_n_seeds; ++is)
  {
    std::vector<MkFitter::IdxChi2List>& hitsToAddForThisSeed = m_hits_to_add[is];
    for (int ih=0;ih<hitsToAddForThisSeed.size() && ih<Config::maxCand;ih++)
    {
      seed_newcand_idx.push_back(std::pair<int,MkFitter::IdxChi2List>(is,hitsToAddForThisSeed[ih]));
    }
  }
  int theEndNewCand = seed_newcand_idx.size();     
  std::vector<std::vector<Track> > cands_for_next_lay(m_n_seeds);
  for (int iseed=0; iseed < cands_for_next_lay.size(); ++iseed)
  {
    cands_for_next_lay[iseed].reserve(Config::maxCand);
  }
  //2b) vectorized loop
  for (int itrack = 0; itrack < theEndNewCand; itrack += NN)
  {

    int end = std::min(itrack + NN, theEndNewCand);

    ////m_fitter->SetNhits(ilay);//here again assuming one hit per layer
	       
    m_fitter->InputTracksAndHitIdx(mp_etabin_of_comb_candidates->m_candidates, seed_newcand_idx, itrack, end);
	       
    //propagate to layer
#ifdef DEBUG
    std::cout << "propagate to lay=" << ilay+1 << " start from x=" << mkfp->getPar(0, 0, 0) << " y=" << mkfp->getPar(0, 0, 1) << " z=" << mkfp->getPar(0, 0, 2)<< " r=" << getHypot(mkfp->getPar(0, 0, 0), mkfp->getPar(0, 0, 1))
              << " px=" << mkfp->getPar(0, 0, 3) << " py=" << mkfp->getPar(0, 0, 4) << " pz=" << mkfp->getPar(0, 0, 5) << " pT=" << getHypot(mkfp->getPar(0, 0, 3), mkfp->getPar(0, 0, 4)) << std::endl;
#endif
    m_fitter->PropagateTracksToR(4.*(m_layer + 1));//fixme: doesn't need itrack, end?

    //now we need to update with the hit in bunch_of_hits.m_hits[ hitsToAddForThisSeed[ih].hitIdx ]
    //also fill the temp vector of candidates
    m_fitter->UpdateWithHit(*mp_bunch_of_hits, seed_newcand_idx, cands_for_next_lay, m_start_seed, itrack, end);

  }
  //3) now swap the temp vector with the input candidates so that they are used as input for the next layer
  for (int is = 0; is < cands_for_next_lay.size(); ++is)
  {
    if (cands_for_next_lay[is].size() > 0)
    {
      mp_etabin_of_comb_candidates->m_candidates[m_start_seed+is].swap(cands_for_next_lay[is]);
      cands_for_next_lay[is].clear();
    }
    else 
    {
      //we do nothing in the SM version here, I think we should put these in the output and avoid keeping looping over them
    }
  }
  //----------------------- END CLONE ENGINE -----------------------//

}
