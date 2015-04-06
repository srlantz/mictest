#include "MkFitter.h"

#include "PropagationMPlex.h"
#include "KalmanUtilsMPlex.h"

void MkFitter::CheckAlignment()
{
  printf("MkFitter alignment check:\n");
  Matriplex::align_check("  Err[0]   =", &Err[0].fArray[0]);
  Matriplex::align_check("  Err[1]   =", &Err[1].fArray[0]);
  Matriplex::align_check("  Par[0]   =", &Par[0].fArray[0]);
  Matriplex::align_check("  Par[1]   =", &Par[1].fArray[0]);
  Matriplex::align_check("  msErr[0] =", &msErr[0].fArray[0]);
  Matriplex::align_check("  msPar[0] =", &msPar[0].fArray[0]);
}

void MkFitter::PrintPt(int idx)
{
  for (int i = 0; i < NN; ++i)
  {
    printf("%5.2f  ", hipo(Par[idx].At(i, 3, 0), Par[idx].At(i, 4, 0)));
  }
}

//==============================================================================

void MkFitter::InputTracksAndHits(std::vector<Track>& tracks, int beg, int end)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Track &trk = tracks[i];

    SeedIdx(itrack, 0, 0) = trk.label();

    Err[iC].CopyIn(itrack, trk.errors().Array());
    Par[iC].CopyIn(itrack, trk.parameters().Array());

    Chg(itrack, 0, 0) = trk.charge();
    Chi2(itrack, 0, 0) = trk.chi2();

    for (int hi = 0; hi < Nhits; ++hi)
    {

      Hit &hit = trk.hitsVector()[hi];

      msErr[hi].CopyIn(itrack, hit.error().Array());
      msPar[hi].CopyIn(itrack, hit.parameters().Array());
    }
  }
}

void MkFitter::InputTracksAndHitIdx(std::vector<Track>& tracks, int beg, int end)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {

    Track &trk = tracks[i];

    SeedIdx(itrack, 0, 0) = trk.label();

    Err[iC].CopyIn(itrack, trk.errors().Array());
    Par[iC].CopyIn(itrack, trk.parameters().Array());

    Chg(itrack, 0, 0) = trk.charge();
    Chi2(itrack, 0, 0) = trk.chi2();

    for (int hi = 0; hi < Nhits; ++hi)
    {

      HitsIdx[hi](itrack, 0, 0) = trk.getHitIdx(hi);//dummy value for now

    }
  }
}

void MkFitter::InputTracksAndHitIdx(std::vector<std::vector<Track> >& tracks, std::vector<std::pair<int,int> >& idxs, int beg, int end)
{
  // Assign track parameters to initial state and copy hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {

    Track &trk = tracks[idxs[i].first][idxs[i].second];

    SeedIdx(itrack, 0, 0) = idxs[i].first;

    Err[iC].CopyIn(itrack, trk.errors().Array());
    Par[iC].CopyIn(itrack, trk.parameters().Array());

    Chg(itrack, 0, 0) = trk.charge();
    Chi2(itrack, 0, 0) = trk.chi2();

    for (int hi = 0; hi < Nhits; ++hi)
    {

      HitsIdx[hi](itrack, 0, 0) = trk.getHitIdx(hi);//dummy value for now

    }
  }
}

int MkFitter::countInvalidHits(int itrack)
{
  int result = 0;
  for (int hi = 0; hi < Nhits; ++hi)
    {
      if (HitsIdx[hi](itrack, 0, 0)<0) result++;
    }
  return result;
}

void MkFitter::InputTracksOnly(std::vector<Track>& tracks, int beg, int end)
{
  // Assign track parameters to initial state, do NOT copy hit values.
  // Used for benchmarking the fitting with less "copy-in" load.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Track &trk = tracks[i];

    Err[iC].CopyIn(itrack, trk.errors().Array());
    Par[iC].CopyIn(itrack, trk.parameters().Array());

    Chg(itrack, 0, 0) = trk.charge();
    Chi2(itrack, 0, 0) = trk.chi2();
    SeedIdx(itrack, 0, 0) = trk.label();
  }
}

void MkFitter::InputHitsOnly(std::vector<Hit>& hits, int beg, int end)
{
  // Push hit values in.

  // This might not be true for the last chunk!
  // assert(end - beg == NN);

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Hit &hit = hits[itrack];

    msErr[Nhits].CopyIn(itrack, hit.error().Array());
    msPar[Nhits].CopyIn(itrack, hit.parameters().Array());
  }
  Nhits++;
}

void MkFitter::FitTracks()
{
  // Fitting loop.

  for (int hi = 0; hi < Nhits; ++hi)
  {
    // Note, charge is not passed (line propagation).
    // propagateLineToRMPlex(Err[iC], Par[iC], msErr[hi], msPar[hi],
    //                       Err[iP], Par[iP]);

    propagateHelixToRMPlex(Err[iC], Par[iC], Chg, msPar[hi],
                           Err[iP], Par[iP]);

    updateParametersMPlex(Err[iP], Par[iP], msErr[hi], msPar[hi],
                          Err[iC], Par[iC]);
  }

  // XXXXX What's with chi2?
}

void MkFitter::OutputTracks(std::vector<Track>& tracks, int beg, int end, int iCP)
{
  // Copies last track parameters (updated) into Track objects.
  // The tracks vector should be resized to allow direct copying.

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Err[iCP].CopyOut(itrack, tracks[i].errors().Array());
    Par[iCP].CopyOut(itrack, tracks[i].parameters().Array());

    tracks[i].charge() = Chg(itrack, 0, 0);

    // XXXXX chi2 is not set (also not in SMatrix fit, it seems)
    tracks[i].setChi2(Chi2(itrack, 0, 0));
    tracks[i].setLabel(SeedIdx(itrack, 0, 0));
  }
}

void MkFitter::OutputFittedTracksAndHits(std::vector<Track>& tracks, int beg, int end)
{
  // Copies last track parameters (updated) into Track objects and up to Nhits.
  // The tracks vector should be resized to allow direct copying.

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Err[iC].CopyOut(itrack, tracks[i].errors().Array());
    Par[iC].CopyOut(itrack, tracks[i].parameters().Array());

    tracks[i].charge() = Chg(itrack, 0, 0);
    tracks[i].setChi2(Chi2(itrack, 0, 0));
    tracks[i].setLabel(SeedIdx(itrack, 0, 0));

    // XXXXX chi2 is not set (also not in SMatrix fit, it seems)

    tracks[i].resetHits();
    for (int hi = 0; hi < Nhits; ++hi)
    {
      Hit hit;

      msErr[hi].CopyOut(itrack, hit.error().Array());
      msPar[hi].CopyOut(itrack, hit.parameters().Array());

      tracks[i].addHit(hit,0.);
      tracks[i].addHitIdx(HitsIdx[hi](itrack, 0, 0),0.);
    }
  }
}


void MkFitter::OutputFittedTracksAndHitIdx(std::vector<Track>& tracks, int beg, int end)
{
  // Copies last track parameters (updated) into Track objects and up to Nhits.
  // The tracks vector should be resized to allow direct copying.

  int itrack = 0;
  for (int i = beg; i < end; ++i, ++itrack)
  {
    Err[iC].CopyOut(itrack, tracks[i].errors().Array());
    Par[iC].CopyOut(itrack, tracks[i].parameters().Array());

    tracks[i].charge() = Chg(itrack, 0, 0);
    tracks[i].setChi2(Chi2(itrack, 0, 0));
    tracks[i].setLabel(SeedIdx(itrack, 0, 0));

    // XXXXX chi2 is not set (also not in SMatrix fit, it seems)

    tracks[i].resetHits();
    for (int hi = 0; hi < Nhits; ++hi)
    {
      tracks[i].addHitIdx(HitsIdx[hi](itrack, 0, 0),0.);
    }

  }
}

void MkFitter::PropagateTracksToR(float R)
{

    propagateHelixToRMPlex(Err[iC], Par[iC], Chg, R,
                           Err[iP], Par[iP]);

}

//fixme: do it properly with phi segmentation
void MkFitter::AddBestHit(std::vector<Hit>& lay_hits, int firstHit, int lastHit, int beg, int end)
{

  //fixme solve ambiguity NN vs beg-end
  float minChi2[NN];
  std::fill_n(minChi2, NN, 9999.);
  int bestHit[NN];
  std::fill_n(bestHit, NN, -1);

  //outer loop over hits, so that tracks can be vectorized
  int ih = 0;
  for (int layhit = firstHit; layhit<lastHit; ++ih, ++layhit) {

    Hit &hit = lay_hits[layhit];
#ifdef DEBUG
    std::cout << "consider hit #" << ih << " of " << lay_hits.size() << std::endl;
    std::cout << "hit x=" << hit.position()[0] << " y=" << hit.position()[1] << std::endl;      
#endif

    //create a dummy matriplex with N times the same hit
    //fixme this is just because I'm lazy and do not want to rewrite the chi2 calculation for simple vectors mixed with matriplex
    MPlexHS msErr_oneHit;
    MPlexHV msPar_oneHit;
    int itrack = 0;
    //fixme: please vectorize me...
    //#pragma simd
    for (int i = beg; i < end; ++i, ++itrack)
    {
      msErr_oneHit.CopyIn(itrack, hit.error().Array());
      msPar_oneHit.CopyIn(itrack, hit.parameters().Array());
    }

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    computeChi2MPlex(Err[iP], Par[iP],msErr_oneHit, msPar_oneHit, outChi2);

    //update best hit in case chi2<minChi2
    itrack = 0;
    //fixme: please vectorize me...
#pragma simd
    for (int i = beg; i < end; ++i, ++itrack)
    {
      float chi2 = fabs(outChi2[itrack]);//fixme negative chi2 sometimes...
#ifdef DEBUG
      std::cout << "chi2=" << chi2 << " minChi2[itrack]=" << minChi2[itrack] << std::endl;      
#endif
      if (chi2<minChi2[itrack]) 
      {
	minChi2[itrack]=chi2;
	bestHit[itrack]=layhit;
      }
    }

  }//end loop over hits

  //copy in MkFitter the hit with lowest chi2
  int itrack = 0;
  //fixme: please vectorize me...
#pragma simd
  for (int i = beg; i < end; ++i, ++itrack)
  {
    //fixme decide what to do in case no hit found
    if (bestHit[itrack] >= 0)
    {
      Hit   &hit  = lay_hits[ bestHit[itrack] ];
      float &chi2 = minChi2[itrack];

#ifdef DEBUG
      std::cout << "ADD BEST HIT FOR TRACK #" << i << std::endl;
      std::cout << "prop x=" << Par[iP].ConstAt(itrack, 0, 0) << " y=" << Par[iP].ConstAt(itrack, 1, 0) << std::endl;      
      std::cout << "copy in hit #" << bestHit[itrack] << " x=" << hit.position()[0] << " y=" << hit.position()[1] << std::endl;    
#endif

      msErr[Nhits].CopyIn(itrack, hit.error().Array());
      msPar[Nhits].CopyIn(itrack, hit.parameters().Array());
      Chi2(itrack, 0, 0) += chi2;
      HitsIdx[Nhits](itrack, 0, 0) = bestHit[itrack];//fixme should add the offset
    }
    else
    {
#ifdef DEBUG
      std::cout << "ADD FAKE HIT FOR TRACK #" << i << std::endl;
#endif

      msErr[Nhits].SetDiagonal3x3(itrack, 666);
      msPar[Nhits](itrack,0,0) = Par[iP](itrack,0,0);
      msPar[Nhits](itrack,1,0) = Par[iP](itrack,1,0);
      msPar[Nhits](itrack,2,0) = Par[iP](itrack,2,0);
      HitsIdx[Nhits](itrack, 0, 0) = -1;

      // Don't update chi2
    }
  }

  //now update the track parameters with this hit (note that some calculations are already done when computing chi2... not sure it's worth caching them?)
#ifdef DEBUG
  std::cout << "update parameters" << std::endl;
#endif
  updateParametersMPlex(Err[iP], Par[iP], msErr[Nhits], msPar[Nhits],
			Err[iC], Par[iC]);

}

void MkFitter::FindCandidates(std::vector<Hit>& lay_hits, int firstHit, int lastHit, int beg, int end, std::vector<std::vector<Track> >& tmp_candidates, int offset)
{

  //outer loop over hits, so that tracks can be vectorized
  int ih = firstHit;
  for ( ; ih<lastHit; ++ih)
  {

    Hit &hit = lay_hits[ih];
#ifdef DEBUG
    std::cout << "consider hit #" << ih << " of " << lay_hits.size() << std::endl;
    std::cout << "hit x=" << hit.position()[0] << " y=" << hit.position()[1] << std::endl;      
#endif

    //create a dummy matriplex with N times the same hit
    //fixme this is just because I'm lazy and do not want to rewrite the chi2 calculation for simple vectors mixed with matriplex
    MPlexHS msErr_oneHit;
    MPlexHV msPar_oneHit;
    int itrack = 0;
    //fixme: please vectorize me...
    for (int i = beg; i < end; ++i, ++itrack)
    {
      msErr_oneHit.CopyIn(itrack, hit.error().Array());
      msPar_oneHit.CopyIn(itrack, hit.parameters().Array());
    }

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    computeChi2MPlex(Err[iP], Par[iP],msErr_oneHit, msPar_oneHit, outChi2);

    //now update the track parameters with this hit (note that some calculations are already done when computing chi2, to be optimized)
    //this is not needed for candidates the hit is not added to, but it's vectorized so doing it serially below should take the same time
    //still it's a waste of time in case the hit is not added to any of the candidates, so check beforehand that at least one cand needs update
    bool oneCandPassCut = false;
    itrack = 0;
    for (int i = beg; i < end; ++i, ++itrack)
      {
	float chi2 = fabs(outChi2[itrack]);//fixme negative chi2 sometimes...
#ifdef DEBUG
	std::cout << "chi2=" << chi2 << std::endl;
#endif
	if (chi2<Config::chi2Cut)
	  {
	    oneCandPassCut = true;
	    break;
	  }
      }

    if (oneCandPassCut) { 

      updateParametersMPlex(Err[iP], Par[iP], msErr_oneHit, msPar_oneHit,
			    Err[iC], Par[iC]);
#ifdef DEBUG
      std::cout << "update parameters" << std::endl;
      std::cout << "propagated track parameters x=" << Par[iP].ConstAt(itrack, 0, 0) << " y=" << Par[iP].ConstAt(itrack, 1, 0) << std::endl;
      std::cout << "               hit position x=" << msPar[iP].ConstAt(itrack, 0, 0) << " y=" << msPar[iP].ConstAt(itrack, 1, 0) << std::endl;
      std::cout << "   updated track parameters x=" << Par[iC].ConstAt(itrack, 0, 0) << " y=" << Par[iC].ConstAt(itrack, 1, 0) << std::endl;
#endif

      //create candidate with hit in case chi2<Config::chi2Cut
      itrack = 0;
      //fixme: please vectorize me... (not sure it's possible in this case)
      for (int i = beg; i < end; ++i, ++itrack)
	{
	  float chi2 = fabs(outChi2[itrack]);//fixme negative chi2 sometimes...
#ifdef DEBUG
	  std::cout << "chi2=" << chi2 << std::endl;      
#endif
	  if (chi2<Config::chi2Cut)
	    {
#ifdef DEBUG
	      std::cout << "chi2 cut passed, creating new candidate" << std::endl;
#endif
	      //create a new candidate and fill the reccands_tmp vector
	      Track newcand;
	      newcand.resetHits();//probably not needed
	      newcand.charge() = Chg(itrack, 0, 0);
	      newcand.setChi2(Chi2(itrack, 0, 0));
	      for (int hi = 0; hi < Nhits; ++hi)
		{
		  newcand.addHitIdx(HitsIdx[hi](itrack, 0, 0),0.);//this should be ok since we already set the chi2 above
		}
#ifdef DEBUG
	      std::cout << "output new hit with x=" << hit.position()[0] << std::endl;
#endif

	      newcand.addHitIdx(ih,chi2);
	      //set the track state to the updated parameters
	      Err[iC].CopyOut(itrack, newcand.errors().Array());
	      Par[iC].CopyOut(itrack, newcand.parameters().Array());
	      
#ifdef DEBUG
	      std::cout << "updated track parameters x=" << newcand.parameters()[0] << " y=" << newcand.parameters()[1] << std::endl;
#endif
	      
	      tmp_candidates[SeedIdx(itrack, 0, 0)-offset].push_back(newcand);
	    }
	}
    }//end if (oneCandPassCut)

  }//end loop over hits

  //now add invalid hit
  int itrack = 0;
  //fixme: please vectorize me...
  for (int i = beg; i < end; ++i, ++itrack)
    {
      if (countInvalidHits(itrack)>0) continue;//check this is ok for vectorization //fixme not optimal
      Track newcand;
      newcand.resetHits();//probably not needed
      newcand.charge() = Chg(itrack, 0, 0);
      newcand.setChi2(Chi2(itrack, 0, 0));
      for (int hi = 0; hi < Nhits; ++hi)
	{
	  newcand.addHitIdx(HitsIdx[hi](itrack, 0, 0),0.);//this should be ok since we already set the chi2 above
	}
      newcand.addHitIdx(-1,0.);
      //set the track state to the propagated parameters
      Err[iP].CopyOut(itrack, newcand.errors().Array());
      Par[iP].CopyOut(itrack, newcand.parameters().Array());	      
      tmp_candidates[SeedIdx(itrack, 0, 0)-offset].push_back(newcand);
    }

}


void MkFitter::GetHitRange(std::vector<std::vector<BinInfo> >& segmentMapLay_, int beg, int end,
                           const float etaDet, int& firstHit, int& lastHit)
{

    int itrack = 0;


    const float eta_predx = Par[iP].ConstAt(itrack, 0, 0);
    const float eta_predy = Par[iP].ConstAt(itrack, 1, 0);
    const float eta_predz = Par[iP].ConstAt(itrack, 2, 0);
    float eta = getEta(eta_predx,eta_predy,eta_predz);
    //protect against anomalous eta (should go into getEtaPartition maybe?)
    if (fabs(eta) > etaDet) eta = (eta>0 ? etaDet*0.99 : -etaDet*0.99);
    unsigned int etabin = getEtaPartition(eta,etaDet);

    //cannot be vectorized, I think
    for (int i = beg; i < end; ++i, ++itrack)
      {

        const float predx = Par[iP].ConstAt(itrack, 0, 0);
        const float predy = Par[iP].ConstAt(itrack, 1, 0);
        const float predz = Par[iP].ConstAt(itrack, 2, 0);

	float phi = getPhi(predx,predy);

	const float px2py2 = predx*predx+predy*predy; // predicted radius^2

	const float dphidx = -predy/px2py2;
	const float dphidy =  predx/px2py2;
	const float dphi2  = dphidx*dphidx*(Err[iP].ConstAt(itrack, 0, 0) /*propState.errors.At(0,0)*/) +
	  dphidy*dphidy*(Err[iP].ConstAt(itrack, 1, 1) /*propState.errors.At(1,1)*/) +
	  2*dphidx*dphidy*(Err[iP].ConstAt(itrack, 0, 1) /*propState.errors.At(0,1)*/);
  
	const float dphi   =  sqrt(std::fabs(dphi2));//how come I get negative squared errors sometimes?
	const float nSigmaDphi = std::min(std::max(Config::nSigma*dphi,(float) Config::minDPhi), float(M_PI/1.));//fixme
	//const float nSigmaDphi = Config::nSigma*dphi;

	//if (nSigmaDphi>0.3) std::cout << "window MX: " << predx << " " << predy << " " << predz << " " << Err[iP].ConstAt(itrack, 0, 0) << " " << Err[iP].ConstAt(itrack, 1, 1) << " " << Err[iP].ConstAt(itrack, 0, 1) << " " << nSigmaDphi << std::endl;

	const float dphiMinus = normalizedPhi(phi-nSigmaDphi);
	const float dphiPlus  = normalizedPhi(phi+nSigmaDphi);
  
	const auto phiBinMinus = getPhiPartition(dphiMinus);
	const auto phiBinPlus  = getPhiPartition(dphiPlus);
    
	BinInfo binInfoMinus = segmentMapLay_[etabin][int(phiBinMinus)];
	BinInfo binInfoPlus  = segmentMapLay_[etabin][int(phiBinPlus)];

	//fixme: temporary to avoid wrapping
	if (binInfoMinus > binInfoPlus)
        {
	  unsigned int phibin = getPhiPartition(phi);
	  binInfoMinus = segmentMapLay_[etabin][phibin];
	  binInfoPlus  = segmentMapLay_[etabin][phibin];
	}

	//fixme: if more than one eta bin we are looping over a huge range (need to make sure we are in the same eta bin)

#ifdef DEBUG
	std::cout << "propagated track parameters eta=" << eta << " bin=" << etabin << " begin=" << binInfoMinus.first << " end=" << binInfoPlus.first+binInfoPlus.second << std::endl;
#endif
	if (firstHit==-1 || binInfoMinus.first<firstHit) firstHit =  binInfoMinus.first;
	if (lastHit==-1 || (binInfoPlus.first+binInfoPlus.second)>lastHit) lastHit = binInfoPlus.first+binInfoPlus.second;
      }
#ifdef DEBUG
    std::cout << "found range firstHit=" << firstHit << " lastHit=" << lastHit << std::endl;
#endif
}


// ======================================================================================
// MT methods
// ======================================================================================

void MkFitter::SelectHitRanges(BunchOfHits &bunch_of_hits)
{
  // must store hit vector into a data member so it can be used in hit selection.
  // or ... can have it passed into other functions.
  // somewhat yucky, either way.

  // Also, must store two ints per Matriplex elements ... first index and size.
  // These are XPos and XSize

  // vecorized for
#pragma simd
  for (int itrack = 0; itrack < NN; ++itrack)
  {
    // Hmmh ... this should all be solved by partitioning ... let's try below ...
    //
    // float eta = getEta(eta_predx,eta_predy,eta_predz);
    // //protect against anomalous eta (should go into getEtaPartition maybe?)
    // if (fabs(eta) > etaDet) eta = (eta>0 ? etaDet*0.99 : -etaDet*0.99);
    // unsigned int etabin = getEtaPartition(eta,etaDet);

    const float predx = Par[iP].ConstAt(itrack, 0, 0);
    const float predy = Par[iP].ConstAt(itrack, 1, 0);
    const float predz = Par[iP].ConstAt(itrack, 2, 0);

    float phi = getPhi(predx,predy);

    const float px2py2 = predx*predx+predy*predy; // predicted radius^2

    const float dphidx = -predy/px2py2;
    const float dphidy =  predx/px2py2;
    const float dphi2  =     dphidx*dphidx*(Err[iP].ConstAt(itrack, 0, 0) /*propState.errors.At(0,0)*/) +
                             dphidy*dphidy*(Err[iP].ConstAt(itrack, 1, 1) /*propState.errors.At(1,1)*/) +
                         2 * dphidx*dphidy*(Err[iP].ConstAt(itrack, 0, 1) /*propState.errors.At(0,1)*/);

    const float dphi       = sqrtf(std::fabs(dphi2));//how come I get negative squared errors sometimes? MT -- how small?
    const float nSigmaDphi = std::min(std::max(Config::nSigma*dphi,(float) Config::minDPhi), float(M_PI/1.));//fixme
    //const float nSigmaDphi = Config::nSigma*dphi;

    //if (nSigmaDphi>0.3) 
    //std::cout << "window MX: " << predx << " " << predy << " " << predz << " " << Err[iP].ConstAt(itrack, 0, 0) << " " << Err[iP].ConstAt(itrack, 1, 1) << " " << Err[iP].ConstAt(itrack, 0, 1) << " " << nSigmaDphi << std::endl;

    const float dphiMinus = normalizedPhi(phi-nSigmaDphi);
    const float dphiPlus  = normalizedPhi(phi+nSigmaDphi);

#ifdef DEBUG
    std::cout << "dphi = " << dphi  << ", dphi2 = " << dphi2 << ", nSigmaDphi = " << nSigmaDphi << ", nSigma = " << Config::nSigma << std::endl;
    std::cout << "phiMinus = " << dphiMinus << ", phiPlus = " << dphiPlus << std::endl;
#endif

    int   phiBinMinus = getPhiPartition(dphiMinus);
    int   phiBinPlus  = getPhiPartition(dphiPlus);

#ifdef DEBUG
    std::cout << "phiBinMinus = " << phiBinMinus << ", phiBinPlus = " << phiBinPlus << std::endl;
#endif

    phiBinMinus = std::max(0,phiBinMinus);
    phiBinMinus = std::min(Config::nPhiPart-1,phiBinMinus);
    phiBinPlus = std::max(0,phiBinPlus);
    phiBinPlus = std::min(Config::nPhiPart-1,phiBinPlus);


    BinInfo binInfoMinus = bunch_of_hits.m_phi_bin_infos[int(phiBinMinus)];
    BinInfo binInfoPlus  = bunch_of_hits.m_phi_bin_infos[int(phiBinPlus)];

    //fixme: temporary to avoid wrapping
    if (binInfoMinus > binInfoPlus)
    {
      int phibin = getPhiPartition(phi);
      phibin = std::max(0,phibin);
      phibin = std::min(Config::nPhiPart-1,phibin);
      binInfoMinus = bunch_of_hits.m_phi_bin_infos[phibin];
      binInfoPlus  = bunch_of_hits.m_phi_bin_infos[phibin];
    }

    // fixme: if more than one eta bin we are looping over a huge range (need to make sure we are in the same eta bin)
    // MT -- this is being fixed here. Plus improving other things ... and screwing up others :)

#ifdef DEBUG
    std::cout << "bin info begin=" << binInfoMinus.first << " end=" << binInfoPlus.first+binInfoPlus.second << std::endl;
#endif

    XHitBegin.At(itrack, 0, 0) = binInfoMinus.first;
    XHitEnd  .At(itrack, 0, 0) = binInfoPlus .first + binInfoPlus.second;

#ifdef DEBUG
    std::cout << "found range firstHit=" << XHitBegin.At(itrack, 0, 0) << " lastHit=" << XHitEnd  .At(itrack, 0, 0) << std::endl;
#endif

  }
}


//==============================================================================
// AddBestHit()
//==============================================================================

void MkFitter::AddBestHit(BunchOfHits &bunch_of_hits)
{
  //fixme solve ambiguity NN vs beg-end
  float minChi2[NN];
  std::fill_n(minChi2, NN, 100.); // XXXX MT was 9999
  int bestHit[NN];
  std::fill_n(bestHit, NN, -1);

  const char *varr      = (char*) bunch_of_hits.m_hits;

  const int   off_error = (char*) bunch_of_hits.m_hits[0].error().Array()      - varr;
  const int   off_param = (char*) bunch_of_hits.m_hits[0].parameters().Array() - varr;

  // int idx[NN]      __attribute__((aligned(64)));
  int *idx  = m_pf_idx;
  int idx_chew[NN] __attribute__((aligned(64)));

  int maxSize = -1;

  // Determine maximum number of hits for tracks in the collection.
  // At the same time prefetch the first set of hits to L1 and the second one to L2.
  for (int it = 0; it < NN; ++it)
  {
    int off = XHitBegin.At(it, 0, 0) * sizeof(Hit);

    _mm_prefetch(varr + off, _MM_HINT_T0);
    _mm_prefetch(varr + sizeof(Hit) + off, _MM_HINT_T1);

    idx[it]      = off;
    idx_chew[it] = it*sizeof(Hit);

    if (XHitEnd.At(it, 0, 0) - XHitBegin.At(it, 0, 0) > maxSize)
    {
      maxSize = (XHitEnd.At(it, 0, 0) - XHitBegin.At(it, 0, 0));
    }
  }

  // XXXX MT Uber hack to avoid tracks with like 300 hits to process.
  if (maxSize > 25) maxSize = 25;

#if defined(MIC_INTRINSICS)
  //__m512i vi = _mm512_setr_epi32(idx[ 0], idx[ 1], idx[ 2], idx[ 3], idx[ 4], idx[ 5], idx[ 6], idx[ 7],
  //                               idx[ 8], idx[ 9], idx[10], idx[11], idx[12], idx[13], idx[14], idx[15]);
  __m512i vi      = _mm512_load_epi32(idx);
  __m512i vi_chew = _mm512_load_epi32(idx_chew);
#endif

// Has basically no effect, it seems.
//#pragma noprefetch
  for (int hit_cnt = 0; hit_cnt < maxSize; ++hit_cnt, varr += sizeof(Hit))
  {
    //fixme what if size is zero???

    // {
    //   std::unique_lock<std::mutex> lk(m_moo);
    //   m_pf_arr  = varr;
    //   m_pf_done = false;
    //   m_cnd.notify_one();
    // }
    // Prefetch to L2 the hits we'll process after two loops iterations.
    // Ideally this would be initiated before coming here, for whole bunch_of_hits.m_hits vector.
    // for (int itrack = 0; itrack < NN; ++itrack)
    // {
    //   _mm_prefetch(varr + 2*sizeof(Hit) + idx[itrack], _MM_HINT_T1);
    // }

#if defined(MIC_INTRINSICS)
    msErr[Nhits].SlurpIn(varr + off_error, vi);
    msPar[Nhits].SlurpIn(varr + off_param, vi);

    //char arr_chew[NN*sizeof(Hit)] __attribute__((aligned(64)));

    // This is a hack ... we know sizeof(Hit) = 64 = cache line = vector width.
    // ChewIn runs about 2% slower than SlurpIn().
    //msErr[Nhits].ChewIn(varr, off_error, idx, arr_chew, vi_chew);
    //msPar[Nhits].ChewIn(varr, off_param, idx, arr_chew, vi_chew);

    // This is a hack ... we know sizeof(Hit) = 64 = cache line = vector width.
    // Contaginate + Plexify runs just about as fast as SlurpIn().
    //msErr[Nhits].Contaginate(varr, idx, arr_chew);
    //msErr[Nhits].Plexify(arr_chew + off_error, vi_chew);
    //msPar[Nhits].Plexify(arr_chew + off_param, vi_chew);
#else
    msErr[Nhits].SlurpIn(varr + off_error, idx);
    msPar[Nhits].SlurpIn(varr + off_param, idx);
#endif

    //now compute the chi2 of track state vs hit
    MPlexQF outChi2;
    computeChi2MPlex(Err[iP], Par[iP], msErr[Nhits], msPar[Nhits], outChi2);

    // Prefetch to L1 the hits we'll process next ...
    for (int itrack = 0; itrack < NN; ++itrack)
    {
      _mm_prefetch(varr + sizeof(Hit) + idx[itrack], _MM_HINT_T0);
    }

    //update best hit in case chi2<minChi2
#pragma simd
    for (int itrack = 0; itrack < NN; ++itrack)
    {
      float chi2 = fabs(outChi2[itrack]);//fixme negative chi2 sometimes...
#ifdef DEBUG
      std::cout << "chi2=" << chi2 << " minChi2[itrack]=" << minChi2[itrack] << std::endl;      
#endif
      if (chi2 < minChi2[itrack]) 
      {
        minChi2[itrack]=chi2;
        bestHit[itrack]=hit_cnt;
      }
    }

  } // end loop over hits

  //copy in MkFitter the hit with lowest chi2
  for (int itrack = 0; itrack < NN; ++itrack)
  {
    _mm_prefetch((const char*) & bunch_of_hits.m_hits[XHitBegin.At(itrack, 0, 0) + bestHit[itrack]], _MM_HINT_T0);
  }
    
#pragma simd
  for (int itrack = 0; itrack < NN; ++itrack)
  {
    //fixme decide what to do in case no hit found
    if (bestHit[itrack] >= 0)
    {
      Hit   &hit  = bunch_of_hits.m_hits[ XHitBegin.At(itrack, 0, 0) + bestHit[itrack] ];
      float &chi2 = minChi2[itrack];
	  
#ifdef DEBUG
      std::cout << "ADD BEST HIT FOR TRACK #" << itrack << std::endl;
      std::cout << "prop x=" << Par[iP].ConstAt(itrack, 0, 0) << " y=" << Par[iP].ConstAt(itrack, 1, 0) << std::endl;      
      std::cout << "copy in hit #" << bestHit[itrack] << " x=" << hit.position()[0] << " y=" << hit.position()[1] << std::endl;    
#endif
	  
      msErr[Nhits].CopyIn(itrack, hit.error().Array());
      msPar[Nhits].CopyIn(itrack, hit.parameters().Array());
      Chi2(itrack, 0, 0) += chi2;
      HitsIdx[Nhits](itrack, 0, 0) = XHitBegin.At(itrack, 0, 0) + bestHit[itrack];
    }
    else
    {
#ifdef DEBUG
      std::cout << "ADD FAKE HIT FOR TRACK #" << itrack << std::endl;
#endif
	  
      msErr[Nhits].SetDiagonal3x3(itrack, 666);
      msPar[Nhits](itrack,0,0) = Par[iP](itrack,0,0);
      msPar[Nhits](itrack,1,0) = Par[iP](itrack,1,0);
      msPar[Nhits](itrack,2,0) = Par[iP](itrack,2,0);
      HitsIdx[Nhits](itrack, 0, 0) = -1;
	  
      // Don't update chi2
    }
  }
  
  //now update the track parameters with this hit (note that some calculations are already done when computing chi2... not sure it's worth caching them?)
#ifdef DEBUG
  std::cout << "update parameters" << std::endl;
#endif
  updateParametersMPlex(Err[iP], Par[iP], msErr[Nhits], msPar[Nhits],
			Err[iC], Par[iC]);

  //std::cout << "Par[iP](0,0,0)=" << Par[iP](0,0,0) << " Par[iC](0,0,0)=" << Par[iC](0,0,0)<< std::endl;
}


//==============================================================================
// Wierd stuff
//==============================================================================

#include <sched.h>

namespace
{
  void pin_this_thread(int loc_id)
  {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(loc_id, &cpuset);
    sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);
  }
}

void MkFitter::PinThisThreadAndSpawnPrefetcher(int cpuid, int cpuid_pref)
{
  pin_this_thread(cpuid);

  {
    std::unique_lock<std::mutex> lk(m_moo);

    m_prefetcher = std::thread([=] { this->Prefetcher(cpuid_pref); });

    m_cnd.wait(lk);
  }
}

void MkFitter::JoinPrefetcher()
{
  // printf("MkFitter::JoinPrefetcher in cpuid %d\n", sched_getcpu());
  {
    std::unique_lock<std::mutex> lk(m_moo);
    m_pf_exit = true;
    m_cnd.notify_one();
  }
  m_prefetcher.join();
}

void MkFitter::Prefetcher(int cpuid)
{
  pin_this_thread(cpuid);

  std::unique_lock<std::mutex> lk(m_moo);

  m_cnd.notify_one();

  while (true)
  {
    m_cnd.wait(lk);

    if (m_pf_exit) break;

    // do work
    // Prefetch to L2 the hits we'll process after two loops iterations.
    // for (int itrack = 0; itrack < NN; ++itrack)
    // {
    //   _mm_prefetch(m_pf_arr + 2*sizeof(Hit) + m_pf_idx[itrack], _MM_HINT_T1);
    //}
    // Prefetch to L1 the hits we'll process in the next loop iteration.
    // for (int itrack = 0; itrack < NN; ++itrack)
    // {
    //   _mm_prefetch(m_pf_arr + sizeof(Hit) + m_pf_idx[itrack], _MM_HINT_T0);
    // }

    for (int i = 0; i < m_bunch_of_hits->m_fill_index; ++i)
    {
      _mm_prefetch((char*) & m_bunch_of_hits->m_hits[i], _MM_HINT_T1);
    }

    
    m_pf_done = true;
  }

  // printf("MkFitter::Prefetcher being joined in cpuid %d\n", sched_getcpu());
}
