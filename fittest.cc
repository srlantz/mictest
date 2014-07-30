#include "fittest.h"
#include "Track.h"
#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"

#ifndef NO_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

#include <iostream>

double runFittingTest(bool saveTree, unsigned int Ntracks)
{
  float pt_mc=0.,pt_fit=0.,pt_err=0.; 
#ifndef NO_ROOT
  TFile* f=0;
  TTree *tree=0;
  if (saveTree) {
    f=TFile::Open("validationtree.root", "recreate");
    tree = new TTree("tree","tree");
    tree->Branch("pt_mc",&pt_mc,"pt_mc");
    tree->Branch("pt_fit",&pt_fit,"pt_fit");
    tree->Branch("pt_err",&pt_err,"pt_err");
  }
#endif

  //these matrices are dummy and can be optimized without multriplying by zero all the world...
  SMatrix36 projMatrix36;
  projMatrix36(0,0)=1.;
  projMatrix36(1,1)=1.;
  projMatrix36(2,2)=1.;
  SMatrix63 projMatrix36T = ROOT::Math::Transpose(projMatrix36);

  std::vector<Track> simtracks;
  for (unsigned int itrack=0;itrack<Ntracks;++itrack) {
    //create the track
    SVector3 pos;
    SVector3 mom;
    SMatrixSym66 covtrk;
    std::vector<Hit> hits;
    int q=0;//set it in setup function
    float pt = 0.5 + g_unif(g_gen) * 9.5;//this input, 0.5<pt<10 GeV  (below ~0.5 GeV does not make 10 layers)
    setupTrackByToyMC(pos,mom,covtrk,hits,q,pt);
    Track simtrk(q,pos,mom,covtrk,hits,0.);
    simtracks.push_back(simtrk);
  }

  double time = dtime();

  for (unsigned int itrack=0;itrack<simtracks.size();++itrack) {

    Track& trk = simtracks[itrack];

    // std::cout << std::endl;
    // std::cout << "processing track #" << itrack << std::endl;
    
    // std::cout << "init x: " << trk.parameters()[0] << " " << trk.parameters()[1] << " " << trk.parameters()[2] << std::endl;
    // std::cout << "init p: " << trk.parameters()[3] << " " << trk.parameters()[4] << " " << trk.parameters()[5] << std::endl;
    // std::cout << "init e: " << std::endl;
    // dumpMatrix(trk.errors());

    std::vector<Hit>& hits = trk.hitsVector();

    TrackState initState = trk.state();
    //make a copy since initState is used at the end to fill the tree
    TrackState updatedState = initState;
    
    bool dump = false;
    
    for (std::vector<Hit>::iterator hit=hits.begin();hit!=hits.end();++hit) {
      
      //for each hit, propagate to hit radius and update track state with hit measurement
      TrackState       propState = propagateHelixToR(updatedState,hit->r());
      MeasurementState measState = hit->measurementState();
      updatedState = updateParameters(propState, measState,projMatrix36,projMatrix36T);
      //updateParameters66(propState, measState, updatedState);//updated state is now modified
      
      if (dump) {
	std::cout << std::endl;
	std::cout << "processing hit #" << hit-hits.begin() << std::endl;
	
	std::cout << "propState.parameters (helix propagation)" << std::endl;
	std::cout << "x: " << propState.parameters[0] << " " << propState.parameters[1] << " " << propState.parameters[2] << std::endl;
	std::cout << "p: " << propState.parameters[3] << " " << propState.parameters[4] << " " << propState.parameters[5] << std::endl;
	std::cout << "propState.errors" << std::endl;
	dumpMatrix(propState.errors);
	
	//TrackState propStateHelix_test = propagateHelixToR_test(updatedState,hit->r());
	//std::cout << "propStateHelix_test.parameters (helix propagation)" << std::endl;
	//std::cout << "x: " << propStateHelix_test.parameters[0] << " " << propStateHelix_test.parameters[1] << " " << propStateHelix_test.parameters[2] << std::endl;
	//std::cout << "p: " << propStateHelix_test.parameters[3] << " " << propStateHelix_test.parameters[4] << " " << propStateHelix_test.parameters[5] << std::endl;
	//std::cout << "propStateHelix_test.errors" << std::endl;
	//dumpMatrix(propStateHelix_test.errors);
	
	// TrackState propStateLine = propagateLineToR(updatedState,hit->r());
	// std::cout << "propStateLine.parameters (line propagation)" << std::endl;
	// std::cout << "x: " << propStateLine.parameters[0] << " " << propStateLine.parameters[1] << " " << propStateLine.parameters[2] << std::endl;
	// std::cout << "p: " << propStateLine.parameters[3] << " " << propStateLine.parameters[4] << " " << propStateLine.parameters[5] << std::endl;
	// std::cout << "propStateLine.errors" << std::endl;
	// dumpMatrix(propStateLine.errors);
	// TrackState propState = propStateLine;
	
	std::cout << "measState.parameters" << std::endl;
	std::cout << "x: " << measState.parameters[0] << " " << measState.parameters[1] << " " << measState.parameters[2] << std::endl;
	std::cout << "measState.errors" << std::endl;
	dumpMatrix(measState.errors);
	
	std::cout << "updatedState" << std::endl;
	std::cout << "x: " << updatedState.parameters[0] << " " << updatedState.parameters[1] << " " << updatedState.parameters[2] << std::endl;
	std::cout << "p: " << updatedState.parameters[3] << " " << updatedState.parameters[4] << " " << updatedState.parameters[5] << std::endl;
	std::cout << "updatedState.errors" << std::endl;
	dumpMatrix(updatedState.errors);	
      }
      
    }
    
    // std::cout << "updatedState" << std::endl;
    // std::cout << "x: " << updatedState.parameters[0] << " " << updatedState.parameters[1] << " " << updatedState.parameters[2] << std::endl;
    // std::cout << "p: " << updatedState.parameters[3] << " " << updatedState.parameters[4] << " " << updatedState.parameters[5] << std::endl;
    // std::cout << "updatedState.errors" << std::endl;
    // dumpMatrix(updatedState.errors);

#ifndef NO_ROOT
    if (saveTree) {
      pt_mc = sqrt(initState.parameters[3]*initState.parameters[3]+initState.parameters[4]*initState.parameters[4]);
      pt_fit = sqrt(updatedState.parameters[3]*updatedState.parameters[3]+updatedState.parameters[4]*updatedState.parameters[4]);
      pt_err = sqrt( updatedState.errors[3][3]*updatedState.parameters[3]*updatedState.parameters[3] +
		     updatedState.errors[4][4]*updatedState.parameters[4]*updatedState.parameters[4] + 
		     2*updatedState.errors[3][4]*updatedState.parameters[3]*updatedState.parameters[4] )/pt_fit;
      tree->Fill();
    }
#endif
  }

  return dtime() - time;

#ifndef NO_ROOT
  if (saveTree) {
    f->Write();
    f->Close();
  }
#endif
}


//==============================================================================
// runFittingTestPlex
//==============================================================================

#ifndef __APPLE__

double runFittingTestPlex(bool saveTree, int Ntracks)
{
  float pt_mc=0.,pt_fit=0.,pt_err=0.; 
#ifndef NO_ROOT
  TFile* f=0;
  TTree *tree=0;
  if (saveTree) {
    f=TFile::Open("validationtree_plex.root", "recreate");
    tree = new TTree("tree","tree");
    tree->Branch("pt_mc",&pt_mc,"pt_mc");
    tree->Branch("pt_fit",&pt_fit,"pt_fit");
    tree->Branch("pt_err",&pt_err,"pt_err");
  }
#endif

  //these matrices are dummy and can be optimized without multriplying by zero all the world...
  SMatrix36 projMatrix36;
  projMatrix36(0,0)=1.;
  projMatrix36(1,1)=1.;
  projMatrix36(2,2)=1.;
  SMatrix63 projMatrix36T = ROOT::Math::Transpose(projMatrix36);

  std::vector<Track> simtracks;
  for (int itrack=0;itrack<Ntracks;++itrack)
  {
    //create the track
    SVector3 pos;
    SVector3 mom;
    SMatrixSym66 covtrk;
    std::vector<Hit> hits;
    int q=0;//set it in setup function
    float pt = 0.5+g_unif(g_gen)*9.5;//this input, 0.5<pt<10 GeV  (below ~0.5 GeV does not make 10 layers)
    setupTrackByToyMC(pos,mom,covtrk,hits,q,pt);
    Track simtrk(q,pos,mom,covtrk,hits,0.);
    simtracks.push_back(simtrk);
  }

  const int Nhits = 10;

  std::vector<TrackState> initStateV(Ntracks);
  //make a copy since initState is used at the end to fill the tree


  MPlexSS psErr(Ntracks);  MPlexMV psPar(Ntracks);
  MPlexSS outErr(Ntracks); MPlexMV outPar(Ntracks);

  std::vector<MPlexSS> Err(2, Ntracks);
  std::vector<MPlexMV> Par(2, Ntracks);

  MPlexQI Chg(Ntracks);

  int iC = 0; // current
  int iP = 1; // propagated

  for (int itrack=0;itrack<Ntracks;++itrack)
  {
    Track &trk = simtracks[itrack];
    Err[iC].Assign(itrack, trk.errors().Array());
    Par[iC].Assign(itrack, trk.parameters().Array());

    initStateV[itrack] = trk.state();
    Chg[0, 0, itrack]  = trk.charge();
  }

  std::vector<MPlexSS> msErr(Nhits, Ntracks);
  std::vector<MPlexMV> msPar(Nhits, Ntracks);

  for (int hi = 0; hi < Nhits; ++hi)
  {
    for (int itrack=0;itrack<Ntracks;++itrack)
    {
      Track &trk = simtracks[itrack];
      Hit   &hit = trk.hitsVector()[hi];

      msErr[hi].Assign(itrack, hit.error().Array());
      msPar[hi].Assign(itrack, hit.parameters().Array());
    }
  }

  updateParametersContext updateCtx(Ntracks);

  double time = dtime();

  for (int hi = 0; hi < Nhits; ++hi)
  {
    //for (int itrack=0;itrack<Ntracks;++itrack)
#ifdef COMMENT_OUT
    {
      Track &trk = simtracks[itrack];
      Hit   &hit = trk.hitsVector()[hi];

    // std::cout << std::endl;
    // std::cout << "processing track #" << itrack << std::endl;
    
    // std::cout << "init x: " << trk.parameters()[0] << " " << trk.parameters()[1] << " " << trk.parameters()[2] << std::endl;
    // std::cout << "init p: " << trk.parameters()[3] << " " << trk.parameters()[4] << " " << trk.parameters()[5] << std::endl;
    // std::cout << "init e: " << std::endl;
    // dumpMatrix(trk.errors());

      TrackState       updatedState;
      outErr.SetArray(itrack, updatedState.errors.Array());
      outPar.SetArray(itrack, updatedState.parameters.Array());
      updatedState.charge = trk.charge();

	// std::cout << "updatedState" << std::endl;
	// std::cout << "x: " << updatedState.parameters[0] << " " << updatedState.parameters[1] << " " << updatedState.parameters[2] << std::endl;
	// std::cout << "p: " << updatedState.parameters[3] << " " << updatedState.parameters[4] << " " << updatedState.parameters[5] << std::endl;
	// std::cout << "updatedState.errors" << std::endl;
	// dumpMatrix(updatedState.errors);	

      float x = msPar[hi].At(0, 0, itrack);
      float y = msPar[hi].At(1, 0, itrack);
      float r = sqrt(x*x + y*y);

      TrackState       propState = propagateHelixToR(updatedState, r);

	// std::cout << "propState.parameters (helix propagation)" << std::endl;
	// std::cout << "x: " << propState.parameters[0] << " " << propState.parameters[1] << " " << propState.parameters[2] << std::endl;
	// std::cout << "p: " << propState.parameters[3] << " " << propState.parameters[4] << " " << propState.parameters[5] << std::endl;
	// std::cout << "propState.errors" << std::endl;
	// dumpMatrix(propState.errors);
		
	// std::cout << "measState.parameters" << std::endl;
	// std::cout << "x: " << measState.parameters[0] << " " << measState.parameters[1] << " " << measState.parameters[2] << std::endl;
	// std::cout << "measState.errors" << std::endl;
	// dumpMatrix(measState.errors);

      psErr.Assign(itrack, propState.errors.Array());
      psPar.Assign(itrack, propState.parameters.Array());

    }
#endif

    // XXXX Note, charge is not passed (line propagation). Could be part of ctxt, too.

    propagateLineToRMPlex(Err[iC], Par[iC], msErr[hi], msPar[hi],
                           Err[iP], Par[iP],
                           updateCtx);

    updateParametersMPlex(Err[iP], Par[iP], msErr[hi], msPar[hi],
                          Err[iC], Par[iC],
                          updateCtx);

    // TrackState  updatedState;
    // outErr.SetArray(itrack, updatedState.errors.Array());
    // outPar.SetArray(itrack, updatedState.parameters.Array());

    // bool dump = false;
    //   if (dump) {
    //     std::cout << std::endl;
    //     std::cout << "processing hit #" << hi << std::endl;
	
    //     std::cout << "propState.parameters (helix propagation)" << std::endl;
    //     std::cout << "x: " << propState.parameters[0] << " " << propState.parameters[1] << " " << propState.parameters[2] << std::endl;
    //     std::cout << "p: " << propState.parameters[3] << " " << propState.parameters[4] << " " << propState.parameters[5] << std::endl;
    //     std::cout << "propState.errors" << std::endl;
    //     dumpMatrix(propState.errors);
		
    //     std::cout << "measState.parameters" << std::endl;
    //     std::cout << "x: " << measState.parameters[0] << " " << measState.parameters[1] << " " << measState.parameters[2] << std::endl;
    //     std::cout << "measState.errors" << std::endl;
    //     dumpMatrix(measState.errors);

	
    //     std::cout << "updatedState" << std::endl;
    //     std::cout << "x: " << updatedState.parameters[0] << " " << updatedState.parameters[1] << " " << updatedState.parameters[2] << std::endl;
    //     std::cout << "p: " << updatedState.parameters[3] << " " << updatedState.parameters[4] << " " << updatedState.parameters[5] << std::endl;
    //     std::cout << "updatedState.errors" << std::endl;
    //     dumpMatrix(updatedState.errors);	
    //   }

    // }

  }
    
#ifndef NO_ROOT
  if (saveTree)
  {
    for (int itrack=0;itrack<Ntracks;++itrack)
    {
      TrackState &initState = initStateV[itrack];
      TrackState  updatedState;
      outErr.SetArray(itrack, updatedState.errors.Array());
      outPar.SetArray(itrack, updatedState.parameters.Array());

	std::cout << "updatedState" << std::endl;
	std::cout << "x: " << updatedState.parameters[0] << " " << updatedState.parameters[1] << " " << updatedState.parameters[2] << std::endl;
	std::cout << "p: " << updatedState.parameters[3] << " " << updatedState.parameters[4] << " " << updatedState.parameters[5] << std::endl;
	std::cout << "updatedState.errors" << std::endl;
	dumpMatrix(updatedState.errors);	
 
      pt_mc = sqrt(initState.parameters[3]*initState.parameters[3]+initState.parameters[4]*initState.parameters[4]);
      pt_fit = sqrt(updatedState.parameters[3]*updatedState.parameters[3]+updatedState.parameters[4]*updatedState.parameters[4]);
      pt_err = sqrt( updatedState.errors[3][3]*updatedState.parameters[3]*updatedState.parameters[3] +
		     updatedState.errors[4][4]*updatedState.parameters[4]*updatedState.parameters[4] + 
		     2*updatedState.errors[3][4]*updatedState.parameters[3]*updatedState.parameters[4] )/pt_fit;
      tree->Fill();
    }

    f->Write();
    f->Close();
  }

  return dtime() - time;

#endif
}

#endif
