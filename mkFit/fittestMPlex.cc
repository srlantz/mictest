#include "fittest.h"

#include "Matrix.h"
#include "KalmanUtils.h"
#include "Propagation.h"
#include "Simulation.h"
#include "Geometry.h"

#include "MkFitter.h"

#ifndef NO_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

// #include <omp.h>
namespace
{
  int omp_get_thread_num() { return 0; }
  int omp_get_num_threads() { return 1; }
}

#include <iostream>
#include <memory>

//==============================================================================
void initGeom(Geometry& geom)
{
  std::cout << "Constructing SimpleGeometry Cylinder geometry" << std::endl;

  // NB: we currently assume that each node is a layer, and that layers
  // are added starting from the center
  // NB: z is just a dummy variable, VUSolid is actually infinite in size.  *** Therefore, set it to the eta of simulation ***
  float eta = 2.0; // can tune this to whatever geometry required (one can make this layer dependent as well)
  for (int l = 0; l < 10; l++) {
    float r = (l+1)*4.;
    VUSolid* utub = new VUSolid(r, r+.01);
    float z = r / std::tan(2.0*std::atan(std::exp(-eta))); // calculate z extent based on eta, r
    geom.AddLayer(utub, r, z);
  }
}


void generateTracks(std::vector<Track>& simtracks, int Ntracks)
{

  Geometry geom;
  initGeom(geom);

   g_gen.seed(742331);

   simtracks.resize(Ntracks);

   // double tstart = dtime();

#pragma omp parallel for num_threads(NUM_THREADS_SIM)
   for (int itrack = 0; itrack < Ntracks; ++itrack)
   {
      //create the track
      SVector3 pos;
      SVector3 mom;
      SMatrixSym66 covtrk;
      std::vector<Hit> hits;

      // Fixed pT / charge for testing
      // int   q  = 1;
      // float pt = 9.99;

      int   q  = 0;                         // set it in setup function
      float pt = 0.5 + g_unif(g_gen) * 9.5; // this input, 0.5<pt<10 GeV (below ~0.5 GeV does not make 10 layers)

      setupTrackByToyMC(pos,mom,covtrk,hits,itrack,q,pt,&geom);
      Track simtrk(q,pos,mom,covtrk,hits,0.);
      simtracks[itrack] = simtrk;
      simtracks[itrack].setLabel(itrack);
   }

   // printf("gen time = %f\n", dtime() - tstart);
}

void make_validation_tree(const char         *fname,
                          std::vector<Track> &simtracks,
                          std::vector<Track> &rectracks)
{
#ifndef NO_ROOT

   assert(simtracks.size() == rectracks.size());

   float pt_mc, pt_fit, pt_err, chg;

   TFile *file = TFile::Open(fname, "recreate");
   TTree *tree = new TTree("T", "Validation Tree for simple Kalman fitter");;

   tree->Branch("pt_mc",  &pt_mc,  "pt_mc");
   tree->Branch("pt_fit", &pt_fit, "pt_fit");
   tree->Branch("pt_err", &pt_err, "pt_err");
   tree->Branch("chg",    &chg,    "chg");

   const int NT = simtracks.size();
   for (int i = 0; i < NT; ++i)
   {
      SVector6     &simp   = simtracks[i].parameters();
      SVector6     &recp   = rectracks[i].parameters();
      SMatrixSym66 &recerr = rectracks[i].errors();

      pt_mc  = sqrt(simp[3]*simp[3] + simp[4]*simp[4]);
      pt_fit = sqrt(recp[3]*recp[3] + recp[4]*recp[4]);
      pt_err = sqrt(recerr[3][3]*recp[3]*recp[3] +
                    recerr[4][4]*recp[4]*recp[4] +
                    recerr[3][4]*recp[3]*recp[4] * 2) / pt_fit;
      chg = simtracks[i].charge();

      tree->Fill();
   }

   file->Write();
   file->Close();
   delete file;
#endif
}

//==============================================================================
// runFittingTest
//==============================================================================

double runFittingTest(std::vector<Track>& simtracks, std::vector<Track>& rectracks)
{
   // Standard fitting test using SMatrix

   //these matrices are dummy and can be optimized without multriplying by zero all the world...
   SMatrix36 projMatrix36;
   projMatrix36(0,0)=1.;
   projMatrix36(1,1)=1.;
   projMatrix36(2,2)=1.;
   SMatrix63 projMatrix36T = ROOT::Math::Transpose(projMatrix36);


   double time = dtime();

   for (unsigned int itrack=0;itrack<simtracks.size();++itrack)
   {
      Track trk = simtracks[itrack];

      // std::cout << std::endl;
      // std::cout << "processing track #" << itrack << std::endl;
    
      // std::cout << "init x: " << trk.parameters()[0] << " " << trk.parameters()[1] << " " << trk.parameters()[2] << std::endl;
      // std::cout << "init p: " << trk.parameters()[3] << " " << trk.parameters()[4] << " " << trk.parameters()[5] << std::endl;
      // std::cout << "init e: " << std::endl;
      // dumpMatrix(trk.errors());

      std::vector<Hit>& hits = trk.hitsVector();

      // Make a copy since initState is used at the end to fill the tree.
      // Hmmh, not any more, it seems.
      // TrackState initState = trk.state();
      TrackState& updatedState = trk.state();
    
      bool dump = false;
    
      for (std::vector<Hit>::iterator hit=hits.begin();hit!=hits.end();++hit)
      {
         //for each hit, propagate to hit radius and update track state with hit measurement
         TrackState propState = propagateHelixToR(updatedState, hit->r());

         MeasurementState measState = hit->measurementState();
         updatedState = updateParameters(propState, measState,projMatrix36,projMatrix36T);
         //updateParameters66(propState, measState, updatedState);//updated state is now modified
      
         if (dump)
         {
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
      rectracks.push_back(trk);
#endif
   }

   return dtime() - time;
}


//==============================================================================
// runFittingTestPlex
//==============================================================================

double runFittingTestPlex(std::vector<Track>& simtracks, std::vector<Track>& rectracks)
{
   const int Nhits = MAX_HITS;
   // XXX What if there's a missing / double layer?
   // Eventually, should sort track vector by number of hits!
   // And pass the number in on each "setup" call.
   // Reserves should be made for maximum possible number (but this is just
   // measurments errors, params).

   // NOTE: MkFitter *MUST* be on heap, not on stack!
   // Standard operator new screws up alignment of ALL MPlex memebrs of MkFitter,
   // even if one adds attr(aligned(64)) thingy to every possible place.

   // MkFitter *mkfp = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Nhits);

   MkFitter *mkfp_arr[NUM_THREADS];

   for (int i = 0; i < NUM_THREADS; ++i)
   {
     mkfp_arr[i] = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(Nhits);
   }

   int theEnd = simtracks.size();

   double time = dtime();

#pragma omp parallel for num_threads(NUM_THREADS)
   for (int itrack = 0; itrack < theEnd; itrack += NN)
   {
      int end = std::min(itrack + NN, theEnd);

      MkFitter *mkfp = mkfp_arr[omp_get_thread_num()];

      mkfp->InputTracksAndHits(simtracks, itrack, end);

      mkfp->FitTracks();

#ifndef NO_ROOT
      mkfp->OutputFittedTracks(rectracks, itrack, end);
#endif
   }

   time = dtime() - time;

   for (int i = 0; i < NUM_THREADS; ++i)
   {
     _mm_free(mkfp_arr[i]);
   }
   //_mm_free(mkfp);

   return time;
}
