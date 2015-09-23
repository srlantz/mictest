#include "testPropagation.h"
#include <iostream>

#include "PropagationMPlex.h"
#include "MkFitter.h"

#include <omp.h>

void testPropagation(Event& ev) {

  std::cout << "testPropagation" << std::endl;

  std::vector<Track>& simtracks = ev.simTracks_;
  
  MkFitter *mkfp = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(10);
  
  for (int itrack = 0; itrack < simtracks.size(); itrack += NN)
    {
      int end = std::min(itrack + NN, int(simtracks.size()));
      
      mkfp->InputTracksAndHits(simtracks, ev.layerHits_, itrack, end);

      mkfp->TestPropagation();      
    }
  
  return;
  
}
