#ifndef _simulation_
#define _simulation_

#include "Track.h"
#include "Matrix.h"
#include "Propagation.h"
#include "Geometry.h"

void setupTrackByToyMC(unsigned int itrack, float pt,
  SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, 
  MCHitVec& mcHits, HitVec& hits, int& charge, const Geometry& geom);
#endif
