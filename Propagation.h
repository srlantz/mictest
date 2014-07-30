#ifndef _propagation_
#define _propagation_

#include "Track.h"
#include "Matrix.h"

// line propagation from state radius to hit radius
// assuming radial direction (i.e. origin at (0,0))
TrackState propagateLineToR(TrackState& inputState, float r);

// helix propagation in steps along helix trajectory. 
// each step travels for a path lenght equal to delta r between the current position and the target radius. 
// for track with pT>=1 GeV this converges to the correct path lenght in <5 iterations
// derivatives need to be updated at each iteration
TrackState propagateHelixToR(TrackState& inputState, float r);

//test towards a helix propagation without iterative approach
//version below solves the equation for the angular path at which x^2+y^2=r^2
//problems: 1. need first order approximation of sin and cos, 
//2. there are 2 numerical solutions, 3. need to propagate uncertainties throgh the 2nd order equation
TrackState propagateHelixToR_test(TrackState& inputState, float r);



class updateParametersContext;

void propagateLineToRMPlex(const MPlexSS &psErr,  const MPlexMV& psPar,
                           const MPlexSS &msErr,  const MPlexMV& msPar,
                                 MPlexSS &outErr,       MPlexMV& outPar,
                                 updateParametersContext &ctx);

#endif
