#include <cmath>

#include "Simulation.h"
#include "Debug.h"

#include "MultipleScattering.h"

//#define SOLID_SMEAR

void setupTrackByToyMC(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, HitVec& hits, unsigned int itrack,
                       int& charge, float pt, const Geometry& geom, HitVec& initHits)
{
#ifdef DEBUG
  bool debug = false;
#endif
  //assume beam spot width 1mm in xy and 1cm in z
  pos=SVector3(0.1*g_gaus(g_gen), 0.1*g_gaus(g_gen), 1.0*g_gaus(g_gen));

  dprint("pos x=" << pos[0] << " y=" << pos[1] << " z=" << pos[2]);

  if (charge==0) {
    if (g_unif(g_gen) > 0.5) charge = -1;
    else charge = 1;
  }

  //float phi = 0.5*TMath::Pi()*(1-g_unif(g_gen)); // make an angle between 0 and pi/2 //fixme
  float phi = 0.5*TMath::Pi()*g_unif(g_gen); // make an angle between 0 and pi/2
  float px = pt * cos(phi);
  float py = pt * sin(phi);

  if (g_unif(g_gen)>0.5) px*=-1.;
  if (g_unif(g_gen)>0.5) py*=-1.;
  float pz = pt*(2.3*(g_unif(g_gen)-0.5));//so that we have -1<eta<1

  mom=SVector3(px,py,pz);
  covtrk=ROOT::Math::SMatrixIdentity();
  //initial covariance can be tricky
  for (unsigned int r=0; r<6; ++r) {
    for (unsigned int c=0; c<6; ++c) {
      if (r==c) {
      if (r<3) covtrk(r,c)=pow(1.0*pos[r],2);//100% uncertainty on position
      else covtrk(r,c)=pow(1.0*mom[r-3],2);  //100% uncertainty on momentum
      } else {
        covtrk(r,c)=0.;                   //no covariance
      }
    }
  }

  dprint("track with px: " << px << " py: " << py << " pz: " << pz << " pt: " << sqrt(px*px+py*py) << " p: " << sqrt(px*px+py*py+pz*pz) << std::endl);

  const float hitposerrXY = 0.01;//assume 100mum uncertainty in xy coordinate
  const float hitposerrZ = 0.1;//assume 1mm uncertainty in z coordinate
  const float hitposerrR = hitposerrXY/10.;

  const float varXY  = hitposerrXY*hitposerrXY;
  const float varZ   = hitposerrZ*hitposerrZ;

  TrackState initState;
  initState.parameters=SVector6(pos[0],pos[1],pos[2],mom[0],mom[1],mom[2]);
  initState.errors=covtrk;
  initState.charge=charge;

  TrackState tmpState = initState;

  // useful info for loopers/overlaps

  unsigned int nLayers = geom.CountLayers();
  unsigned int layer_counts[nLayers];
  for (unsigned int ilayer=0;ilayer<nLayers;++ilayer){
    layer_counts[ilayer]=0;
  }

  unsigned int nTotHit = nLayers; // can tune this number!
  // to include loopers, and would rather add a break on the code if layer ten exceeded
  // if block BREAK if hit.Layer == theGeom->CountLayers() 
  // else --> if (NMAX TO LOOPER (maybe same as 10?) {break;} else {continue;}
  
  unsigned int simLayer = 0;

  hits.reserve(nTotHit);
  initHits.reserve(nTotHit);

  for (unsigned int ihit=0;ihit<nTotHit;++ihit) {  // go to first layer in radius using propagation.h
    //TrackState propState = propagateHelixToR(tmpState,4.*float(ihit+1));//radius of 4*ihit
    auto propState = propagateHelixToNextSolid(tmpState,geom);

    float initX   = propState.parameters.At(0);
    float initY   = propState.parameters.At(1);
    float initZ   = propState.parameters.At(2);

    UVector3 init_point(initX,initY,initZ);
    simLayer = geom.LayerIndex(init_point);

#ifdef SCATTERING
    const float z1 = g_gaus(g_gen);
    const float z2 = g_gaus(g_gen);
    const float phismear = g_unif(g_gen)*TMath::TwoPi(); // random rotation of scattering plane

    SVector6 deltaPar;
    applyMultipleScattering(geom, propState, z1, z2, phismear, deltaPar);

    propState.parameters += deltaPar;
    
    // PW END

    const float scatteredZ   = propState.parameters[2];
    const float scatteredPhi = atan2(propState.parameters[1],propState.parameters[0]);
    const float scatteredRad = sqrt(propState.parameters[0]*propState.parameters[0]+propState.parameters[1]*propState.parameters[1]);

    float hitZ    = hitposerrZ*g_gaus(g_gen)+scatteredZ;
    float hitPhi  = ((hitposerrXY/scatteredRad)*g_gaus(g_gen))+scatteredPhi;
    float hitRad  = (hitposerrR)*g_gaus(g_gen)+scatteredRad;
#else // SCATTERING
    float initPhi = atan2(initY,initX);
    float initRad = sqrt(initX*initX+initY*initY);
    float hitZ    = hitposerrZ*g_gaus(g_gen)+initZ;
    float hitPhi  = ((hitposerrXY/initRad)*g_gaus(g_gen))+initPhi;
    float hitRad  = (hitposerrR)*g_gaus(g_gen)+initRad;
#endif // SCATTERING

#ifdef SOLID_SMEAR
    UVector3 scattered_point(scatteredX,scatteredY,scatteredZ);
    const auto theScatteredSolid = geom.InsideWhat(scattered_point);
    if ( ! theScatteredSolid ) {
      std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find solid AFTER scatter." << std::endl;
      std::cerr << "itrack = " << itrack << ", ihit = " << ihit << ", r = " << sqrt(scatteredX*scatteredX + scatteredY*scatteredY) << ", r*4cm = " << 4*ihit << ", phi = " << atan2(scatteredY,scatteredX) << std::endl;
      std::cerr << "initX = " << initX << ", initY = " << initY << ", initZ = " << initZ << std::endl;
      std::cerr << "scatteredX = " << scatteredX << ", scatteredY = " << scatteredY << ", scatteredZ = " << scatteredZ << std::endl << std::endl;       
      //    continue;
    }

    UVector3 scattered_xprime; // normal on surface
    bool scattered_good = theScatteredSolid->Normal(scattered_point, scattered_xprime);
      
    if ( ! scattered_good ) {
      std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find normal vector at " << scattered_point <<std::endl;
    }
    assert(std::abs(scattered_xprime[2])<1e-10); // in this geometry normal vector is in xy plane
    
    // smear along scattered yprime (where yprime = z_prime x x_prime, zprime = (0,0,1)

    float xyres_smear = g_gaus(g_gen)*hitposerrXY;

    float hitX = scatteredX - (scattered_xprime[1] * xyres_smear); 
    float hitY = scatteredY + (scattered_xprime[0] * xyres_smear); 
#else
    float hitRad2 = hitRad*hitRad;
    float hitX    = hitRad*cos(hitPhi);
    float hitY    = hitRad*sin(hitPhi);

    float varPhi = varXY/hitRad2;
    float varR   = hitposerrR*hitposerrR;
#endif

#ifdef DEBUG
    if (debug) {
      UVector3 hit_point(hitX,hitY,hitZ);
      const auto theHitSolid = geom.InsideWhat(hit_point);
      if ( ! theHitSolid ) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find solid AFTER scatter+smear." << std::endl;
        std::cerr << "itrack = " << itrack << ", ihit = " << ihit << ", r = " << sqrt(hitX*hitX + hitY*hitY) << ", r*4cm = " << 4*ihit << ", phi = " << atan2(hitY,hitX) << std::endl;
        std::cerr << "initX = " << initX << ", initY = " << initY << ", initZ = " << initZ << std::endl;
#ifdef SCATTERING
        std::cerr << "scatteredX = " << scatteredX << ", scatteredY = " << scatteredY << ", scatteredZ = " << scatteredZ << std::endl;
#endif
        std::cerr << "hitX = " << hitX << ", hitY = " << hitY << ", hitZ = " << hitZ << std::endl << std::endl; 
      }
    }
#endif

    SVector3 x1(hitX,hitY,hitZ);
    SMatrixSym33 covXYZ = ROOT::Math::SMatrixIdentity();

#ifdef SOLID_SMEAR
    covXYZ(0,0) = varXY; // yn^2 / (xn^2 + yn^2) * delx^2 + xn^2 / (xn^2 + yn^2) * dely^2
    covXYZ(1,1) = varXY; // xn^2 / (xn^2 + yn^2) * delx^2 + yn^2 / (xn^2 + yn^2) * dely^2
    // covXYZ(0,1) -> -(xn * yn) / (xn^2 + yn^2) * delx^2 + (xn * yn) / (xn^2 + yn^2) * dely^2 
    // covXYZ(1,0)  = covXYZ(0,1)    
    covXYZ(2,2) = varZ;
#else
    covXYZ(0,0) = hitX*hitX*varR/hitRad2 + hitY*hitY*varPhi;
    covXYZ(1,1) = hitX*hitX*varPhi + hitY*hitY*varR/hitRad2;
    covXYZ(2,2) = varZ;
    covXYZ(0,1) = hitX*hitY*(varR/hitRad2 - varPhi);
    covXYZ(1,0) = covXYZ(0,1);

    dprint("initPhi: " << initPhi << " hitPhi: " << hitPhi << " initRad: " << initRad  << " hitRad: " << hitRad << std::endl
        << "initX: " << initX << " hitX: " << hitX << " initY: " << initY << " hitY: " << hitY << " initZ: " << initZ << " hitZ: " << hitZ << std::endl 
        << "cov(0,0): " << covXYZ(0,0) << " cov(1,1): " << covXYZ(1,1) << " varZ: " << varZ << " cov(2,2): " << covXYZ(2,2) << std::endl 
        << "cov(0,1): " << covXYZ(0,1) << " cov(1,0): " << covXYZ(1,0) << std::endl);
#endif

    SVector3 initVecXYZ(initX,initY,initZ);
    Hit initHitXYZ(initVecXYZ,covXYZ,itrack,simLayer,layer_counts[simLayer]); 
    initHits.push_back(initHitXYZ);

    Hit hit1(x1,covXYZ,initHitXYZ.mcHitInfo());
    hits.push_back(hit1);
    tmpState = propState;

    dprint("initHitId: " << initHitXYZ.hitID() << " hit1Id: " << hit1.hitID() <<std::endl
                         << "ihit: " << ihit << " layer: " << simLayer << " counts: " << layer_counts[simLayer]);

    ++layer_counts[simLayer]; // count the number of times passed into layer

  } // end loop over nHitsPerTrack
}
