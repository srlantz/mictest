#include <cmath>

#include "Simulation.h"

void setupTrackByToyMC(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, std::vector<Hit>& hits, int& charge, float pt) {

  unsigned int nTotHit = 10;

  //assume beam spot width 1mm in xy and 1cm in z
  pos=SVector3(0.1*g_gaus(g_gen), 0.1*g_gaus(g_gen), 1.0*g_gaus(g_gen));

  //std::cout << "pos x=" << pos[0] << " y=" << pos[1] << " z=" << pos[2] << std::endl;

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
  for (unsigned int r=0;r<6;++r) {
    for (unsigned int c=0;c<6;++c) {
      if (r==c) {
	if (r<3) covtrk(r,c)=pow(1.0*pos[r],2);//100% uncertainty on position
	else covtrk(r,c)=pow(1.0*mom[r-3],2);  //100% uncertainty on momentum
      } else covtrk(r,c)=0.;                   //no covariance
    }
  }

  //std::cout << "track with p=" << px << " " << py << " " << pz << " pt=" << sqrt(px*px+py*py) << " p=" << sqrt(px*px+py*py+pz*pz) << std::endl;

  float hitposerrXY = 0.01;//assume 100mum uncertainty in xy coordinate//fixme
  float hitposerrZ = 0.1;//assume 1mm uncertainty in z coordinate

  TrackState initState;
  initState.parameters=SVector6(pos[0],pos[1],pos[2],mom[0],mom[1],mom[2]);
  initState.errors=covtrk;
  initState.charge=charge;

  TrackState tmpState = initState;

  //do 4 cm in radius using propagation.h
  for (unsigned int nhit=1;nhit<=nTotHit;++nhit) {

    //GIUSEPPE
    TrackState propState = propagateHelixToR_old(tmpState,4.*float(nhit));//radius of 4*nhit

    //if (nhit==1) std::cout << "crossing #" << nhit << " " << propState.parameters.At(0) << " " << propState.parameters.At(1) << " " << propState.parameters.At(2) << std::endl;

    float hitx = hitposerrXY*g_gaus(g_gen)+propState.parameters.At(0);
    float hity = hitposerrXY*g_gaus(g_gen)+propState.parameters.At(1);
    //float hity = sqrt((pos.At(0) + k*(px*sinAP-py*(1-cosAP)))*(pos.At(0) + k*(px*sinAP-py*(1-cosAP)))+
    //          	(pos.At(1) + k*(py*sinAP+px*(1-cosAP)))*(pos.At(1) + k*(py*sinAP+px*(1-cosAP)))-
    //	   	        hitx*hitx);//try to get the fixed radius
    float hitz = hitposerrZ*g_gaus(g_gen)+propState.parameters.At(2);

    //if (nhit==1) std::cout << "hit#" << nhit << " " << hitx << " " << hity << " " << hitz << std::endl;
    SVector3 x1(hitx,hity,hitz);
    SMatrixSym33 covx1 = ROOT::Math::SMatrixIdentity();
    covx1(0,0)=hitposerrXY*hitposerrXY; 
    covx1(1,1)=hitposerrXY*hitposerrXY;
    covx1(2,2)=hitposerrZ*hitposerrZ;
    Hit hit1(x1,covx1);    
    hits.push_back(hit1);  
    tmpState = propState;

    /*
    //KEVIN
    TrackState propState = propagateHelixToR(tmpState,4.*float(nhit));//radius of 4*ihit
    float initX   = propState.parameters.At(0);
    float initY   = propState.parameters.At(1);
    float initZ   = propState.parameters.At(2);
    float initPhi = atan2(initY,initX);
    float rad     = sqrt(initX*initX+initY*initY);
    float rad2    = rad*rad;

    float radfix  = rad;//fixme rad 4;

    float hitZ = hitposerrZ*g_gaus(g_gen)+initZ;
    // Comment these lines for phi smearing if using xy smear
    //float hitX = hitposerrXY*g_gaus(g_gen)+initX;
    //float hitY = hitposerrXY*g_gaus(g_gen)+initY;
    //float hitPhi = atan2(hitY,hitX);
    //   //float hity = sqrt((pos.At(0) + k*(px*sinAP-py*(1-cosAP)))*(pos.At(0) + k*(px*sinAP-py*(1-cosAP)))+
    //          (pos.At(1) + k*(py*sinAP+px*(1-cosAP)))*(pos.At(1) + k*(py*sinAP+px*(1-cosAP)))-
    //           hitx*hitx);//try to get the fixed radius
    //std::cout << "hit#" << nhit << " " << hitx << " " << hity << " " << hitz << std::endl;
    // Comment out these lines if not using xy smear

    float hitPhi  = ((hitposerrXY/radfix)*g_gaus(g_gen))+initPhi;
    //float rand_const = g_gaus(g_gen);
    //rand_const+=hitPhi;
    float hitX    = rad*cos(hitPhi);
    float hitY    = rad*sin(hitPhi);

    float radHit  = sqrt(hitX*hitX+hitY*hitY);
    float radHit2 = radHit*radHit;

    SVector3 vecXYZ(hitX,hitY,hitZ);

    float varXY = hitposerrXY * hitposerrXY;
    float varZ  = hitposerrZ * hitposerrZ;
    SMatrixSym33 covXYZ = ROOT::Math::SMatrixIdentity();
    //fixme
    covXYZ(0,0) = varXY*(hitY*hitY)/rad2;
    covXYZ(1,1) = varXY*(hitX*hitX)/rad2;
    covXYZ(2,2) = varZ;
    covXYZ(0,1) = -varXY*(hitX*hitY)/rad2/5.;//fixme
    covXYZ(1,0) = covXYZ(0,1);

    // covXYZ(0,0) = varXY*(initY*initY)/rad2;
    // covXYZ(1,1) = varXY*(initX*initX)/rad2;
    // covXYZ(2,2) = varZ;
    // // covXYZ(0,1) = -varXY*(initX*initY)/rad2;//fixme
    // // covXYZ(1,0) = covXYZ(0,1);

    //covXYZ(0,0) = varXY;
    //covXYZ(1,1) = varXY;
    //covXYZ(2,2) = varZ;

    
    //if (dump) 
    std::cout << "rad:   " << rad   << " rad2: " << rad2 << " radH:  " << radHit << " radH2: " << radHit2 << std::endl
	      << "initPhi: " << initPhi << " hitPhi: " << hitPhi << std::endl
	      << "initX: " << initX << " hitX: " << hitX << " initY: " << initY << " hitY: " << hitY << " initZ: " << initZ << " hitZ: " << hitZ << std::endl 
	      << "varXY: " << varXY << " cov(0,0): " << covXYZ(0,0) << " cov(1,1): " << covXYZ(1,1) << " varZ: " << varZ << " cov(2,2): " << covXYZ(2,2) << std::endl 
	      << "cov(0,1): " << covXYZ(0,1) << " cov(1,0): " << covXYZ(1,0) << std::endl;
    dumpMatrix(covXYZ);
    std::cout << std::endl;
    
    Hit hitXYZ(vecXYZ,covXYZ);    
    hits.push_back(hitXYZ);  

    // SVector3 initVecXYZ(initX,initY,initZ);
    // Hit initHitXYZ(initVecXYZ,covXYZ);
    // initHits.push_back(initHitXYZ);
    
    tmpState = propState;
    */
  }
  
  /*
  //do 4 cm along path
  for (unsigned int nhit=1;nhit<=nTotHit;++nhit) {
    float distance = 4.*float(nhit);//~4 cm distance along curvature between each hit
    float angPath = distance/curvature;
    float cosAP=cos(angPath);
    float sinAP=sin(angPath);
    float hitx = gRandom->Gaus(0,hitposerr)+(pos.At(0) + k*(px*sinAP-py*(1-cosAP)));
    float hity = gRandom->Gaus(0,hitposerr)+(pos.At(1) + k*(py*sinAP+px*(1-cosAP)));
    //float hity = sqrt((pos.At(0) + k*(px*sinAP-py*(1-cosAP)))*(pos.At(0) + k*(px*sinAP-py*(1-cosAP)))+
    //          	(pos.At(1) + k*(py*sinAP+px*(1-cosAP)))*(pos.At(1) + k*(py*sinAP+px*(1-cosAP)))-
    //	   	        hitx*hitx);//try to get the fixed radius
    float hitz = gRandom->Gaus(0,hitposerr)+(pos.At(2) + distance*ctgTheta);    
    //std::cout << "hit#" << nhit << " " << hitx << " " << hity << " " << hitz << std::endl;
    SVector3 x1(hitx,hity,hitz);
    SMatrixSym33 covx1 = ROOT::Math::SMatrixIdentity();
    covx1(0,0)=hitposerr*hitposerr; 
    covx1(1,1)=hitposerr*hitposerr;
    covx1(2,2)=hitposerr*hitposerr;
    Hit hit1(x1,covx1);    
    hits.push_back(hit1);
  }
  */

}
