#include "MultipleScattering.h"
#include "Debug.h"

//#define SCATTER_XYZ

void applyMultipleScattering(const Geometry& geom, const TrackState& propState, const float z1, const float z2, const float phismear, SVector6& deltaPar) {

  // PW START
  // multiple scattering. Follow discussion in PDG Ch 32.3
  // this generates a scattering distance yplane and an angle theta_plane in a coordinate
  // system normal to the initial direction of the incident particle ...
  // question: do I need to generate two scattering angles in two planes (theta_plane) and add those or 
  // can I do what I did, generate one (theta_space) and rotate it by another random numbers
  // const float z1 = g_gaus(g_gen);
  // const float z2 = g_gaus(g_gen);
  // const float phismear = g_unif(g_gen)*TMath::TwoPi(); // random rotation of scattering plane
  const float X0 = 9.370; // cm, from http://pdg.lbl.gov/2014/AtomicNuclearProperties/HTML/silicon_Si.html
  //const float X0 = 0.5612; // cm, for Pb
  float x = 0.1; // cm  -assumes radial impact. This is bigger than what we have in main
  // will update for tilt down a few lines
  const float p = sqrt(propState.parameters.At(3)*propState.parameters.At(3)+
		       propState.parameters.At(4)*propState.parameters.At(4)+
		       propState.parameters.At(5)*propState.parameters.At(5));
  UVector3 pvec(propState.parameters.At(3)/p, propState.parameters.At(4)/p, propState.parameters.At(5)/p);

  // now we need to transform to the right coordinate system
  // in this coordinate system z -> z'=z i.e., unchanged (this is a freedom I have since the 
  // definition of the plane defined by nhat is invariant under rotations about nhat
  // x, y change
  // so the transformation is a rotation about z and a translation
  // if the normal vector is given by n == x' = (x1,y1,0) [NB no component along z here]
  // then the 
  // y axis in the new coordinate system is given by y' = z' x x' = (-y1,x1,0)

  const float& initX   = propState.parameters[0];
  const float& initY   = propState.parameters[1];
  const float& initZ   = propState.parameters[2];
  UVector3 init_point(initX,initY,initZ);

  const auto theInitSolid = geom.InsideWhat(init_point);
  if ( ! theInitSolid ) {
    std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find solid." <<std::endl;
    float initPhi = atan2(initY,initX);
    float initRad = sqrt(initX*initX+initY*initY);
    std::cerr << "itrack = " << /*itrack*/0 << ", ihit = " << /*ihit*/0 << ", r = " << initRad << ", r*4cm = " << /*4*ihit*/0 << ", phi = " << initPhi << std::endl;
    std::cerr << "initX = " << initX << ", initY = " << initY << ", initZ = " << initZ << std::endl;
    float pt = sqrt(propState.parameters[3]*propState.parameters[3]+
		    propState.parameters[4]*propState.parameters[4]);
    std::cerr << "pt = " << pt << ", pz = " << propState.parameters[5] << std::endl;

    return;
  }
  UVector3 init_xprime; // normal on surface
  bool init_good = theInitSolid->Normal(init_point, init_xprime);
      
  if ( ! init_good ) {
    std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find normal vector at " << init_point <<std::endl;
    return;
  }
  assert(std::abs(init_xprime[2])<1e-10); // in this geometry normal vector is in xy plane

  x /= std::abs(init_xprime.Dot(pvec)); // take care of dip angle
  const float betac = sqrt(p*p+(.135*.135))/(p*p); // m=130 MeV, pion
  const float theta_0 = 0.0136/(betac*p)*sqrt(x/X0)*(1+0.038*log(x/X0));// eq 32.14
  const float y_plane = z1*x*theta_0/sqrt(12.)+ z2*x*theta_0/2.;
  const float theta_plane = z2*theta_0;
  const float theta_space = sqrt(2)*theta_plane;
  dprint("yplane, theta_space = " << y_plane << ", " << theta_space);

  UVector3 yprime(-init_xprime[1],init_xprime[0],0); // result of cross product with zhat
  //const double phi0 = atan2(xpime[1], init_xprime[0]);

#ifdef SCATTER_XYZ
  const float scatteredX = initX + y_plane *(-init_xprime[1]*cos(phismear)); 
  const float scatteredY = initY + y_plane *(+init_xprime[0]*cos(phismear));
  const float scatteredZ = initZ + y_plane *(           sin(phismear));
#else 
  const float scatteredX = initX;
  const float scatteredY = initY;
  const float scatteredZ = initZ;
#endif  // SCATTER_XYZ
  UVector3 pvecprime;

  const float v0 = sqrt(2+pow((pvec[1]+pvec[2])/pvec[0],2));
  const float v1 = sqrt(2-pow((pvec[1]-pvec[2]),2));
  const float a = pvec[0]; 
  const float b = pvec[1]; 
  const float c = pvec[2];

  auto sign = [] (const float a) { 
    if ( a > 0 ) 
      return 1;
    else if ( a < 0 ) 
      return -1;
    else
      return 0;
  };
                    
  pvecprime[0] = a*cos(theta_space) + ((b + c)*cos(phismear)*sin(theta_space))/(a*v0) + 
    (a*(b - c)*sin(theta_space)*sin(phismear))/(v1*sign(1 + (b - c)*c));
  pvecprime[1] = b*cos(theta_space) - (cos(phismear)*sin(theta_space))/v0 + 
    ((-1 + pow(b,2) - b*c)*sin(theta_space)*sin(phismear))/(v1*sign(1 + (b - c)*c));
  pvecprime[2] = c*cos(theta_space) - (cos(phismear)*sin(theta_space))/v0 + 
    (std::abs(1 + (b - c)*c)*sin(theta_space)*sin(phismear))/v1; 
  assert(pvecprime.Mag()<=1.0001);

  dprint("theta_space, phismear = " << theta_space << ", " << phismear << std::endl
	 << "init_xprime = " << init_xprime << std::endl
	 << "yprime = " << yprime << std::endl
	 << "phat      = " << pvec << "\t" << pvec.Mag() << std::endl
	 << "pvecprime = " << pvecprime << "\t" << pvecprime.Mag() << std::endl
	 << "angle     = " << pvecprime.Dot(pvec) << "(" << cos(theta_space) << ")" << std::endl
	 << "pt, before and after: " << pvec.Perp()*p << ", "<< pvecprime.Perp()*p);

  pvecprime.Normalize();

  // now update propstate with the new scattered results
  deltaPar[0] = scatteredX - propState.parameters[0];
  deltaPar[1] = scatteredY - propState.parameters[1];
  deltaPar[2] = scatteredZ - propState.parameters[2];
  deltaPar[3] = pvecprime[0]*p - propState.parameters[3];
  deltaPar[4] = pvecprime[1]*p - propState.parameters[4];
  deltaPar[5] = pvecprime[2]*p - propState.parameters[5];
    
}
