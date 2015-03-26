#ifndef _hit_
#define _hit_

#include <cmath>

#include "Matrix.h"

namespace Config
{
  static constexpr const int nPhiPart   = 63;
  static constexpr const int nPhiFactor = 10;  // nPhiPart/2pi
  static constexpr const int nEtaPart   = 10;

  static const int   nEtaBin   = 2*nEtaPart - 1;
  static const float fEtaDet   = 1;
  static const float fEtaFull  = 2 * fEtaDet;
  static const float fEtaOff   = fEtaDet - fEtaDet / (2 * nEtaPart);
  static const float fEtaFac   = nEtaBin / (fEtaFull - fEtaFull / (2 * nEtaPart));

  static const float fEtaOffB1 = fEtaDet;
  static const float fEtaFacB1 = nEtaPart / fEtaFull;
  static const float fEtaOffB2 = fEtaDet - fEtaFull / (2 * nEtaPart);
  static const float fEtaFacB2 = (nEtaPart - 1) / (fEtaFull - fEtaFull / nEtaPart);

  // This is for extra bins narrower ... thinking about this some more it
  // seems it would be even better to have two more exta bins, hanging off at
  // both ends.
  //
  // Anyway, it doesn't matter ... as with wide vertex region this eta binning
  // won't make much sense -- will have to be done differently for different
  // track orgin hypotheses. In about a year or so.

  inline int getEtaBin(float eta)
  {
    int bin = (eta + fEtaOff) * fEtaFac;
    // Return -1 if outside ... this back-fitting causes trouble.
    // if      (bin < 0)        bin = 0;
    // else if (bin >= nEtaBin) bin = nEtaBin - 1;
    if (bin < 0 || bin >= nEtaBin) return -1;
    return bin;
  }

  inline int getBothEtaBins(float eta, int& b1, int& b2)
  {
    b1 = b2 = -1;

    if (eta < -fEtaDet || eta > fEtaDet)
    {
      return 0;
    }

    int b1p = std::floor((eta + fEtaOffB1) * fEtaFacB1);
    int b2p = std::floor((eta + fEtaOffB2) * fEtaFacB2);

    // printf("b1' = %d   b2' = %d\n", b1p, b2p);

    int cnt = 0;
    if (b1p >= 0 && b1p < nEtaPart)
    {
      b1 = 2 * b1p;
      ++cnt;
    }
    if (b2p >= 0 && b2p < nEtaPart - 1)
    {
      b2 = 2 * b2p + 1;
      ++cnt;
    }

    // printf("b1  = %d   b2  = %d\n", b1, b2);

    return cnt;
  }
};

class HitsOnLayer
{
  
public:
};


inline int getPhiPartition(float phi)
{
  //assume phi is between -PI and PI
  //  if (!(fabs(phi)<TMath::Pi())) std::cout << "anomalous phi=" << phi << std::endl;
  float phiPlusPi  = phi + TMath::Pi();
  int bin = phiPlusPi*10;
  return bin;
}

inline int getEtaPartition(float eta, float etaDet)
{
  float etaPlusEtaDet  = eta + etaDet;
  float twiceEtaDet    = 2.0 * etaDet;
  int bin     = (etaPlusEtaDet * Config::nEtaPart) / twiceEtaDet;  // ten bins for now ... update if changed in Event.cc
  return bin;
}

inline float getPhi(float x, float y)
{
  return std::atan2(y,x); 
}

inline float getPhiFast(float x, float y)
{
  //from http://math.stackexchange.com/questions/1098487/atan2-faster-approximation
  //need to check this is actually faster
  //need to check what is the precision
  float a = std::min(fabs(x), fabs(y)) / std::max(fabs(x), fabs(y));
  float s = a * a;
  float r = ((-0.0464964749 * s + 0.15931422) * s - 0.327622764) * s * a + a;
  if (fabs(y) > fabs(x)) r = 1.57079637 - r;
  if (x < 0) r = 3.14159274 - r;
  if (y < 0) r = -r;
  return r;
}

inline float getEta(float x, float y, float z)
{
  float theta = atan2( std::sqrt(x*x+y*y), z );
  return -1. * log( tan(theta/2.) );
}

inline float getHypot(float x, float y)
{
  return sqrtf(x*x + y*y);
}


struct MCHitInfo {
  MCHitInfo() : mcHitID_(++mcHitIDCounter_) {}
  MCHitInfo(int track, int layer, int ithlayerhit)
    : mcTrackID_(track), layer_(layer), ithLayerHit_(ithlayerhit), mcHitID_(++mcHitIDCounter_) {}

  int mcTrackID_;
  int layer_;
  int ithLayerHit_;
  int mcHitID_;

  static int mcHitIDCounter_;
};

struct MeasurementState
{
public:
  SVector3 parameters;
  SMatrixSym33 errors;
};

class Hit
{
public:
  Hit(){}

  Hit(MeasurementState state)
  {
    state_=state;
  }

  Hit(SVector3 position, SMatrixSym33 error)
  {
    state_.parameters=position;
    state_.errors=error;
  }

  Hit(SVector3 position, SMatrixSym33 error, int itrack, int ilayer, int ithLayerHit){
    mcHitInfo_.mcTrackID_ = itrack;
    mcHitInfo_.layer_ = ilayer;
    mcHitInfo_.ithLayerHit_ = ithLayerHit;
    state_.parameters=position;
    state_.errors=error;
  }

  Hit(SVector3 position, SMatrixSym33 error, const MCHitInfo& mcHitInfo)
    : mcHitInfo_(mcHitInfo)
  {
    state_.parameters=position;
    state_.errors=error;
  }

  ~Hit(){}

  SVector3&  position() {return state_.parameters;}
  SMatrixSym33& error() {return state_.errors;}
  SVector3& parameters() {return state_.parameters;}
  const SVector3&  position() const {return state_.parameters;}
  const SMatrixSym33& error() const {return state_.errors;}
  const SVector3& parameters() const {return state_.parameters;}

  float r() {
    return sqrt(state_.parameters.At(0)*state_.parameters.At(0) +
                state_.parameters.At(1)*state_.parameters.At(1));
  }
  float z() const {
    return state_.parameters.At(2);
  }
  float phi() const{
    return getPhi(state_.parameters.At(0),state_.parameters.At(1));
  }
  float eta() const{
    return getEta(state_.parameters.At(0),state_.parameters.At(1),state_.parameters.At(2));
  }

  MeasurementState measurementState() {
    return state_;
  }

  const MCHitInfo& mcHitInfo() const {return mcHitInfo_;}
  int mcTrackID() const {return mcHitInfo_.mcTrackID_;}
  int layer() const {return mcHitInfo_.layer_;}
  int ithLayerHit() const {return mcHitInfo_.ithLayerHit_;}
  int hitID() const {return mcHitInfo_.mcHitID_;}

private:
  MeasurementState state_;
  MCHitInfo mcHitInfo_;
};

typedef std::vector<Hit> HitVec;

#endif
