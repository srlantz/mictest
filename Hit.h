#ifndef _hit_
#define _hit_

#include <cmath>
#include <climits>
#include <cstdint>

#include "Matrix.h"
#include <atomic>

namespace Config
{
  static constexpr float    PI = 3.14159265358979323846;
  static constexpr float TwoPI = 6.28318530717958647692;
  
  static constexpr int   nPhiPart   = 1260; // 63;
  static constexpr float nPhiFactor = nPhiPart / TwoPI;
  static constexpr int   nEtaPart   = 11; // 10;

  static const int   nEtaBin   = 2*nEtaPart - 1;
  static const float fEtaDet   = 1;
  static const float fEtaFull  = 2 * fEtaDet;
  static const float lEtaPart  = fEtaFull/float(nEtaPart);
  static const float lEtaBin   = lEtaPart/2.;

  static const float fEtaOffB1 = fEtaDet;
  static const float fEtaFacB1 = nEtaPart / fEtaFull;
  static const float fEtaOffB2 = fEtaDet - fEtaFull / (2 * nEtaPart);
  static const float fEtaFacB2 = (nEtaPart - 1) / (fEtaFull - fEtaFull / nEtaPart);
};

// This is for extra bins narrower ... thinking about this some more it
// seems it would be even better to have two more exta bins, hanging off at
// both ends.
//
// Anyway, it doesn't matter ... as with wide vertex region this eta binning
// won't make much sense -- will have to be done differently for different
// track orgin hypotheses. In about a year or so.

inline float normalizedPhi(float phi)
{
  // Return phi between -pi and +pi.
  //   return std::fmod(phi, (float) M_PI);

  while ( phi < -Config::PI ) phi += Config::TwoPI;
  while ( phi >  Config::PI ) phi -= Config::TwoPI;
  return phi;
}

inline int getEtaBin(float eta)
{
  using namespace Config;

  //in this case we are out of bounds
  if (std::abs(eta)>fEtaDet) return -1;

  //first and last bin have extra width
  if (eta<(lEtaBin-fEtaDet)) return 0;
  if (eta>(fEtaDet-lEtaBin)) return nEtaBin-1;

  //now we can treat all bins as if they had same size
  return int( (eta+fEtaDet-lEtaBin/2.)/lEtaBin );
}

inline int getBothEtaBins(float eta, int& b1, int& b2)
{
  using namespace Config;

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

inline int getPhiPartition(float phi)
{
  //assume phi is between -PI and PI
  //  if (!(fabs(phi)<TMath::Pi())) std::cout << "anomalous phi=" << phi << std::endl;
  float phiPlusPi  = phi + Config::PI;
  int bin = phiPlusPi * Config::nPhiFactor;
  return bin;
}

inline int getEtaPartition(float eta, float etaDet = Config::fEtaFull)
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

inline float getEta(float x, float y, float z)
{
  float theta = atan2( std::sqrt(x*x+y*y), z );
  return -1. * log( tan(theta/2.) );
}

inline float getHypot(float x, float y)
{
  return std::sqrt(x*x + y*y);
}

inline unsigned int binNumber(unsigned int etaPart, unsigned int phiPart)
{
 return Config::nPhiPart*etaPart + phiPart;
}

struct HitID {
  HitID() : layer_(0), index_(0) {}
  static constexpr int MCLayerID = -1;
  static constexpr int InvHitID = -1;
  HitID(int layer, int index) : layer_(layer), index_(index) {}
  HitID(int index) : layer_(MCLayerID), index_(index) {}
  int layer_;
  int index_;
};

inline bool operator==(const HitID& a, const HitID& b) {
  return a.layer_ == b.layer_ && a.index_ == b.index_;
}

struct MeasurementState
{
public:
  MeasurementState() {}
  MeasurementState(const SVector3& p, const SVector6& e)
    : pos_(p), err_(e) {}
  MeasurementState(const SVector3& p, const SMatrixSym33& e)
    : pos_(p) {
      for (int i=0;i<6;++i) err_[i] = e.Array()[i];
    }
  const SVector3& parameters() const { return pos_; }
  SMatrixSym33 errors() const { 
    SMatrixSym33 result;
    for (int i=0;i<6;++i) result.Array()[i]=err_[i];
    return result; 
  }
  SVector3 pos_;
  SVector6 err_;
};

class Hit
{
public:
  Hit() : mcHitID_(HitID::InvHitID) {}
  Hit(const MeasurementState& state) : state_(state) {}
  Hit(const SVector3& position, const SMatrixSym33& error)
    : state_(position, error) {}
  Hit(const SVector3& position, const SMatrixSym33& error, int mcHitID)
    : state_(position, error), mcHitID_(mcHitID) {}

  ~Hit(){}

  const SVector3&     position()   const {return state_.pos_;}
  const SVector3&     parameters() const {return state_.pos_;}
  const SMatrixSym33  error()      const {return state_.errors();}

  const float* posArray() const {return state_.pos_.Array();}
  const float* errArray() const {return state_.err_.Array();}

  // Non-const versions needed for CopyOut of Matriplex.
  SVector3&     parameters_nc() {return state_.pos_;}
  SVector6&     error_nc()      {return state_.err_;}

  float r() const {
    return std::sqrt(state_.parameters().At(0)*state_.parameters().At(0) +
                     state_.parameters().At(1)*state_.parameters().At(1));
  }
  float x() const {
    return state_.parameters().At(0);
  }
  float y() const {
    return state_.parameters().At(1);
  }
  float z() const {
    return state_.parameters().At(2);
  }
  float phi() const {
    return getPhi(state_.parameters().At(0), state_.parameters().At(1));
  }
  float eta() const {
    return getEta(state_.parameters().At(0), state_.parameters().At(1), state_.parameters().At(2));
  }
  int phiPart() const { return getPhiPartition(phi()); }
  int etaPart() const { return getEtaPartition(eta()); }
  int binIndex() const { return binNumber(etaPart(), phiPart()); }
  const MeasurementState& measurementState() const {
    return state_;
  }
  int mcHitID() const { return mcHitID_; }
  void setMCHitID(int id) { mcHitID_ = id; }

private:
  MeasurementState state_;
  int mcHitID_;
};

struct MCHit : public Hit
{
  MCHit() : hitID_(HitID::MCLayerID, HitID::InvHitID) {}
  MCHit(const Hit& hit, int track, int layer, int ithlayerhit)
    : Hit(hit), mcTrackID_(track), layer_(layer), ithLayerHit_(ithlayerhit)
      { setMCHitID(mcHitIDCounter_++); }
  int layer() const { return layer_; }
  HitID hitID() const { return hitID_; }
  void setHitID(HitID hitID) { hitID_ = hitID; }
  int mcTrackID() const { return mcTrackID_; }

  int mcTrackID_;
  int layer_;
  int ithLayerHit_;
  HitID hitID_;
  static std::atomic<int> mcHitIDCounter_;
};

typedef std::vector<Hit> HitVec;
typedef std::vector<MCHit> MCHitVec;
typedef std::vector<HitID> HitIDVec;

#endif
