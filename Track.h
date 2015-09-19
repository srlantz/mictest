#ifndef _track_
#define _track_

#include "Hit.h"
#include "Matrix.h"
#include <vector>
#include <queue>

namespace Config {
  static constexpr const unsigned int nlayers_per_seed = 3;
  static constexpr const unsigned int maxCand = 10;
  static constexpr const float chi2Cut = 15.;
  static constexpr const float nSigma = 3.;
  static constexpr const float minDPhi = 0.;
};

typedef std::pair<unsigned int,unsigned int> SimTkIDInfo;

struct TrackState
{
public:
  TrackState() : valid(true) {}
  SVector6 parameters;
  SMatrixSym66 errors;
  int8_t charge;
  bool valid;
  float r() const { return std::sqrt(parameters[0]*parameters[0] + parameters[1]*parameters[1]); }
};

class Track
{
public:
  Track(int reserve = 10) : chi2_(0.0) { hitIDs_.reserve(reserve); }
  Track(const TrackState& state, const HitIDVec& hitIDs, float chi2) : state_(state), hitIDs_(hitIDs), chi2_(chi2) {}
  Track(int charge, const SVector3& position, const SVector3& momentum, const SMatrixSym66& errors, const HitIDVec& hitIDs, float chi2) 
    : hitIDs_(hitIDs), chi2_(chi2) 
  {
    state_.charge=charge;
    state_.errors=errors;
    state_.parameters = SVector6(position.At(0),position.At(1),position.At(2),momentum.At(0),momentum.At(1),momentum.At(2));
    state_.valid = true;
  }
  Track(int charge, const SVector6& parameters, const SMatrixSym66& errors, const HitIDVec& hitIDs, float chi2)
    : hitIDs_(hitIDs), chi2_(chi2) 
  {
    state_.charge=charge;
    state_.errors=errors;
    state_.parameters = parameters;
    state_.valid = true;
  }
  Track(TrackState state, float chi2) :
    state_(state), chi2_(chi2) { hitIDs_.reserve(10); }
  Track(const TrackState& state, float chi2, int label, int nHits, const HitIDVec& hitIDs) :
      state_(state), hitIDs_(hitIDs), chi2_(chi2), label_(label) {}
  
  ~Track(){}

  const SVector6&     parameters() const {return state_.parameters;}
  const SMatrixSym66& errors()     const {return state_.errors;}
  const TrackState&   state()      const {return state_;}

  // Non-const versions needed for CopyOut of Matriplex.
  SVector6&     parameters_nc() {return state_.parameters;}
  SMatrixSym66& errors_nc()     {return state_.errors;}
  TrackState&   state_nc()      {return state_;}

  SVector3 position() const {return SVector3(state_.parameters[0],state_.parameters[1],state_.parameters[2]);}
  SVector3 momentum() const {return SVector3(state_.parameters[3],state_.parameters[4],state_.parameters[5]);}

  int8_t   charge() const {return state_.charge;}
  float    chi2()   const {return chi2_;}
  int      label()  const {return label_;}

  float posPhi() const { return getPhi(state_.parameters[0],state_.parameters[1]); }
  float momPhi() const { return getPhi(state_.parameters[3],state_.parameters[4]); }
  float posEta() const { return getEta(state_.parameters[0],state_.parameters[1],state_.parameters[2]); }
  float momEta() const { return getEta(state_.parameters[3],state_.parameters[4],state_.parameters[5]); }

  float posR()   const { return getHypot(state_.parameters[0],state_.parameters[1]); }
  float pT()     const { return getHypot(state_.parameters[3],state_.parameters[4]); }

  const HitIDVec& hitIDs() const {return hitIDs_;}

  void addHit(HitID hitID,float chi2) {hitIDs_.push_back(hitID);chi2_+=chi2;}
  void resetHits() { hitIDs_.clear(); }
  int  nHits()   const { return hitIDs_.size(); }
  HitID getHitIdx(int hi) {
    if (hi < hitIDs_.size()) { return hitIDs_[hi]; } else { return HitID(HitID::InvHitID); }
  }
  int nFoundHits() const { return hitIDs_.size(); } // fixme
  int nTotalHits() const { return hitIDs_.size(); } // fixme

  void setCharge(int8_t chg)  {state_.charge=chg;}
  void setChi2(float chi2) {chi2_=chi2;}

  void setLabel(int lbl)   {label_=lbl;}
  void setState(TrackState newState) {state_=newState;}

  Track clone() const {return Track(state_,hitIDs_,chi2_);}
  Track clone_for_io() { return Track(state_,chi2_);}

  void write_out(FILE *fp);
  void read_in  (FILE *fp);

private:
  TrackState state_;
  HitIDVec hitIDs_;
  float chi2_;
  int label_ = -1;
};

inline bool operator<(const Track& cand1, const Track& cand2)
{
  if (cand1.nFoundHits()==cand2.nFoundHits()) return cand1.chi2()<cand2.chi2();
  return cand1.nFoundHits()>cand2.nFoundHits();
}

typedef std::vector<Track> TrackVec;

class TrackQueue : private std::priority_queue<Track> {
public:
  typedef std::priority_queue<Track> parent_type;
  using parent_type::empty;
  using parent_type::size;
  using parent_type::top;
  //using parent_type::push;
  using parent_type::emplace;
  using parent_type::pop;

  template <class... Args> reference emplace_start(Args&&... args)
  {
    if (size() == c.capacity()) pop();
    c.emplace_back(std::forward<Args>(args)...); return c.back();
  }
  void emplace_finish() { std::push_heap(c.begin(), c.end(), comp); }
  void clear() { c.clear(); }
  void swap(container_type& x) { c.swap(x); }
  void reserve(size_type n) { c.reserve(n); }
  void maybe_push(const value_type& x) {
    if (size() < c.capacity()) {
      push(x);
    } else {
      if (x < top()) {
        pop();
        push(x);
      }
    }
  }
  bool addp(float chi2, int nhits) {
    return size() < c.capacity() || nhits > top().nFoundHits() || (nhits == top().nFoundHits() && chi2 < top().chi2());
  }
};

typedef std::vector<TrackQueue> CandVec;
#endif
