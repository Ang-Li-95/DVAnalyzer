#ifndef DVAnalyzer_TrackRescaler_H
#define DVAnalyzer_TrackRescaler_H

#include "DataFormats/TrackReco/interface/Track.h"

// borrowed from https://github.com/jordantucker/cmssw-usercode/blob/76e289e858e08c43d1de9d24ad24abd30e53a78e/Tools/src/TrackRescaler.cc
// currently only using JetHT2018B
// just set era and which to 0 since there are not used :)
namespace DVAna{
  class TrackRescaler {
    struct scales_t {
      scales_t() : dxyerr(1), dszerr(1), dxydszcov(1) {}
      double dxyerr;
      double dszerr;
      double dxydszcov;
    } scales_;

    bool enable_;
    int era_;
    int which_; 
  public:
    TrackRescaler() : enable_(false), era_(0), which_(w_JetHT) {}
    void setup(bool enable, int era, int which) { enable_ = enable; era_ = era; which_ = which; }

    enum { w_JetHT, w_max };

    void set_JetHT2017B(double pt, double eta);
    void set_JetHT2017C(double pt, double eta);
    void set_JetHT2017DE(double pt, double eta);
    void set_JetHT2017F(double pt, double eta);
    void set_JetHT2018A(double pt, double eta);
    void set_JetHT2018B(double pt, double eta);
    void set_JetHT2018C(double pt, double eta);
    void set_JetHT2018D(double pt, double eta);

    void set(double era, int which, double pt, double eta);

    struct ret_t {
      reco::Track tk;
      reco::Track rescaled_tk;
    };

    ret_t scale(const reco::Track& tk);
  };
}

#endif
