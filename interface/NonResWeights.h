#ifndef _NonResWeights_H
#define _NonResWeights_H

#include <iostream>
#include <string>
#include "TFile.h"
#include "TH2F.h"

class NonResWeights {
 public:
  NonResWeights();
  virtual ~NonResWeights();

  void LoadHists(std::string f1, std::string f2);
  float GetWeight(UInt_t n, Float_t m, Float_t c);
    
 private:
  TFile *_w1507grid, *_w12bench;
  TH2F *_NR_Wei_Hists[1519];

};


#endif
