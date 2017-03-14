#include "flashgg/bbggTools/interface/NonResWeights.h"

NonResWeights::NonResWeights(){}

NonResWeights::~NonResWeights(){}

void NonResWeights::LoadHists(std::string f1, std::string f2){

  // This file should incude the 1507 hists for full grid re-weighting
  _w1507grid = new TFile(f1.c_str(), "OPEN");
  
  // This file should contain 12 benchmark hists of weights.
  _w12bench = new TFile(f2.c_str(), "OPEN");


  if (_w1507grid->IsZombie() || _w12bench->IsZombie() ){
    std::cout<<" Input file does not exist!"<<std::endl;
    exit(1);
  }
  _w1507grid->Print();
  _w12bench->Print();

  TList *histList1 = _w1507grid->GetListOfKeys();
  for (UInt_t n=0; n<1507; n++){
    if (histList1->Contains(Form("point_%i_weights",n)))
      _NR_Wei_Hists[n] = (TH2F*)_w1507grid->Get(Form("point_%i_weights",n));
    else
      std::cout<<"This one does not existe pas: "<<n<<std::endl;
  }


  TList *histList2 = _w12bench->GetListOfKeys();
  for (UInt_t n=0; n<12; n++){
    if (histList2->Contains(Form("point_%i_weights",n)))
      _NR_Wei_Hists[1507+n] = (TH2F*)_w12bench->Get(Form("point_%i_weights",n));
    else
      std::cout<<"This one does not existe pas: "<<1507+n<<std::endl;
  }


}

float NonResWeights::GetWeight(UInt_t n, Float_t mhh, Float_t cs){
  float w = 0;
  if (n==324 || n==910 || n==985 || n==990){
    // The points above do not exist in the input file provided by Alexandra (and wont ever be added)
    //cout<<"This one was not existing in the input file: "<<n<<endl;
    w=0;
  }
  else {
    UInt_t binNum = _NR_Wei_Hists[n]->FindBin(mhh, fabs(cs));
    w = _NR_Wei_Hists[n]->GetBinContent(binNum);
  }
  return w;
}
