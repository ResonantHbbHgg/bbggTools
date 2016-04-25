//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 15 10:49:49 2015 by ROOT version 6.02/05
// from TTree bbggSelectionTree/Flat tree for HH->bbgg analyses (after pre selection)
// found on file: bbgg2DFitter_DoubleEG_37.root
//////////////////////////////////////////////////////////

#ifndef bbgg2DFitter_h
#define bbgg2DFitter_h

// C++ headers
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <boost/program_options.hpp>
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
// RooFit headers
#include <RooWorkspace.h>
#include <RooFitResult.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooArgSet.h>
#include "RooStats/HLFactory.h"
#include <RooDataSet.h>
#include <RooFormulaVar.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooAbsPdf.h>
#include <RooBernstein.h>
#include <RooExtendPdf.h>
#include <RooMinimizer.h>
#include "RooStats/RooStatsUtils.h"
#include <RooProdPdf.h>
#include <RooExponential.h>
#include <RooPolynomial.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <cmath>
#include <Math/LorentzVector.h>
#include <algorithm>
#include <string>
#include <utility>

// namespaces
using namespace std;
using namespace RooFit;
using namespace RooStats;
namespace po = boost::program_options;

class bbgg2DFitter {
public :
//   Parameters
   Bool_t _doblinding;
   Int_t _NCAT;
   Int_t _sigMass;
   bool _addHiggs;
   float _lumi;
   TString _cut;
   float _minMgg;
   float _maxMgg;
   float _minMjj;
   float _maxMjj;
   float _minMggHig;
   float _maxMggHig;
   float _minMjjHig;
   float _maxMjjHig;
   std::string _lumiStr;
   std::string _energyStr;
   
   //Workspace
   RooWorkspace* _w;

   bbgg2DFitter(Int_t SigMass, float Lumi, Bool_t doBlinding, Int_t nCat, bool AddHiggs );// { doblinding = false; NCAT = 0; addHiggs = true; cardName = cardname; Init();}
   virtual ~bbgg2DFitter() { }
   void SetCut(TString cut) {_cut = cut;}
   RooArgSet* defineVariables(); //DONE
   void AddSigData(float mass, TString signalfile); //DONE
   void AddHigData(float mass, TString signalfile, int higgschannel); //DONE
   void AddBkgData(TString datafile); //DONE
   void SigModelFit(float mass); //DONE
   void HigModelFit(float mass, int higgschannel); //DONE
   RooFitResult* BkgModelFit(Bool_t); //DONE
   void MakePlots(float mass); //DONE
   void MakePlotsHiggs(float mass); //DONE
   void MakeSigWS(const char* filename); //DONE
   void MakeHigWS(const char* filename, int higgschannel); //DONE
   void MakeBkgWS(const char* filename); //DONE
   void MakeDataCard(const char* filename, const char* filename1,
                     const char* fileHiggsNameggh, const char* fileHiggsNametth, 
                     const char* fileHiggsNamevbf, const char* fileHiggsNamevh, 
                     const char* fileHiggsNamebbh, Bool_t); //DONE
   void MakeDataCardonecatnohiggs(const char* filename1, const char* filename2, Bool_t useSigTheoryUnc); //DONE
   void SetConstantParams(const RooArgSet* params); //DONE
   void PrintWorkspace() {_w->Print("v");}
   void SetMinMaxMggMjj(float minmgg, float maxmgg, float minmjj, float maxmjj) { _minMgg = minmgg; _maxMgg = maxmgg; _minMjj = minmjj; _maxMjj = maxmjj;}
   void SetMinMaxMggMjjHig(float minmgg, float maxmgg, float minmjj, float maxmjj) { _minMggHig = minmgg; _maxMggHig = maxmgg; _minMjjHig = minmjj; _maxMjjHig = maxmjj;}
   void SetWorkspace( RooWorkspace* workspace ) { _w = new RooWorkspace(*workspace); }
   std::string MakeModelCard( std::string file );
   void SetEnergyAndLumi( std::string energy, std::string lumi) { _lumiStr = lumi; _energyStr = energy;}
   TStyle * style(); //DONE
   
   ClassDef(bbgg2DFitter,0);
};

#endif

#ifdef bbgg2DFitter_cxx
bbgg2DFitter::bbgg2DFitter(Int_t SigMass, float Lumi, Bool_t doBlinding = false, 
                                Int_t nCat = 0, bool AddHiggs = true)
{
    _doblinding = doBlinding;
    _NCAT = nCat;
    _sigMass = SigMass;
    _addHiggs = AddHiggs;
    _lumi = Lumi;
    _cut = "1";
   _minMgg = 100.;
   _maxMgg = 180.;
   _minMjj = 60.;
   _maxMjj = 180.;
   _minMggHig = 115.;
   _maxMggHig = 135.;
   _minMjjHig = 60.;
   _maxMjjHig = 180.;
    
}

#endif // #ifdef bbgg2DFitter_cxx
