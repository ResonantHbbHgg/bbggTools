// C++ headers
#include "flashgg/bbggTools/interface/bbgg2DFitter.h"

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

// namespaces
using namespace std;
using namespace RooFit;
using namespace RooStats;
namespace po = boost::program_options;

//Important options first
//Bool_t doblinding = false; //True if you want to blind

// this one is for 2D fit


int main(int argc, const char* argv[])
{
  Float_t mass = 1;
  Float_t lumi = 1;
  Bool_t doBands = 1;
  int version = 1;
  string analysisType = "a";
  string nonresFile = "b";
  Bool_t useSigTheoryUnc = 1;
  int sigMass = 1;
  Bool_t doblinding = 1;
  Int_t NCAT = 1;
  bool addHiggs =1;
 /* 
  try
    {
      po::options_description desc("Allowed options");
      desc.add_options()
	("help,h", "produce help message")
    ("lumi,l", po::value<float>(&lumi)->default_value(1900.0), "Integrated luminosity analyzed.")
    ("blind,b", po::value<bool>(&doblinding)->default_value(0), "Do blinded analysis.")
    ("addHiggs", po::value<bool>(&addHiggs)->default_value(true), "addHiggs?")
	("Hmass", po::value<float>(&mass)->default_value(125.02), "Mass of SM Higgs. Default is 125.02.")
	("doBands", po::value<bool>(&doBands)->default_value(true), "Option to calculate and show 1,2 sigma bands on bkg fit.")
	("version,v", po::value<int>(&version)->default_value(42), "Version for limit trees.")
	("ncat,n", po::value<int>(&NCAT)->default_value(2), "Number of categories to fit")
	("sigMass", po::value<int>(&sigMass)->default_value(0), "Mass of signal. 0 is for nonresonant.")
	("analysisType", po::value<string>(&analysisType)->default_value("fitTo2D_nonresSearch_withKinFit"), "Can choose among fitTo{2D,FTR14001}_{nonres,res}Search_with{RegKin,Kin}Fit")
	("nonresFile", po::value<string>(&nonresFile)->default_value("Lam_1d0_Yt_1d0_c2_0d0"), "nonres signal to run in the case sigMass is 0. default is the SM value.")
	("useSigTheoryUnc", po::value<bool>(&useSigTheoryUnc)->default_value(false), "option to add an uncertainty to the datacard for the SM diHiggs theory uncertainty. Default is off.")
        ;
      po::variables_map vm;
      po::store(po::parse_command_line(argc, argv, desc), vm);
      po::notify(vm);
      if (vm.count("help")) {
	cout << desc << "\n";
	return 1;
      }
    } catch(exception& e) {
    cerr << "error: " << e.what() << "\n";
    return 1;
  } catch(...) {
    cerr << "Exception of unknown type!\n";
  }
  // end of argument parsing
  */

  TString fileBaseName = TString::Format("hgg.mH%.1f_8TeV", mass);
  TString fileHiggsNameggh = TString::Format("hgg.hig.mH%.1f_8TeV.ggh", mass);
  TString fileHiggsNametth = TString::Format("hgg.hig.mH%.1f_8TeV.tth", mass);
  TString fileHiggsNamevbf = TString::Format("hgg.hig.mH%.1f_8TeV.vbf", mass);
  TString fileHiggsNamevh = TString::Format("hgg.hig.mH%.1f_8TeV.vh", mass);
  TString fileHiggsNamebbh = TString::Format("hgg.hig.mH%.1f_8TeV.bbh", mass);
  TString fileBkgName = "hgg.inputbkg_8TeV";
  TString card_name = "models_2D.rs"; // put the model parameters here!
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
  
  //Object
  bbgg2DFitter TheFitter = bbgg2DFitter( w, sigMass, lumi, doblinding, NCAT, addHiggs );
  TheFitter.style();
  
  RooFitResult* fitresults;

  // the limit trees to be addeed
  //
  TString dir = TString::Format("/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v%d/v%d_%s",version,version,analysisType.c_str());

  TString hhiggsggh = TString::Format("%s/ggh_m125_powheg_8TeV_m%d.root",dir.Data(),sigMass);
  TString hhiggstth = TString::Format("%s/tth_m125_8TeV_m%d.root",dir.Data(),sigMass);;
  TString hhiggsvbf = TString::Format("%s/vbf_m125_8TeV_m%d.root",dir.Data(),sigMass);;
  TString hhiggsvh =  TString::Format("%s/wzh_m125_8TeV_zh_m%d.root",dir.Data(),sigMass);;
  TString hhiggsbbh = TString::Format("%s/bbh_m125_8TeV_m%d.root",dir.Data(),sigMass);;
  //
  TString ddata = TString::Format("%s/Data_m%d.root",dir.Data(),sigMass);
  TString ssignal;
  if (sigMass == 260) ssignal = TString::Format("%s/MSSM_m260_8TeV_m260.root",dir.Data());
  else if (sigMass >= 270) ssignal = TString::Format("%s/Radion_m%d_8TeV_m%d.root",dir.Data(),sigMass,sigMass);
  else ssignal = TString::Format("%s/ggHH_%s_8TeV_m0.root",dir.Data(),nonresFile.c_str());

  cout<<"Signal: "<<ssignal<<endl;
  cout<<"Data: "<<ddata<<endl;
  //

  TheFitter.AddSigData( mass,ssignal);
  cout<<"SIGNAL ADDED"<<endl;
  TheFitter.SigModelFit( mass); // constructing signal pdf
  TheFitter.MakeSigWS( fileBaseName);
  TheFitter.MakePlots( mass);
  cout<<" did signal WS's"<<endl;

  //
  cout<<"Higgs: "<<hhiggsggh<<endl;
  TheFitter.AddHigData( mass,hhiggsggh,0);
  TheFitter.HigModelFit( mass,0); // constructing higgs pdf
  TheFitter.MakeHigWS( fileHiggsNameggh,0);
  //
  cout<<"Higgs: "<<hhiggstth<<endl;
  TheFitter.AddHigData( mass,hhiggstth,1);
  TheFitter.HigModelFit( mass,1); // constructing higgs pdf
  TheFitter.MakeHigWS( fileHiggsNametth,1);
  //
  cout<<"Higgs: "<<hhiggsvbf<<endl;
  TheFitter.AddHigData( mass,hhiggsvbf,2);
  TheFitter.HigModelFit( mass,2); // constructing higgs pdf
  TheFitter.MakeHigWS( fileHiggsNamevbf,2);
  //
  cout<<"=============================== Higgs: ============================= "<<hhiggsvh<<endl;
  TheFitter.AddHigData( mass,hhiggsvh,3);
  TheFitter.HigModelFit( mass,3); // constructing higgs pdf
  TheFitter.MakeHigWS( fileHiggsNamevh,3);
  cout<<"HIGGS ADDED"<<endl;
  //
  cout<<"=============================== Higgs: ============================= "<<hhiggsbbh<<endl;
  TheFitter.AddHigData( mass,hhiggsbbh,4);
  TheFitter.HigModelFit( mass,4); // constructing higgs pdf
  TheFitter.MakeHigWS( fileHiggsNamebbh,4);
  cout<<"HIGGS ADDED"<<endl;
  TheFitter.MakePlotsHiggs( mass);
  //


  TheFitter.AddBkgData(ddata);
  //w->Print("v");
  TheFitter.PrintWorkspace();
  cout<<"BKG ADDED"<<endl;
  bool dobands=true;
  fitresults = TheFitter.BkgModelFit( dobands); // this is berestein 3
  TheFitter.MakeBkgWS( fileBkgName);
  // construct the models to fit
  //
  TheFitter.MakeDataCardonecatnohiggs( fileBaseName, fileBkgName, useSigTheoryUnc);
  TheFitter.MakeDataCard( fileBaseName, fileBkgName, fileHiggsNameggh, fileHiggsNametth, fileHiggsNamevbf, fileHiggsNamevh, fileHiggsNamebbh, useSigTheoryUnc);

  cout<< "here"<<endl;
  fitresults->Print();

  return 0;
} // close runfits
