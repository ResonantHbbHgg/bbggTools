// C++ headers
#include "flashgg/bbggTools/interface/bbgg2DFitter.h"

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <boost/program_options.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
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
 
int main(int argc, const char* argv[])
{
  if(argc < 2) {
	cout << "Please, parse xml configuration file!" << endl;
	return 0;
  }

//Parameters to set
  Float_t mass = 1;
  Float_t lumi = 1;
  Float_t minMgg = 100.;
  Float_t maxMgg = 180.;
  Float_t minMjj = 60.;
  Float_t maxMjj = 180;
  Float_t minMggHig = 115.;
  Float_t maxMggHig = 135.;
  Float_t minMjjHig = 60.;
  Float_t maxMjjHig = 180.;
  Bool_t doBands = 1;
  int version = 44;
  string analysisType = "a";
  string nonresFile = "b";
  Bool_t useSigTheoryUnc = 1;
  int sigMass = 1;
  Bool_t doblinding = 1;
  Int_t NCAT = 1;
  bool addHiggs =1;
  string signalType = "a";
  string signalDir = "a";
  string energy = "13";
  string cardName = "";
  string sigFile = "";
  string dataFile = "";

  //Read config file
//  using boost::property_tree::ptree;
 // ptree configs;
  boost::property_tree::ptree pt;
  boost::property_tree::read_json( argv[1], pt );

  cout << "Reading input configuration file..." << endl;
  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "" ) )
  {

        cout << "Reading " << rowPair.first << " options..." << endl;
        if (rowPair.first == "signal") {
                signalType = rowPair.second.get<std::string>("type");
                signalDir = rowPair.second.get<std::string>("dir");
                sigMass = rowPair.second.get<int>("mass");
		cardName = rowPair.second.get<std::string>("signalModelCard");
		sigFile = rowPair.second.get<std::string>("signalFile");
                cout << "\t Signal type: " << signalType << endl;
                cout << "\t Signal samples location: " << signalDir << endl;
		cout << "\t Signal model card: " << cardName << endl;
                cout << "\t Signal mass: " << mass << endl;
                if (mass == 0 ) {
                        cout << "Mass == 0 means non-resonant analysis, therefore:" << endl;
                        nonresFile = rowPair.second.get<std::string>("nonresFile");
                        cout << "\t Non resonant file: " << nonresFile << endl;
                }
        }
        if (rowPair.first == "other") {
		dataFile = rowPair.second.get<std::string>("dataFile");
                lumi = rowPair.second.get<float>("integratedLumi");
                energy = rowPair.second.get<std::string>("energy");
                mass = rowPair.second.get<float>("higgsMass");
                addHiggs = rowPair.second.get<bool>("addHiggs");
                doblinding = rowPair.second.get<bool>("doBlinding");
                doBands = rowPair.second.get<bool>("doBands");
                NCAT = rowPair.second.get<int>("ncat");
                useSigTheoryUnc = rowPair.second.get<bool>("useSigTheoryUnc");
                analysisType = rowPair.second.get<string>("analysisType");
		minMgg = rowPair.second.get<float>("minMgg");
		maxMgg = rowPair.second.get<float>("maxMgg");
		minMjj = rowPair.second.get<float>("minMjj");
		maxMjj = rowPair.second.get<float>("maxMjj");
		minMggHig = rowPair.second.get<float>("minMggHig");
		maxMggHig = rowPair.second.get<float>("maxMggHig");
		minMjjHig = rowPair.second.get<float>("minMjjHig");
		maxMjjHig = rowPair.second.get<float>("maxMjjHig");
                cout << "Running options: " << endl;
                cout << "\t Integrated luminosity: " << lumi << endl;
                cout << "\t Center of mass energy: " << energy << endl;
                cout << "\t Higgs mass: " << mass << endl;
                cout << "\t Add Higgs: " << addHiggs << endl;
                cout << "\t Do blinded analysis: " << doblinding << endl;
                cout << "\t Calculate and show 1 and 2 sigma bands on background fit: " << doBands << endl;
                cout << "\t Number of categories to fit: " << NCAT << endl;
                cout << "\t Analysis type: " << analysisType << endl;
                cout << "\t Add an uncertainty to the datacard for the SM diHiggs theory uncertainty: " << useSigTheoryUnc << endl;
        }
}


  TString fileBaseName = TString::Format("hgg.mH%.1f_8TeV", mass);
  TString fileHiggsNameggh = TString::Format("hgg.hig.mH%.1f_8TeV.ggh", mass);
  TString fileHiggsNametth = TString::Format("hgg.hig.mH%.1f_8TeV.tth", mass);
  TString fileHiggsNamevbf = TString::Format("hgg.hig.mH%.1f_8TeV.vbf", mass);
  TString fileHiggsNamevh = TString::Format("hgg.hig.mH%.1f_8TeV.vh", mass);
  TString fileHiggsNamebbh = TString::Format("hgg.hig.mH%.1f_8TeV.bbh", mass);
  TString fileBkgName = "hgg.inputbkg_8TeV";

  TString card_name(cardName); // put the model parameters here!
  bbgg2DFitter TheFitter = bbgg2DFitter( sigMass, lumi, doblinding, NCAT, addHiggs );
  TheFitter.SetMinMaxMggMjj(minMgg, maxMgg, minMjj, maxMjj);
  TheFitter.SetMinMaxMggMjjHig(minMggHig, maxMggHig, minMjjHig, maxMjjHig);
  std::string CardName = TheFitter.MakeModelCard(card_name.Data());
  
  if(CardName == "")
	return 0;

  HLFactory hlf("HLFactory", CardName, false);
  RooWorkspace* w = hlf.GetWs();
  
  //Object
  TheFitter.SetWorkspace(w);
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

  ssignal = sigFile;
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


//  ddata = "/tmp/rateixei/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/S15V7_v0/LimitTrees/600_900/LT_data.root";
  ddata = dataFile;
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
