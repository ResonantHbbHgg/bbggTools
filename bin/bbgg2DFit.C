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

//Important options first
//Bool_t doblinding = false; //True if you want to blind

// this one is for 2D fit
/*       Example of xml config file
         <signal>
                <type>          Radion          </type>
                <dir>           /a/b/c/d        </dir>
                <mass>          270             </mass>
                <nonresFile>    Lam_1d0_Yt_1d0_c2_0d0 </nonresFile> <!-- nonres signal to run in the case sigMass is 0. default is the SM value. -->
        </signal>
        <other>
                <integratedLumi> 0.15   </integratedLumi>
                <energy>        13TeV   </energy>
                <higgsMass>     125.0   </higgsMass>
                <addHiggs>      1       </addHiggs>
                <doBlinding>    0       </doBlinding> <!-- do blinded analysis? -->
                <doBands>       1       </doBands> <!-- Option to calculate and show 1,2 sigma bands on bkg fit -->
                <ncat>          2       </ncat> <!-- Number of categories to fit -->
                <analysisType>  fitTo2D_nonresSearch_withKinFit </analysisType> <!-- Can choose among fitTo{2D,FTR14001}_{nonres,res}Search_with{RegKin,Kin}Fit -->
                <useSigTheoryUnc> 0     </useSigTheoryUnc> <!-- option to add an uncertainty to the datacard for the SM diHiggs theory uncertainty. Default is off -->
        </other>
*/
 
int main(int argc, const char* argv[])
{
  if(argc < 2) {
	cout << "Please, parse xml configuration file!" << endl;
	return 0;
  }

//Parameters to set
  Float_t mass = 1;
  Float_t lumi = 1;
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
  string energy = "13TeV";
  string cardName = "";


  //Read config file
//  using boost::property_tree::ptree;
 // ptree configs;
  boost::property_tree::ptree pt;
  boost::property_tree::read_json( argv[1], pt );

  cout << "Reading input configuration file..." << endl;
//  read_xml(argv[1], configs);
 // const ptree & formats = configs.get_child("configs", empty_ptree());
//  BOOST_FOREACH( boost::property_tree::ptree::value_type const& v, formats )
  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "" ) )
  {

        cout << "Reading " << rowPair.first << " options..." << endl;
        if (rowPair.first == "signal") {
                signalType = rowPair.second.get<std::string>("type");
                signalDir = rowPair.second.get<std::string>("dir");
                sigMass = rowPair.second.get<int>("mass");
		cardName = rowPair.second.get<std::string>("signalModelCard");
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
                lumi = rowPair.second.get<float>("integratedLumi");
                energy = rowPair.second.get<std::string>("energy");
                mass = rowPair.second.get<float>("higgsMass");
                addHiggs = rowPair.second.get<bool>("addHiggs");
                doblinding = rowPair.second.get<bool>("doBlinding");
                doBands = rowPair.second.get<bool>("doBands");
                NCAT = rowPair.second.get<int>("ncat");
                useSigTheoryUnc = rowPair.second.get<bool>("useSigTheoryUnc");
                analysisType = rowPair.second.get<string>("analysisType");
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

/*
	std::cout << rowPair.first << ": " << std::endl;
	cout << "brababa " << rowPair.second.get<string>("type") << endl;
	BOOST_FOREACH( boost::property_tree::ptree::value_type const& itemPair, rowPair.second ) 
        {
            std::cout << "\t" << itemPair.first << " " << itemPair.second.data() << endl;
 
            BOOST_FOREACH( boost::property_tree::ptree::value_type const& node, itemPair.second ) 
            {
                std::cout << node.first << " ";
            }
 
            std::cout << std::endl;
 
        }
 
        std::cout << std::endl;


	const ptree & attributes = v.second.get_child("<xmlattr>", empty_ptree());
	cout << "v first " <<  v.first << endl;
	BOOST_FOREACH( boost::property_tree::ptree::value_type const& f, attributes ) {
		cout << f.first.data() << " " << f.second.data() << endl;
	}
*/	
/*
	cout << "Reading " << v.first << " options..." << endl;
	if (v.first == "signal") {
		signalType = v.second.get<std::string>("type");
		signalDir = v.second.get<std::string>("dir");
		sigMass = v.second.get<int>("mass");
		cout << "\t Signal type: " << signalType << endl;
		cout << "\t Signal samples location: " << signalDir << endl;
		cout << "\t Signal mass: " << mass << endl;
		if (mass == 0 ) {
			cout << "Mass == 0 means non-resonant analysis, therefore:" << endl;
			nonresFile = v.second.get<std::string>("nonresFile");
			cout << "\t Non resonant file: " << nonresFile << endl;
		}
	}
	if (v.first == "other") {
		lumi = v.second.get<float>("integratedLumi");
		energy = v.second.get<std::string>("energy");
		mass = v.second.get<float>("higgsMass");
		addHiggs = v.second.get<bool>("addHiggs");
		doblinding = v.second.get<bool>("doBlinding");
		doBands = v.second.get<bool>("doBands");
		NCAT = v.second.get<int>("ncat");
		useSigTheoryUnc = v.second.get<bool>("useSigTheoryUnc");
		analysisType = v.second.get<string>("analysisType");
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
	}*/
  }
//  return 0;

/*
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
 */
/*
  try
    {
      boost::program_options::options_description desc("Allowed options");
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
      boost::program_options::variables_map vm;
      boost::program_options::store(po::parse_command_line(argc, argv, desc), vm);
      boost::program_options::notify(vm);
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
//  TString card_name = "LimitModels/models_2D.rs"; // put the model parameters here!
  TString card_name(cardName); // put the model parameters here!
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
