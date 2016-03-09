//////////////////////////////////////////////////////////////////////////
//
// 'ORGANIZATION AND SIMULTANEOUS FITS' RooFit tutorial macro #505
// 
// Reading and writing ASCII configuration files
//
//
//
// 07/2008 - Wouter Verkerke 
// 
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "TString.h"
#include "TChain.h"
#include "TLegend.h"
#include "RooMinuit.h"
#include "TH1F.h"

//Boost
#include <boost/program_options.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

using namespace RooFit;
namespace po = boost::program_options;


int main(int argc, const char* argv[])
{
  // C r e a t e  w o r k s p a c e 
  // -------------------------------
 
  if(argc < 2) {
    std::cout << "Please, parse JSON configuration file!" << std::endl;
    return 0;
  }

  //Parameters
  std::string rootName;
  std::string doBlinding;
  std::string blindingCut;
  std::string prefix;
  std::string plotsDir;
  std::string error;
  std::vector<std::string> vars;
  std::vector<std::string> varNames;
  std::vector<std::string> nbins;
  std::vector<std::string> functions;
  std::vector<std::string> functionsToFit;
  std::vector<std::string> categories;
  std::vector<std::string> categoriesCuts;
  std::vector<TH1F> funcHistograms;

/*
enum EColor { kWhite =0,   kBlack =1,   kGray=920,
              kRed   =632, kGreen =416, kBlue=600, kYellow=400, kMagenta=616, kCyan=432,
              kOrange=800, kSpring=820, kTeal=840, kAzure =860, kViolet =880, kPink=900 };
*/

  int colors[] = {632, 416, 600, 400, 616, 432, 800, 820, 840, 860};

  //Read json
  boost::property_tree::ptree pt;
  boost::property_tree::read_json( argv[1], pt );
  rootName = pt.get_child("file").data();
  prefix = pt.get_child("prefix").data();
  plotsDir = pt.get_child("plotsDir").data();
  error = pt.get_child("error").data();

  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "vars" ) ){
    std::cout << rowPair.first << "\t" << rowPair.second.data() << std::endl;
    vars.push_back(rowPair.first);
    varNames.push_back(rowPair.second.data());
  }

  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "nbins" ) ){
    nbins.push_back(rowPair.second.data());
  }

  if(nbins.size() != vars.size()){
    std::cout << "Size of nbins vector must be the same as variables vector..." << std::endl;
    return 0;
  }

  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "functions" ) ){
    functions.push_back(rowPair.second.data());
  }

  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "functionsToFit" ) ){
    functionsToFit.push_back(rowPair.second.data());
  }
  
  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "categories" ) ){
    categories.push_back( rowPair.first );
    categoriesCuts.push_back(rowPair.second.data());
  }

  BOOST_FOREACH( boost::property_tree::ptree::value_type const& rowPair, pt.get_child( "blinding" ) ){
    blindingCut = rowPair.second.data();
    doBlinding = rowPair.first;
    break;
  }
  
  for( unsigned int i = 0; i < vars.size(); i++){
    for ( unsigned int cat = 0; cat < categories.size(); cat++) {
      RooWorkspace* w = new RooWorkspace("w");

      // A d d  V a r i a b l e
      // ----------------------

      // Initialize fit variable
      w->factory( vars[i].c_str() );
      std::string sVar = ( (TObjString *) ( (TObjArray *) (TString(vars[i]).Tokenize("[")) )->At(0) )->String().Data();
      std::cout << "Variable to be fit: " << sVar << std::endl;
      w->var( sVar.c_str() )->SetTitle(varNames[i].c_str());
      w->var( sVar.c_str() )->setUnit("GeV");

      //Make plot frame
      RooPlot* frame = w->var(sVar.c_str())->frame(Bins(TString( nbins[i] ).Atoi() ));

      //Initialize all needed PDFs
      for ( unsigned int f = 0; f < functions.size(); f++){
        TString thisFunction(functions[f]);
        thisFunction.ReplaceAll( "$VAR", sVar);
        std::cout << "Function added: " << thisFunction << std::endl;
        w->factory( thisFunction.Data() );
      }

      //Get distribution from file
      TChain* treeTemp = new TChain("TCVARS");
      treeTemp->AddFile(TString(rootName));
      TTree* tree = (TTree*) treeTemp->CopyTree( categoriesCuts[cat].c_str() );
      std::cout << "Reading file: " << rootName << "\n\t with " << tree->GetEntries() << " entries \n\t Making the following cut: " << categoriesCuts[cat] << std::endl;

      RooDataSet* data = new RooDataSet( "data", "data", tree, *(w->var(sVar.c_str())) );
      w->import(*data);

      //Plot data on frame
      data->plotOn(frame);
      //Data histogram
      TH1F* dataHist = (TH1F*) data->createHistogram("dataHist", *(w->var(sVar.c_str())) );

      //TLegend
      TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);

      //Start fitting businness - fit for each fit function
      std::vector<double> kolmos;
      for ( unsigned int ff = 0; ff < functionsToFit.size(); ff++) {
	const char* modelName = functionsToFit[ff].c_str();
	std::cout << "Fitting model: " << modelName << std::endl;
	if(w->pdf(modelName) == 0 ){
	  std::cout << "Model " << modelName << " not found in workspace! Are you sure you have added it to the list of funtions?" << std::endl;
	  continue;
	}
	RooAbsReal* nll = w->pdf( modelName )->createNLL(*data, NumCPU(8));
	RooMinuit minuit(*nll);
	minuit.migrad();
	minuit.hesse();
	RooFitResult * fitResult = minuit.save();
	TH1F* histo = (TH1F*) w->pdf( modelName )->createHistogram( TString::Format("%s_Histo", modelName), *(w->var(sVar.c_str())) );
	histo->Scale(dataHist->Integral());
	float kolmo = histo->KolmogorovTest(dataHist);
	std::cout << "############## kolmo: " <<  histo->KolmogorovTest(dataHist) << std::endl;
	kolmos.push_back(kolmo);

	w->pdf( modelName )->plotOn(frame, LineColor( colors[ff] ), Name(modelName), MoveToBack());
	if(error == "1")
	  w->pdf( modelName )->plotOn(frame, FillColor( kGray ), VisualizeError(*fitResult, 2), MoveToBack());
	
	leg->AddEntry( frame->findObject(modelName), TString::Format("%s KS Test = %f", modelName, kolmo), "l" );

	delete nll;
	delete histo;
      }

      w->Print();

      //Save plots
      TCanvas* c = new TCanvas("a", "a", 1200, 1000);
      frame->Draw();
      leg->Draw("same");
      TString plotName = TString::Format("%s/%s_%s_%s.png", plotsDir.c_str(), prefix.c_str(), categories[cat].c_str(), sVar.c_str());
      c->SaveAs(plotName);
      plotName = TString::Format("%s/%s_%s_%s.pdf", plotsDir.c_str(), prefix.c_str(), categories[cat].c_str(), sVar.c_str());
      c->SaveAs(plotName);


      delete w;
      delete tree;
      delete data;
      delete frame;
      delete treeTemp;
      delete dataHist;
    }
  }


  return 0;
}
