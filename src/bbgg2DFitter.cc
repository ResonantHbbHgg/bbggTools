#define bbgg2DFitter_cxx
#include "flashgg/bbggTools/interface/bbgg2DFitter.h"

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
#include <RooAbsMoment.h>
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

RooArgSet* bbgg2DFitter::defineVariables()
{
  RooRealVar* mgg = new RooRealVar("mgg","M(#gamma#gamma)",100,180,"GeV");
  RooRealVar* mtot = new RooRealVar("mtot","M(#gamma#gammajj)",200,1600,"GeV");
  RooRealVar* mjj = new RooRealVar("mjj","M(jj)",60,180,"GeV");
  RooRealVar* evWeight = new RooRealVar("evWeight","HqT x PUwei",0,100,"");
  RooCategory* cut_based_ct = new RooCategory("cut_based_ct","event category 4") ;
  //
  cut_based_ct->defineType("cat4_0",0);
  cut_based_ct->defineType("cat4_1",1);
  cut_based_ct->defineType("cat4_2",2);
  cut_based_ct->defineType("cat4_3",3);
  //
  RooArgSet* ntplVars = new RooArgSet(*mgg, *mjj, *cut_based_ct, *evWeight);
  ntplVars->add(*mgg);
  ntplVars->add(*mtot);
  ntplVars->add(*mjj);
  ntplVars->add(*cut_based_ct);
  return ntplVars;
}

void bbgg2DFitter::AddSigData(float mass, TString signalfile)
{
    cout << "================= Add Signal==============================" << endl;
    //Luminosity
    RooRealVar lumi("lumi","lumi", _lumi);
    _w->import(lumi);
    //Define variables
    RooArgSet* ntplVars = bbgg2DFitter::defineVariables();
    //Signal file & tree
    TFile sigFile(signalfile);
    TTree* sigTree = (TTree*) sigFile.Get("TCVARS");
    //Data set
    RooDataSet sigScaled(
    		       "sigScaled",
    		       "dataset",
    		       sigTree,
    		       *ntplVars,
    		       _cut,
    		       "evWeight");
    cout << "======================================================================" << endl;
    RooDataSet* sigToFit[_NCAT];
    TString cut0 = " && 1>0";
    for ( int i=0; i<_NCAT; ++i){
        sigToFit[i] = (RooDataSet*) sigScaled.reduce(
    						 RooArgList(*_w->var("mgg"),*_w->var("mjj")),
    						 _cut+TString::Format(" && cut_based_ct==%d ",i)+cut0); //This defines each category
        _w->import(*sigToFit[i],Rename(TString::Format("Sig_cat%d",i)));
    }
  // Create full signal data set without categorization
    RooDataSet* sigToFitAll = (RooDataSet*) sigScaled.reduce(
  							   RooArgList(*_w->var("mgg"),*_w->var("mjj")),
  							   _cut);
    cout << "======================================================================" << endl;
    _w->import(*sigToFitAll,Rename("Sig"));
  // here we print the number of entries on the different categories
    cout << "========= the number of entries on the different categories ==========" << endl;
    cout << "---- one channel: " << sigScaled.sumEntries() << endl;
    for (int c = 0; c < _NCAT; ++c) {
      Float_t nExpEvt = sigToFit[c]->sumEntries();
      cout << TString::Format("nEvt exp. cat%d : ",c) << nExpEvt
  	 << TString::Format(" eff x Acc cat%d : ",c)
  	 << "%"
  	 << endl;
    } // close ncat
    cout << "======================================================================" << endl;
    sigScaled.Print("v");
    return;
}

void bbgg2DFitter::AddHigData(float mass, TString signalfile, int higgschannel)
{
  RooArgSet* ntplVars = defineVariables();
  TFile higFile(signalfile);
  TTree* higTree = (TTree*) higFile.Get("TCVARS");
  RooDataSet higScaled(
		       "higScaled1",
		       "dataset",
		       higTree, // all variables of RooArgList
		       *ntplVars,
		       _cut,
		       "evWeight");
  //
  RooDataSet* higToFit[_NCAT];
  TString cut0 = "&& 1>0";
  //
  // we take only mtot to fit to the workspace, we include the cuts
  for ( int i=0; i<_NCAT; ++i){
    higToFit[i] = (RooDataSet*) higScaled.reduce(
						 RooArgList(*_w->var("mgg"),*_w->var("mjj")),
						 _cut+TString::Format(" && cut_based_ct==%d ",i)+cut0);
    _w->import(*higToFit[i],Rename(TString::Format("Hig_%d_cat%d",higgschannel,i)));
  }
  // Create full signal data set without categorization
  RooDataSet* higToFitAll = (RooDataSet*) higScaled.reduce(
							    RooArgList(*_w->var("mgg"),*_w->var("mjj")),
							    _cut);
  _w->import(*higToFitAll,Rename("Hig"));
  // here we print the number of entries on the different categories
  cout << "========= the number of entries on the different categories (Higgs data) ==========" << endl;
  cout << "---- one channel: " << higScaled.sumEntries() << endl;
  for (int c = 0; c < _NCAT; ++c) {
    Float_t nExpEvt = higToFit[c]->sumEntries();
    cout << TString::Format("nEvt exp. cat%d : ",c) << nExpEvt
	 << TString::Format(" eff x Acc cat%d : ",c)
	 << "%"
	 << endl;
  }
  cout << "======================================================================" << endl;
  higScaled.Print("v");
  return;
}

void bbgg2DFitter::AddBkgData(TString datafile)
{
    //Define variables
    RooArgSet* ntplVars = bbgg2DFitter::defineVariables();
    RooRealVar weightVar("weightVar","",1,0,1000);
    weightVar.setVal(1.);
    TFile dataFile(datafile);
    TTree* dataTree = (TTree*) dataFile.Get("TCVARS");
    RooDataSet Data("Data","dataset",dataTree,*ntplVars,"","evWeight");
    RooDataSet* dataToFit[_NCAT];
    RooDataSet* dataToPlot[_NCAT];
    TString cut0 = "&& 1>0";
    TString cut1 = "&& 1>0";
    
    cout<<"================= Add Bkg ==============================="<<endl;

      for( int i=0; i<_NCAT; ++i){
        dataToFit[i] = (RooDataSet*) Data.reduce(
    					     RooArgList(*_w->var("mgg"),*_w->var("mjj")),
    					     _cut+TString::Format(" && cut_based_ct==%d",i));
        if(_doblinding){
          dataToPlot[i] = (RooDataSet*) Data.reduce(
    						RooArgList(*_w->var("mgg"),*_w->var("mjj")),
    						_cut+TString::Format(" && cut_based_ct==%d",i)
    						+cut0);
        }
        else{
          dataToPlot[i] = (RooDataSet*) Data.reduce(
    						RooArgList(*_w->var("mgg"),*_w->var("mjj")),
    						_cut+TString::Format(" && cut_based_ct==%d",i) );
        }
      }
      for (int c = 0; c < _NCAT; ++c) {
        _w->import(*dataToFit[c],Rename(TString::Format("Data_cat%d",c)));
        _w->import(*dataToPlot[c],Rename(TString::Format("Dataplot_cat%d",c)));
      }
      // Create full data set without categorization
      RooDataSet* data = (RooDataSet*) Data.reduce(RooArgList(*_w->var("mgg"),*_w->var("mjj")),_cut);
      _w->import(*data, Rename("Data"));
      data->Print("v");
      return;
}

void bbgg2DFitter::SigModelFit(float mass)
{
  //******************************************//
  // Fit signal with model pdfs
  //******************************************//
  // four categories to fit
  RooDataSet* sigToFit[_NCAT];
  RooAbsPdf* mggSig[_NCAT];
  RooAbsPdf* mjjSig[_NCAT];
  RooProdPdf* SigPdf[_NCAT];
  // fit range
  Float_t minSigFitMgg(115),maxSigFitMgg(135); //This should be an option
  Float_t minSigFitMjj(60),maxSigFitMjj(180); //This should be an option
  RooRealVar* mgg = _w->var("mgg");
  RooRealVar* mjj = _w->var("mjj");
  mgg->setRange("SigFitRange",minSigFitMgg,maxSigFitMgg);
  mjj->setRange("SigFitRange",minSigFitMjj,maxSigFitMjj);

  for (int c = 0; c < _NCAT; ++c) {
    // import sig and data from workspace
    sigToFit[c] = (RooDataSet*) _w->data(TString::Format("Sig_cat%d",c));
    mggSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggSig_cat%d",c));
    mjjSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjSig_cat%d",c));
    SigPdf[c] = new RooProdPdf(TString::Format("SigPdf_cat%d",c),"",RooArgSet(*mggSig[c], *mjjSig[c]));

    ((RooRealVar*) _w->var(TString::Format("mgg_sig_m0_cat%d",c)))->setVal(mass);
    //RooRealVar* peak = w->var(TString::Format("mgg_sig_m0_cat%d",c));
    //peak->setVal(MASS);
    cout << "OK up to now... Mass point: " <<mass<< endl;

    SigPdf[c]->fitTo(*sigToFit[c],Range("SigFitRange"),SumW2Error(kTRUE));
    cout << "old = " << ((RooRealVar*) _w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal() << endl;

    double mPeak = ((RooRealVar*) _w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal()+(mass-125.0); // shift the peak //This should be an option
    ((RooRealVar*) _w->var(TString::Format("mgg_sig_m0_cat%d",c)))->setVal(mPeak); // shift the peak

    cout << "mPeak = " << mPeak << endl;
    cout << "new mPeak position = " << ((RooRealVar*) _w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal() << endl;

    // IMPORTANT: fix all pdf parameters to constant, why?
    RooArgSet sigParams( *_w->var(TString::Format("mgg_sig_m0_cat%d",c)),
			 *_w->var(TString::Format("mgg_sig_sigma_cat%d",c)), // CB sigma
			 *_w->var(TString::Format("mgg_sig_alpha_cat%d",c)),
			 *_w->var(TString::Format("mgg_sig_n_cat%d",c)),
			 *_w->var(TString::Format("mgg_sig_gsigma_cat%d",c)), //gaussian sigma
			 *_w->var(TString::Format("mgg_sig_frac_cat%d",c)));
    sigParams.add(RooArgSet(
			       *_w->var(TString::Format("mjj_sig_m0_cat%d",c)),
		           *_w->var(TString::Format("mjj_sig_sigma_cat%d",c)),
		           *_w->var(TString::Format("mjj_sig_alpha_cat%d",c)),
		           *_w->var(TString::Format("mjj_sig_n_cat%d",c)),
		           *_w->var(TString::Format("mjj_sig_gsigma_cat%d",c)),
		           *_w->var(TString::Format("mjj_sig_frac_cat%d",c))) );

    _w->defineSet(TString::Format("SigPdfParam_cat%d",c), sigParams);
    SetConstantParams(_w->set(TString::Format("SigPdfParam_cat%d",c)));

    _w->import(*SigPdf[c]);

  } // close for ncat
} // close signal model fit

void bbgg2DFitter::HigModelFit(float mass, int higgschannel)
{
  // four categories to fit
//  RooBernstein* mjjHig_Ber0;
//  RooBernstein* mjjHig_Ber1;
//  RooBernstein* mjjHig_Ber2;
  RooPolynomial* mjjHig_pol0;
  RooDataSet* higToFit[_NCAT];
  RooAbsPdf* mggHig[_NCAT];
  RooAbsPdf* mjjHig[_NCAT];
  RooProdPdf* HigPdf[_NCAT];
  // fit range
  Float_t minHigMggFit(115),maxHigMggFit(135);//This should be an option
  Float_t minHigMjjFit(60),maxHigMjjFit(180);//This should be an option
  RooRealVar* mgg = _w->var("mgg");
  RooRealVar* mjj = _w->var("mjj");
  mgg->setRange("HigFitRange",minHigMggFit,maxHigMggFit);
  mjj->setRange("HigFitRange",minHigMjjFit,maxHigMjjFit);

  for (int c = 0; c < _NCAT; ++c) {
    // import sig and data from workspace
    mjjHig_pol0 = new RooPolynomial(TString::Format("mjjHig_pol0_%d_cat%d",higgschannel,c),"",*mjj,RooArgList());
    higToFit[c] = (RooDataSet*) _w->data(TString::Format("Hig_%d_cat%d",higgschannel,c));
    mggHig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggHig_%d_cat%d",higgschannel,c));
    if(higgschannel == 1 || higgschannel == 3) mjjHig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_%d_cat%d",higgschannel,c));
    else{
       mjjHig[c] = mjjHig_pol0;
    }
    HigPdf[c] = new RooProdPdf(TString::Format("HigPdf_%d_cat%d",higgschannel,c),"",RooArgSet(*mggHig[c], *mjjHig[c]));
    //((RooRealVar*) w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->setVal(MASS);
    cout << "OK up to now... Mass point: " <<mass<< endl;
    HigPdf[c]->fitTo(*higToFit[c],Range("HigFitRange"),SumW2Error(kTRUE));
    RooArgSet* paramsMjj;
    paramsMjj = (RooArgSet*) mjjHig[c]->getParameters(*mjj);
    TIterator* iterMjj = paramsMjj->createIterator();
    TObject* tempObjMjj=0;
    while((tempObjMjj=iterMjj->Next())){
                RooRealVar* var = (RooRealVar*)tempObjMjj;
                std::cout << "Variables after fit = " << tempObjMjj->GetName() << " " << var->getVal() << "+/-" << var->getError() << std::endl;
    }

    cout << "old = " << ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getVal() << endl;

    //There are very few events in some fits, so adjust the max by a good amount so the MASS-125.0 shift doesn't touch it.
    ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->setMax( ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getMax()+(mass-125.0) );
    double mPeak = ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getVal()+(mass-125.0); // shift the peak
    ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->setVal(mPeak); // shift the peak

    cout << "mPeak = " << mPeak << endl;
    cout << "new mPeak position = " << ((RooRealVar*) _w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)))->getVal() << endl;

    // IMPORTANT: fix all pdf parameters to constant
    RooArgSet sigParams( *_w->var(TString::Format("mgg_hig_m0_%d_cat%d",higgschannel,c)),
			 *_w->var(TString::Format("mgg_hig_sigma_%d_cat%d",higgschannel,c)),
			 *_w->var(TString::Format("mgg_hig_alpha_%d_cat%d",higgschannel,c)),
			 *_w->var(TString::Format("mgg_hig_n_%d_cat%d",higgschannel,c)),
			 *_w->var(TString::Format("mgg_hig_gsigma_%d_cat%d",higgschannel,c)),
			 *_w->var(TString::Format("mgg_hig_frac_%d_cat%d",higgschannel,c)) );
    if(higgschannel == 1 || higgschannel == 3){
       sigParams.add(RooArgSet(
			   *_w->var(TString::Format("mjj_hig_m0_%d_cat%d",higgschannel,c)),
			   *_w->var(TString::Format("mjj_hig_sigma_%d_cat%d",higgschannel,c)),
			   *_w->var(TString::Format("mjj_hig_alpha_%d_cat%d",higgschannel,c)),
			   *_w->var(TString::Format("mjj_hig_n_%d_cat%d",higgschannel,c)),
			   *_w->var(TString::Format("mjj_hig_gsigma_%d_cat%d",higgschannel,c)),
			   *_w->var(TString::Format("mjj_hig_frac_%d_cat%d",higgschannel,c)) ) );
    }
    _w->defineSet(TString::Format("HigPdfParam_%d_cat%d",higgschannel,c), sigParams);

    SetConstantParams(_w->set(TString::Format("HigPdfParam_%d_cat%d",higgschannel,c)));

    _w->import(*HigPdf[c]);
  } // close for ncat
} // close higgs model fit

void bbgg2DFitter::MakePlots(float mass)
{
  std::vector<TString> catdesc;
  if ( _NCAT == 2 ){
    catdesc.push_back(" High Purity");
    catdesc.push_back(" Med. Purity");
  }
  else{
    catdesc.push_back(" #splitline{High Purity}{High m_{#gamma#gammajj}^{kin}}");
    catdesc.push_back(" #splitline{Med. Purity}{High m_{#gamma#gammajj}^{kin}}");
    catdesc.push_back(" #splitline{High Purity}{Low m_{#gamma#gammajj}^{kin}}");
    catdesc.push_back(" #splitline{Med. Purity}{Low m_{#gamma#gammajj}^{kin}}");
  }
  // retrieve data sets from the workspace
  // RooDataSet* dataAll = (RooDataSet*) w->data("Data");
  RooDataSet* signalAll = (RooDataSet*) _w->data("Sig");
  //RooDataSet* higgsAll = (RooDataSet*) w->data("Hig");
  // blinded dataset
  // RooDataSet* data[ncat];
  RooDataSet* sigToFit[_NCAT];
  RooAbsPdf* mggGaussSig[_NCAT];
  RooAbsPdf* mggCBSig[_NCAT];
  RooAbsPdf* mggSig[_NCAT];
  RooAbsPdf* mjjGaussSig[_NCAT];
  RooAbsPdf* mjjCBSig[_NCAT];
  RooAbsPdf* mjjSig[_NCAT];
  //
  RooAbsPdf* mggBkg[_NCAT];
  RooAbsPdf* mjjBkg[_NCAT];
  vector<double> mean_mjj;
  vector<double> sigma_mjj;
  vector<double> mean_mgg;
  vector<double> sigma_mgg;
  for (int c = 0; c < _NCAT; ++c) {
    // data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c] = (RooDataSet*) _w->data(TString::Format("Sig_cat%d",c));
    mggGaussSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggGaussSig_cat%d",c));
    mggCBSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggCBSig_cat%d",c));
    mggSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggSig_cat%d",c));
    mggBkg[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggBkg_cat%d",c));
    mjjGaussSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjGaussSig_cat%d",c));
    mjjCBSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjCBSig_cat%d",c));
    mjjSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjSig_cat%d",c));
    mjjBkg[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjBkg_cat%d",c));

    double mgg_sigma = (mggSig[c]->sigma(*_w->var("mgg")))->getVal();
    sigma_mgg.push_back(mgg_sigma);

    double mgg_mean = (mggSig[c]->mean(*_w->var("mgg")))->getVal();
    mean_mgg.push_back(mgg_mean);

    double mjj_sigma = (mjjSig[c]->sigma(*_w->var("mjj")))->getVal();
    sigma_mjj.push_back(mjj_sigma);

    double mjj_mean = (mjjSig[c]->mean(*_w->var("mjj")))->getVal();
    mean_mjj.push_back(mjj_mean);

  } // close categories
  RooRealVar* mgg = _w->var("mgg");
  mgg->setUnit("GeV");
  RooAbsPdf* mggGaussSigAll = _w->pdf("mggGaussSig");
  RooAbsPdf* mggCBSigAll = _w->pdf("mggCBSig");
  RooAbsPdf* mggSigAll = _w->pdf("mggSig");
  RooRealVar* mjj = _w->var("mjj");
  mjj->setUnit("GeV");
  RooAbsPdf* mjjGaussSigAll = _w->pdf("mjjGaussSig");
  RooAbsPdf* mjjCBSigAll = _w->pdf("mjjCBSig");
  RooAbsPdf* mjjSigAll = _w->pdf("mjjSig");
  //RooAbsPdf* mggBkgAll = w->pdf("mggBkg_cat1");
  //
  //****************************//
  // Plot mgg Fit results
  //****************************//
  // Set P.D.F. parameter names
  // WARNING: Do not use it if Workspaces are created
  // SetParamNames(w);
  Float_t minSigPlotMgg(115),maxSigPlotMgg(135);
  Float_t minSigPlotMjj(60),maxSigPlotMjj(180);
  mgg->setRange("SigPlotRange",minSigPlotMgg,maxSigPlotMgg);
  mjj->setRange("SigPlotRange",minSigPlotMjj,maxSigPlotMjj);
  Int_t nBinsMass(20); // just need to plot
  RooPlot* plotmggAll = mgg->frame(Range("SigPlotRange"),Bins(nBinsMass));
  signalAll->plotOn(plotmggAll);
  gStyle->SetOptTitle(0);
  TCanvas* c1 = new TCanvas("cMgg","mgg",0,0,500,500);
  c1->cd(1);
  //********************************************//
  // Plot Signal Categories
  //****************************//
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  RooPlot* plotmgg[_NCAT];
  for (int c = 0; c < _NCAT; ++c) {
    plotmgg[c] = mgg->frame(Range("SigPlotRange"),Bins(nBinsMass));
    sigToFit[c]->plotOn(plotmgg[c]);
    mggSig[c] ->plotOn(plotmgg[c]);
    double chi2n = plotmgg[c]->chiSquare(0) ;
    cout << "------------------------- Experimental chi2 = " << chi2n << endl;
    mggSig[c] ->plotOn(
		       plotmgg[c],
		       Components(TString::Format("mggGaussSig_cat%d",c)),
		       LineStyle(kDashed),LineColor(kGreen));
    mggSig[c] ->plotOn(
		       plotmgg[c],
		       Components(TString::Format("mggCBSig_cat%d",c)),
		       LineStyle(kDashed),LineColor(kRed));
    mggSig[c] ->paramOn(plotmgg[c]);
    sigToFit[c] ->plotOn(plotmgg[c]);
    // TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
//    TH1F *hist = new TH1F(TString::Format("histMgg_cat%d",c), "hist", 400, minSigPlotMgg, maxSigPlotMgg);
    //plotmgg[c]->SetTitle("CMS preliminary 19.7/fb ");
    plotmgg[c]->SetMinimum(0.0);
    plotmgg[c]->SetMaximum(1.40*plotmgg[c]->GetMaximum());
    plotmgg[c]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
    TCanvas* ctmp = new TCanvas(TString::Format("ctmpSigMgg_cat%d",c),"Background Categories",0,0,500,500);
    plotmgg[c]->Draw();
    plotmgg[c]->Draw("SAME");
    TLegend *legmc = new TLegend(0.62,0.75,0.99,0.99);
    legmc->AddEntry(plotmgg[c]->getObject(5),"Simulation","LPE");
    legmc->AddEntry(plotmgg[c]->getObject(1),"Parametric Model","L");
    legmc->AddEntry(plotmgg[c]->getObject(2),"Gaussian Outliers","L");
    legmc->AddEntry(plotmgg[c]->getObject(3),"Crystal Ball component","L");
    legmc->SetHeader(" ");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();
    TPaveText *pt = new TPaveText(0.1,0.94,0.7,0.99, "brNDC");
    //pt->SetName("title");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetTextSize(0.035);
    pt->AddText("CMS Preliminary Simulation ");
    pt->Draw();
    // float effS = effSigma(hist);
    TString str_desc;
    if(_sigMass==0)
      str_desc=" Nonresonant HH";
    else
      str_desc=TString::Format(" m_{X} = %d GeV",_sigMass);
    TLatex *lat = new TLatex(
			     minSigPlotMgg+0.5,0.85*plotmgg[c]->GetMaximum(),str_desc);
    lat->Draw();
    TLatex *lat2 = new TLatex(
			      minSigPlotMgg+0.5,0.70*plotmgg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    str_desc=TString::Format(" #mu = %.2f GeV",mean_mgg[c]);
    TLatex *lat3 = new TLatex(
			      minSigPlotMgg+12,0.6*plotmgg[c]->GetMaximum(),str_desc);
    lat3->Draw();
    str_desc=TString::Format(" #sigma = %.2f GeV",sigma_mgg[c]);
    TLatex *lat4 = new TLatex(
			      minSigPlotMgg+12.,0.5*plotmgg[c]->GetMaximum(),str_desc);
    lat4->Draw();
    ///////
    char myChi2buffer[50];
    sprintf(myChi2buffer,"#chi^{2}/ndof = %f",chi2n);
    TLatex* latex = new TLatex(0.52, 0.7, myChi2buffer);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    //latex -> Draw("same");
    ctmp->SaveAs(TString::Format("sigmodelMgg_cat%d.pdf",c));
    ctmp->SaveAs(TString::Format("sigmodelMgg_cat%d.png",c));
    //ctmp->SaveAs(TString::Format("sigmodelMgg_cat%d.C",c));
  } // close categories

  c1 = new TCanvas("cMjj","mgg",0,0,500,500);
  c1->cd(1);
  //********************************************//
  // Plot Signal Categories
  //****************************//
  text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  RooPlot* plotmjj[_NCAT];
  for (int c = 0; c < _NCAT; ++c) {
    plotmjj[c] = mjj->frame(Range("SigPlotRange"),Bins(nBinsMass));
    sigToFit[c]->plotOn(plotmjj[c]);
    mjjSig[c] ->plotOn(plotmjj[c]);
    double chi2n = plotmjj[c]->chiSquare(0) ;
    cout << "------------------------- Experimental chi2 = " << chi2n << endl;
    mjjSig[c] ->plotOn(
		       plotmjj[c],
		       Components(TString::Format("mjjGaussSig_cat%d",c)),
		       LineStyle(kDashed),LineColor(kGreen));
    mjjSig[c] ->plotOn(
		       plotmjj[c],
		       Components(TString::Format("mjjCBSig_cat%d",c)),
		       LineStyle(kDashed),LineColor(kRed));
    mjjSig[c] ->paramOn(plotmjj[c]);
    sigToFit[c] ->plotOn(plotmjj[c]);
    // TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
//    TH1F *hist = new TH1F(TString::Format("histMjj_cat%d",c), "hist", 400, minSigPlotMjj, maxSigPlotMjj);
    //plotmjj[c]->SetTitle("CMS preliminary 19.7/fb ");
    plotmjj[c]->SetMinimum(0.0);
    plotmjj[c]->SetMaximum(1.40*plotmjj[c]->GetMaximum());
    plotmjj[c]->GetXaxis()->SetTitle("m_{jj} (GeV)");
    TCanvas* ctmp = new TCanvas(TString::Format("ctmpSigMjj_cat%d",c),"Background Categories",0,0,500,500);
    plotmjj[c]->Draw();
    plotmjj[c]->Draw("SAME");
    TLegend *legmc = new TLegend(0.62,0.75,0.99,0.99);
    legmc->AddEntry(plotmjj[c]->getObject(5),"Simulation","LPE");
    legmc->AddEntry(plotmjj[c]->getObject(1),"Parametric Model","L");
    legmc->AddEntry(plotmjj[c]->getObject(2),"Gaussian Outliers","L");
    legmc->AddEntry(plotmjj[c]->getObject(3),"Crystal Ball component","L");
    legmc->SetHeader(" ");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();
    TPaveText *pt = new TPaveText(0.1,0.94,0.7,0.99, "brNDC");
    //pt->SetName("title");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetTextSize(0.035);
    pt->AddText("CMS Preliminary Simulation ");
    pt->Draw();
    // float effS = effSigma(hist);
    TString str_desc;
    if(_sigMass==0)
      str_desc=" Nonresonant HH";
    else
      str_desc=TString::Format(" m_{X} = %d GeV",_sigMass);
    TLatex *lat = new TLatex(
			     minSigPlotMjj+0.5,0.85*plotmjj[c]->GetMaximum(),str_desc);
    lat->Draw();
    TLatex *lat2 = new TLatex(
			      minSigPlotMjj+0.5,0.70*plotmjj[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    str_desc=TString::Format(" #mu = %.2f GeV",mean_mjj[c]);
    TLatex *lat3 = new TLatex(
			      minSigPlotMjj+70.5,0.6*plotmjj[c]->GetMaximum(),str_desc);
    lat3->Draw();
    str_desc=TString::Format(" #sigma = %.2f GeV",sigma_mjj[c]);
    TLatex *lat4 = new TLatex(
			      minSigPlotMjj+70.5,0.5*plotmjj[c]->GetMaximum(),str_desc);
    lat4->Draw();
    ///////
    char myChi2buffer[50];
    sprintf(myChi2buffer,"#chi^{2}/ndof = %f",chi2n);
    TLatex* latex = new TLatex(0.52, 0.7, myChi2buffer);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    //latex -> Draw("same");
    ctmp->SaveAs(TString::Format("sigmodelMjj_cat%d.pdf",c));
    ctmp->SaveAs(TString::Format("sigmodelMjj_cat%d.png",c));
    //ctmp->SaveAs(TString::Format("sigmodelMjj_cat%d.C",c));
  } // close categories

  return;
} // close makeplots signal

void bbgg2DFitter::MakePlotsHiggs(float mass)
{
  std::vector<TString> catdesc;
  if ( _NCAT == 2 ){
    catdesc.push_back(" High Purity");
    catdesc.push_back(" Med. Purity");
  }
  else{
    catdesc.push_back(" #splitline{High Purity}{High m_{#gamma#gammajj}^{kin}}");
    catdesc.push_back(" #splitline{Med. Purity}{High m_{#gamma#gammajj}^{kin}}");
    catdesc.push_back(" #splitline{High Purity}{Low m_{#gamma#gammajj}^{kin}}");
    catdesc.push_back(" #splitline{Med. Purity}{Low m_{#gamma#gammajj}^{kin}}");
  }
  // retrieve data sets from the workspace
  // RooDataSet* dataAll = (RooDataSet*) w->data("Data");
  //RooDataSet* signalAll = (RooDataSet*) w->data("Sig");
  //RooDataSet* higgsAll = (RooDataSet*) w->data("Hig");
  // blinded dataset
  // RooDataSet* data[ncat];
  TString component[5] = {"ggH","ttH","VBF","VH", "bbH"};

  //    higToFit[c] = (RooDataSet*) w->data(TString::Format("Hig_%d_cat%d",higgschannel,c));

  for (int d = 0; d < 5; ++d){

    RooDataSet* higToFit[_NCAT];
    RooAbsPdf* mggGaussSig[_NCAT];
    RooAbsPdf* mggCBSig[_NCAT];
    RooAbsPdf* mggSig[_NCAT];
    RooAbsPdf* mjjGaussSig[_NCAT];
    RooAbsPdf* mjjCBSig[_NCAT];
    RooAbsPdf* mjjSig[_NCAT];
    //
    RooAbsPdf* mggBkg[_NCAT];
    RooAbsPdf* mjjBkg[_NCAT];
    for (int c = 0; c < _NCAT; ++c) {
      // data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
      higToFit[c] = (RooDataSet*) _w->data(TString::Format("Hig_%d_cat%d",d,c));
      mggGaussSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggGaussHig_%d_cat%d",d,c));
      mggCBSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggCBHig_%d_cat%d",d,c));
      mggSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggHig_%d_cat%d",d,c));
      mggBkg[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggBkg_%d_cat%d",d,c));
      if(d == 1 || d == 3){
         mjjGaussSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjGaussHig_%d_cat%d",d,c));
         mjjCBSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjCBHig_%d_cat%d",d,c));
         mjjSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_%d_cat%d",d,c));
      }else{
         mjjSig[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_pol0_%d_cat%d",d,c));
      }
      mjjBkg[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjBkg_%d_cat%d",d,c));
    } // close categories
    

    RooRealVar* mgg = _w->var("mgg");
    mgg->setUnit("GeV");
    RooRealVar* mjj = _w->var("mjj");
    mjj->setUnit("GeV");
    //RooAbsPdf* mggBkgAll = w->pdf("mggBkg_cat1");
    //
    //****************************//
    // Plot mgg Fit results
    //****************************//
    // Set P.D.F. parameter names
    // WARNING: Do not use it if Workspaces are created
    // SetParamNames(w);
    Float_t minHigPlotMgg(115),maxHigPlotMgg(135);
    Float_t minHigPlotMjj(60),maxHigPlotMjj(180);
    Int_t nBinsMass(20); // just need to plot
    //RooPlot* plotmggAll = mgg->frame(Range(minSigFit,maxSigFit),Bins(nBinsMass));
    //higgsAll->plotOn(plotmggAll);
    gStyle->SetOptTitle(0);
    TCanvas* c1 = new TCanvas(TString::Format("cMgg_%d",d),"mgg",0,0,500,500);
    c1->cd(1);
    //********************************************//
    // Plot Signal Categories
    //****************************//
    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.04);
    RooPlot* plotmgg[_NCAT];

    for (int c = 0; c < _NCAT; ++c) {
      plotmgg[c] = mgg->frame(Range("HigFitRange"),Bins(nBinsMass));
      higToFit[c]->plotOn(plotmgg[c]);
      mggSig[c] ->plotOn(plotmgg[c]);
      double chi2n = plotmgg[c]->chiSquare(0) ;
      cout << "------------------------- Experimental chi2 = " << chi2n << endl;
      mggSig[c] ->plotOn(
			 plotmgg[c],
			 Components(TString::Format("mggGaussHig_%d_cat%d",d,c)),
			 LineStyle(kDashed),LineColor(kGreen));
      mggSig[c] ->plotOn(
			 plotmgg[c],
			 Components(TString::Format("mggCBHig_%d_cat%d",d,c)),
			 LineStyle(kDashed),LineColor(kRed));
      mggSig[c] ->paramOn(plotmgg[c]);
      higToFit[c] ->plotOn(plotmgg[c]);
      // TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
//      TH1F *hist = new TH1F(TString::Format("histMgg_%d_cat%d",d,c), "hist", 400, minHigPlotMgg, maxHigPlotMgg);
      //plotmgg[c]->SetTitle("CMS Preliminary 19.7/fb ");
      plotmgg[c]->SetMinimum(0.0);
      plotmgg[c]->SetMaximum(1.40*plotmgg[c]->GetMaximum());
      plotmgg[c]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
      TCanvas* ctmp = new TCanvas(TString::Format("ctmpHigMgg_%d_cat%d",d,c),"Background Categories",0,0,500,500);
      plotmgg[c]->Draw();
      plotmgg[c]->Draw("SAME");
      TLegend *legmc = new TLegend(0.62,0.75,0.99,0.99);
      legmc->AddEntry(plotmgg[c]->getObject(5),component[d],"LPE");
      legmc->AddEntry(plotmgg[c]->getObject(1),"Parametric Model","L");
      legmc->AddEntry(plotmgg[c]->getObject(2),"Gaussian Outliers","L");
      legmc->AddEntry(plotmgg[c]->getObject(3),"Crystal Ball component","L");
      legmc->SetHeader(" ");
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();
      TPaveText *pt = new TPaveText(0.1,0.94,0.7,0.99, "brNDC");
      //pt->SetName("title");
      pt->SetBorderSize(0);
      pt->SetFillColor(0);
      pt->SetTextSize(0.035);
      pt->AddText("CMS Preliminary Simulation ");
      pt->Draw();
      // float effS = effSigma(hist);
      TString str_desc;
      if(_sigMass==0)
	      str_desc=" Nonresonant HH";
      else
	      str_desc=TString::Format(" m_{X} = %d GeV",_sigMass);
      TLatex *lat = new TLatex(
			       minHigPlotMgg+0.5,0.85*plotmgg[c]->GetMaximum(),str_desc);
      lat->Draw();
      TLatex *lat2 = new TLatex(
				minHigPlotMgg+0.5,0.70*plotmgg[c]->GetMaximum(),catdesc.at(c));
      lat2->Draw();
      ///////
      char myChi2buffer[50];
      sprintf(myChi2buffer,"#chi^{2}/ndof = %f",chi2n);
      TLatex* latex = new TLatex(0.52, 0.7, myChi2buffer);
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      //latex -> Draw("same");
      ctmp->SaveAs(TString::Format("higmodelMgg_%d_cat%d.pdf",d,c));
      ctmp->SaveAs(TString::Format("higmodelMgg_%d_cat%d.png",d,c));
      //ctmp->SaveAs(TString::Format("sigmodelMgg_cat%d.C",c));
    } // close categories

    c1 = new TCanvas(TString::Format("cMjj_%d",d),"mjj",0,0,500,500);
    c1->cd(1);
    //********************************************//
    // Plot Signal Categories
    //****************************//
    text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.04);
    RooPlot* plotmjj[_NCAT];

    for (int c = 0; c < _NCAT; ++c) {
      plotmjj[c] = mjj->frame(Range("HigFitRange"),Bins(nBinsMass));
      higToFit[c]->plotOn(plotmjj[c]);
      mjjSig[c] ->plotOn(plotmjj[c]);
      double chi2n = plotmjj[c]->chiSquare(0) ;
      cout << "------------------------- Experimental chi2 = " << chi2n << endl;
      if(d == 1 || d == 3){
         mjjSig[c] ->plotOn(
			 plotmjj[c],
			 Components(TString::Format("mjjGaussHig_%d_cat%d",d,c)),
			 LineStyle(kDashed),LineColor(kGreen));
         mjjSig[c] ->plotOn(
			 plotmjj[c],
			 Components(TString::Format("mjjCBHig_%d_cat%d",d,c)),
			 LineStyle(kDashed),LineColor(kRed));
         mjjSig[c] ->paramOn(plotmjj[c]);
         higToFit[c] ->plotOn(plotmjj[c]);
         // TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
//         TH1F *hist = new TH1F(TString::Format("histMjj_%d_cat%d",d,c), "hist", 400, minHigPlotMjj, maxHigPlotMjj);
         //plotmjj[c]->SetTitle("CMS preliminary 19.7/fb ");
         plotmjj[c]->SetMinimum(0.0);
         plotmjj[c]->SetMaximum(1.40*plotmjj[c]->GetMaximum());
         plotmjj[c]->GetXaxis()->SetTitle("m_{jj} (GeV)");
         TCanvas* ctmp = new TCanvas(TString::Format("ctmpHigMjj_%d_cat_%d",d,c),"Background Categories",0,0,500,500);
         plotmjj[c]->Draw();
         plotmjj[c]->Draw("SAME");
         TLegend *legmc = new TLegend(0.62,0.75,0.99,0.99);

         legmc->AddEntry(plotmjj[c]->getObject(5),component[d],"LPE");
         legmc->AddEntry(plotmjj[c]->getObject(1),"Parametric Model","L");
         legmc->AddEntry(plotmjj[c]->getObject(2),"Gaussian Outliers","L");
         legmc->AddEntry(plotmjj[c]->getObject(3),"Crystal Ball component","L");
         legmc->SetHeader(" ");
         legmc->SetBorderSize(0);
         legmc->SetFillStyle(0);
         legmc->Draw();
         TPaveText *pt = new TPaveText(0.1,0.94,0.7,0.99, "brNDC");
         //pt->SetName("title");
         pt->SetBorderSize(0);
         pt->SetFillColor(0);
         pt->SetTextSize(0.035);
         pt->AddText("CMS Preliminary Simulation ");
         pt->Draw();
         // float effS = effSigma(hist);
         TString str_desc;
         if(_sigMass==0)
	         str_desc=" Nonresonant HH";
         else
	         str_desc=TString::Format(" m_{X} = %d GeV",_sigMass);
         TLatex *lat = new TLatex(
			       minHigPlotMjj+0.5,0.85*plotmjj[c]->GetMaximum(),str_desc);
         lat->Draw();
         TLatex *lat2 = new TLatex(
				minHigPlotMjj+0.5,0.70*plotmjj[c]->GetMaximum(),catdesc.at(c));
         lat2->Draw();
         ///////
         char myChi2buffer[50];
         sprintf(myChi2buffer,"#chi^{2}/ndof = %f",chi2n);
         TLatex* latex = new TLatex(0.52, 0.7, myChi2buffer);
         latex -> SetNDC();
         latex -> SetTextFont(42);
         latex -> SetTextSize(0.04);
         //latex -> Draw("same");
         ctmp->SaveAs(TString::Format("higmodelMjj_%d_cat%d.pdf",d,c));
         ctmp->SaveAs(TString::Format("higmodelMjj_%d_cat%d.png",d,c));
      }else{
         //mjjSig[c] ->paramOn(plotmjj[c]);
         higToFit[c] ->plotOn(plotmjj[c]);
         // TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
//         TH1F *hist = new TH1F(TString::Format("histMjj_%d_cat%d",d,c), "hist", 400, minHigPlotMjj, maxHigPlotMjj);
         //plotmjj[c]->SetTitle("CMS preliminary 19.7/fb ");
         plotmjj[c]->SetMinimum(0.0);
         plotmjj[c]->SetMaximum(1.40*plotmjj[c]->GetMaximum());
         plotmjj[c]->GetXaxis()->SetTitle("m_{jj} (GeV)");
         TCanvas* ctmp = new TCanvas(TString::Format("ctmpHigMjj_%d_cat_%d",d,c),"Background Categories",0,0,500,500);
         plotmjj[c]->Draw();
         plotmjj[c]->Draw("SAME");
         TLegend *legmc = new TLegend(0.62,0.75,0.99,0.99);
         legmc->AddEntry(plotmjj[c]->getObject(5),component[d],"LPE");
         legmc->AddEntry(plotmjj[c]->getObject(1),"Parametric Model","L");
         legmc->SetHeader(" ");
         legmc->SetBorderSize(0);
         legmc->SetFillStyle(0);
         legmc->Draw();
         TPaveText *pt = new TPaveText(0.1,0.94,0.7,0.99, "brNDC");
         //pt->SetName("title");
         pt->SetBorderSize(0);
         pt->SetFillColor(0);
         pt->SetTextSize(0.035);
         pt->AddText("CMS Preliminary Simulation ");
         pt->Draw();
         // float effS = effSigma(hist);
         TString str_desc;
         if(_sigMass==0)
	    str_desc=" Nonresonant HH";
         else
	    str_desc=TString::Format(" m_{X} = %d GeV",_sigMass);
         TLatex *lat = new TLatex(
			       minHigPlotMjj+0.5,0.85*plotmjj[c]->GetMaximum(),str_desc);
         lat->Draw();
         TLatex *lat2 = new TLatex(
				minHigPlotMjj+0.5,0.70*plotmjj[c]->GetMaximum(),catdesc.at(c));
         lat2->Draw();
         ///////
         char myChi2buffer[50];
         sprintf(myChi2buffer,"#chi^{2}/ndof = %f",chi2n);
         TLatex* latex = new TLatex(0.52, 0.7, myChi2buffer);
         latex -> SetNDC();
         latex -> SetTextFont(42);
         latex -> SetTextSize(0.04);
         //latex -> Draw("same");
         ctmp->SaveAs(TString::Format("higmodelMjj_%d_cat%d.pdf",d,c));
         ctmp->SaveAs(TString::Format("higmodelMjj_%d_cat%d.png",d,c));
      }
      //ctmp->SaveAs(TString::Format("sigmodelMjj_cat%d.C",c));
    } // close categories

  } // close to higgs component
  return;
} // close makeplots signal

void bbgg2DFitter::MakeSigWS(const char* fileBaseName)
{
  TString wsDir = "workspaces/";
  //**********************************************************************//
  // Write pdfs and datasets into the workspace before to save
  // for statistical tests.
  //**********************************************************************//
  RooAbsPdf* SigPdf[_NCAT];
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  for (int c = 0; c < _NCAT; ++c) {
    SigPdf[c] = (RooAbsPdf*) _w->pdf(TString::Format("SigPdf_cat%d",c));
    wAll->import(*_w->pdf(TString::Format("SigPdf_cat%d",c)));
  }
  // (2) Systematics on energy scale and resolution
  // 1,1,1 statistical to be treated on the datacard
  wAll->factory("CMS_hgg_sig_m0_absShift[1,1,1]");
  wAll->factory("prod::CMS_hgg_sig_m0_cat0(mgg_sig_m0_cat0, CMS_hgg_sig_m0_absShift)");
  wAll->factory("prod::CMS_hgg_sig_m0_cat1(mgg_sig_m0_cat1, CMS_hgg_sig_m0_absShift)");
  if ( _NCAT > 2 ){
    wAll->factory("prod::CMS_hgg_sig_m0_cat2(mgg_sig_m0_cat2, CMS_hgg_sig_m0_absShift)");
    wAll->factory("prod::CMS_hgg_sig_m0_cat3(mgg_sig_m0_cat3, CMS_hgg_sig_m0_absShift)");
  }
  wAll->factory("CMS_hbb_sig_m0_absShift[1,1,1]");
  wAll->factory("prod::CMS_hbb_sig_m0_cat0(mjj_sig_m0_cat0, CMS_hbb_sig_m0_absShift)");
  wAll->factory("prod::CMS_hbb_sig_m0_cat1(mjj_sig_m0_cat1, CMS_hbb_sig_m0_absShift)");
  if ( _NCAT > 2 ){
    wAll->factory("prod::CMS_hbb_sig_m0_cat2(mjj_sig_m0_cat2, CMS_hbb_sig_m0_absShift)");
    wAll->factory("prod::CMS_hbb_sig_m0_cat3(mjj_sig_m0_cat3, CMS_hbb_sig_m0_absShift)");
  }
  // (3) Systematics on resolution
  wAll->factory("CMS_hgg_sig_sigmaScale[1,1,1]");
  wAll->factory("prod::CMS_hgg_sig_sigma_cat0(mgg_sig_sigma_cat0, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_sigma_cat1(mgg_sig_sigma_cat1, CMS_hgg_sig_sigmaScale)");
  if ( _NCAT > 2 ){
    wAll->factory("prod::CMS_hgg_sig_sigma_cat2(mgg_sig_sigma_cat2, CMS_hgg_sig_sigmaScale)");
    wAll->factory("prod::CMS_hgg_sig_sigma_cat3(mgg_sig_sigma_cat3, CMS_hgg_sig_sigmaScale)");
  }
  wAll->factory("prod::CMS_hgg_sig_gsigma_cat0(mgg_sig_gsigma_cat0, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_gsigma_cat1(mgg_sig_gsigma_cat1, CMS_hgg_sig_sigmaScale)");
  if ( _NCAT > 2 ){
    wAll->factory("prod::CMS_hgg_sig_gsigma_cat2(mgg_sig_gsigma_cat2, CMS_hgg_sig_sigmaScale)");
    wAll->factory("prod::CMS_hgg_sig_gsigma_cat3(mgg_sig_gsigma_cat3, CMS_hgg_sig_sigmaScale)");
  }
  wAll->factory("CMS_hbb_sig_sigmaScale[1,1,1]");
  wAll->factory("prod::CMS_hbb_sig_sigma_cat0(mjj_sig_sigma_cat0, CMS_hbb_sig_sigmaScale)");
  wAll->factory("prod::CMS_hbb_sig_sigma_cat1(mjj_sig_sigma_cat1, CMS_hbb_sig_sigmaScale)");
  if ( _NCAT > 2 ){
    wAll->factory("prod::CMS_hbb_sig_sigma_cat2(mjj_sig_sigma_cat2, CMS_hbb_sig_sigmaScale)");
    wAll->factory("prod::CMS_hbb_sig_sigma_cat3(mjj_sig_sigma_cat3, CMS_hbb_sig_sigmaScale)");
  }
  wAll->factory("prod::CMS_hbb_sig_gsigma_cat0(mjj_sig_gsigma_cat0, CMS_hbb_sig_sigmaScale)");
  wAll->factory("prod::CMS_hbb_sig_gsigma_cat1(mjj_sig_gsigma_cat1, CMS_hbb_sig_sigmaScale)");
  if ( _NCAT > 2 ){
    wAll->factory("prod::CMS_hbb_sig_gsigma_cat2(mjj_sig_gsigma_cat2, CMS_hbb_sig_sigmaScale)");
    wAll->factory("prod::CMS_hbb_sig_gsigma_cat3(mjj_sig_gsigma_cat3, CMS_hbb_sig_sigmaScale)");
  }
  // (4) do reparametrization of signal
  for (int c = 0; c < _NCAT; ++c) wAll->factory(
					       TString::Format("EDIT::CMS_sig_cat%d(SigPdf_cat%d,",c,c) +
					       TString::Format(" mgg_sig_m0_cat%d=CMS_hgg_sig_m0_cat%d,", c,c) +
					       TString::Format(" mgg_sig_sigma_cat%d=CMS_hgg_sig_sigma_cat%d,", c,c) +
					       TString::Format(" mgg_sig_gsigma_cat%d=CMS_hgg_sig_gsigma_cat%d,", c,c) +
					       TString::Format(" mjj_sig_m0_cat%d=CMS_hbb_sig_m0_cat%d,", c,c) +
					       TString::Format(" mjj_sig_sigma_cat%d=CMS_hbb_sig_sigma_cat%d,", c,c) +
					       TString::Format(" mjj_sig_gsigma_cat%d=CMS_hbb_sig_gsigma_cat%d)", c,c)
					       );
  TString filename(wsDir+TString(fileBaseName)+".inputsig.root");
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;
  return;
} // close make signal WP

void bbgg2DFitter::MakeHigWS(const char* fileHiggsName,int higgschannel)
{
  TString wsDir = "workspaces/";
  //**********************************************************************//
  // Write pdfs and datasets into the workspace before to save to a file
  // for statistical tests.
  //**********************************************************************//
  RooAbsPdf* HigPdf[_NCAT];
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  for (int c = 0; c < _NCAT; ++c) {
    HigPdf[c] = (RooAbsPdf*) _w->pdf(TString::Format("HigPdf_%d_cat%d",higgschannel,c));
    wAll->import(*_w->pdf(TString::Format("HigPdf_%d_cat%d",higgschannel,c)));
  }
  // (2) Systematics on energy scale and resolution
  // 1,1,1 statistical to be treated on the datacard
  wAll->factory("CMS_hgg_sig_m0_absShift[1,1,1]");
  wAll->factory(TString::Format("prod::CMS_hgg_hig_m0_%d_cat0(mgg_hig_m0_%d_cat0, CMS_hgg_sig_m0_absShift)",higgschannel,higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_m0_%d_cat1(mgg_hig_m0_%d_cat1, CMS_hgg_sig_m0_absShift)",higgschannel,higgschannel));
  if ( _NCAT > 2 ){
       wAll->factory(TString::Format("prod::CMS_hgg_hig_m0_%d_cat2(mgg_hig_m0_%d_cat2, CMS_hgg_sig_m0_absShift)",higgschannel,higgschannel));
       wAll->factory(TString::Format("prod::CMS_hgg_hig_m0_%d_cat3(mgg_hig_m0_%d_cat3, CMS_hgg_sig_m0_absShift)",higgschannel,higgschannel));
  }
  if(higgschannel == 1 || higgschannel == 3){
     wAll->factory("CMS_hbb_sig_m0_absShift[1,1,1]");
     wAll->factory(TString::Format("prod::CMS_hbb_hig_m0_%d_cat0(mjj_hig_m0_%d_cat0, CMS_hbb_sig_m0_absShift)",higgschannel,higgschannel));
     wAll->factory(TString::Format("prod::CMS_hbb_hig_m0_%d_cat1(mjj_hig_m0_%d_cat1, CMS_hbb_sig_m0_absShift)",higgschannel,higgschannel));
     if ( _NCAT > 2 ){
          wAll->factory(TString::Format("prod::CMS_hbb_hig_m0_%d_cat2(mjj_hig_m0_%d_cat2, CMS_hbb_sig_m0_absShift)",higgschannel,higgschannel));
          wAll->factory(TString::Format("prod::CMS_hbb_hig_m0_%d_cat3(mjj_hig_m0_%d_cat3, CMS_hbb_sig_m0_absShift)",higgschannel,higgschannel));
     }
  }
  // (3) Systematics on resolution
  wAll->factory("CMS_hgg_sig_sigmaScale[1,1,1]");
  wAll->factory(TString::Format("prod::CMS_hgg_hig_sigma_%d_cat0(mgg_hig_sigma_%d_cat0, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_sigma_%d_cat1(mgg_hig_sigma_%d_cat1, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
  if ( _NCAT > 2 ){
       wAll->factory(TString::Format("prod::CMS_hgg_hig_sigma_%d_cat2(mgg_hig_sigma_%d_cat2, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
       wAll->factory(TString::Format("prod::CMS_hgg_hig_sigma_%d_cat3(mgg_hig_sigma_%d_cat3, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
  }
  wAll->factory(TString::Format("prod::CMS_hgg_hig_gsigma_%d_cat0(mgg_hig_gsigma_%d_cat0, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
  wAll->factory(TString::Format("prod::CMS_hgg_hig_gsigma_%d_cat1(mgg_hig_gsigma_%d_cat1, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
  if ( _NCAT > 2 ){
       wAll->factory(TString::Format("prod::CMS_hgg_hig_gsigma_%d_cat2(mgg_hig_gsigma_%d_cat2, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
       wAll->factory(TString::Format("prod::CMS_hgg_hig_gsigma_%d_cat3(mgg_hig_gsigma_%d_cat3, CMS_hgg_sig_sigmaScale)",higgschannel,higgschannel));
  }
  if(higgschannel == 1 || higgschannel == 3){
        wAll->factory("CMS_hbb_sig_sigmaScale[1,1,1]");
        wAll->factory(TString::Format("prod::CMS_hbb_hig_sigma_%d_cat0(mjj_hig_sigma_%d_cat0, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
        wAll->factory(TString::Format("prod::CMS_hbb_hig_sigma_%d_cat1(mjj_hig_sigma_%d_cat1, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
        if ( _NCAT > 2 ){
             wAll->factory(TString::Format("prod::CMS_hbb_hig_sigma_%d_cat2(mjj_hig_sigma_%d_cat2, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
             wAll->factory(TString::Format("prod::CMS_hbb_hig_sigma_%d_cat3(mjj_hig_sigma_%d_cat3, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
        }
        wAll->factory(TString::Format("prod::CMS_hbb_hig_gsigma_%d_cat0(mjj_hig_gsigma_%d_cat0, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
        wAll->factory(TString::Format("prod::CMS_hbb_hig_gsigma_%d_cat1(mjj_hig_gsigma_%d_cat1, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
        if ( _NCAT > 2 ){
             wAll->factory(TString::Format("prod::CMS_hbb_hig_gsigma_%d_cat2(mjj_hig_gsigma_%d_cat2, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
             wAll->factory(TString::Format("prod::CMS_hbb_hig_gsigma_%d_cat3(mjj_hig_gsigma_%d_cat3, CMS_hbb_sig_sigmaScale)",higgschannel,higgschannel));
        }
  }
  // (4) do reparametrization of signal
  if(higgschannel == 1 || higgschannel == 3){
     for (int c = 0; c < _NCAT; ++c) wAll->factory(
					       TString::Format("EDIT::CMS_hig_%d_cat%d(HigPdf_%d_cat%d,",higgschannel,c,higgschannel,c) +
					       TString::Format(" mgg_hig_m0_%d_cat%d=CMS_hgg_hig_m0_%d_cat%d,",higgschannel, c,higgschannel,c) +
					       TString::Format(" mgg_hig_sigma_%d_cat%d=CMS_hgg_hig_sigma_%d_cat%d,",higgschannel, c,higgschannel,c) +
					       TString::Format(" mgg_hig_gsigma_%d_cat%d=CMS_hgg_hig_gsigma_%d_cat%d,",higgschannel, c,higgschannel,c) +
					       TString::Format(" mjj_hig_m0_%d_cat%d=CMS_hbb_hig_m0_%d_cat%d,",higgschannel, c,higgschannel,c) +
					       TString::Format(" mjj_hig_sigma_%d_cat%d=CMS_hbb_hig_sigma_%d_cat%d,",higgschannel, c,higgschannel,c) +
					       TString::Format(" mjj_hig_gsigma_%d_cat%d=CMS_hbb_hig_gsigma_%d_cat%d)",higgschannel, c,higgschannel,c)
					       );
  }else{
     for (int c = 0; c < _NCAT; ++c) wAll->factory(
					       TString::Format("EDIT::CMS_hig_%d_cat%d(HigPdf_%d_cat%d,",higgschannel,c,higgschannel,c) +
					       TString::Format(" mgg_hig_m0_%d_cat%d=CMS_hgg_hig_m0_%d_cat%d,",higgschannel, c,higgschannel,c) +
					       TString::Format(" mgg_hig_sigma_%d_cat%d=CMS_hgg_hig_sigma_%d_cat%d,",higgschannel, c,higgschannel,c) +
					       TString::Format(" mgg_hig_gsigma_%d_cat%d=CMS_hgg_hig_gsigma_%d_cat%d)",higgschannel, c,higgschannel,c)
					       );
  }

  TString filename(wsDir+TString(fileHiggsName)+".inputhig.root");
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;
  return;
} // close make higgs WP

void bbgg2DFitter::MakeBkgWS(const char* fileBaseName)
{
  TString wsDir = "workspaces/";

  //**********************************************************************//
  // Write pdfs and datasets into the workspace before to save to a file
  // for statistical tests.
  //**********************************************************************//
  RooDataSet* data[_NCAT];
  RooAbsPdf* BkgPdf[_NCAT];
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  for (int c = 0; c < _NCAT; ++c) {
    cout<<"here"<<endl;
    data[c] = (RooDataSet*) _w->data(TString::Format("Data_cat%d",c));

    //RooDataHist* dataBinned = data[c]->binnedClone(); // Uncomment if you want to use wights in the limits

    BkgPdf[c] = (RooAbsPdf*) _w->pdf(TString::Format("BkgPdf_cat%d",c));
    wAll->import(*data[c], Rename(TString::Format("data_obs_cat%d",c)));// Comment if you want to use wights in the limits

    //wAll->import(*dataBinned, Rename(TString::Format("data_obs_cat%d",c))); // Uncomment if you want to use wights in the limits

    cout<<"here"<<endl;
    wAll->import(*_w->pdf(TString::Format("BkgPdf_cat%d",c)));
    wAll->import(*_w->var(TString::Format("bkg_8TeV_norm_cat%d",c)));
    cout<<"here"<<endl;
    if(_sigMass == 0 || (_sigMass != 0 && c==1)){
     wAll->factory(
		  TString::Format("CMS_bkg_8TeV_cat%d_norm[%g,0.,100000.0]",
				  c, _w->var(TString::Format("bkg_8TeV_norm_cat%d",c))->getVal()));
    cout<<"here: " << c << " " << _w->var(TString::Format("bkg_8TeV_norm_cat%d",c))->getVal() <<endl;
    wAll->factory(TString::Format("CMS_hgg_bkg_8TeV_slope1_cat%d[%g,-100.,100.]",c, _w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c))->getVal()));
    wAll->factory(TString::Format("CMS_hbb_bkg_8TeV_slope1_cat%d[%g,-100.,100.]",c, _w->var(TString::Format("mjj_bkg_8TeV_slope1_cat%d",c))->getVal()));

     cout<<"here"<<endl;
    }
    else if(_sigMass != 0 && c==0){
     wAll->factory(
		  TString::Format("CMS_bkg_8TeV_cat%d_norm[%g 0.,100000.0]",
				  c, _w->var(TString::Format("bkg_8TeV_norm_cat%d",c))->getVal()));
    wAll->factory(TString::Format("CMS_hgg_bkg_8TeV_slope2_cat%d[%g,-100,100]",c, _w->var(TString::Format("mgg_bkg_8TeV_slope2_cat%d",c))->getVal()));
    wAll->factory(TString::Format("CMS_hgg_bkg_8TeV_slope3_cat%d[%g,-100,100]",c, _w->var(TString::Format("mgg_bkg_8TeV_slope3_cat%d",c))->getVal()));
    wAll->factory(TString::Format("CMS_hbb_bkg_8TeV_slope2_cat%d[%g,-100,100]",c, _w->var(TString::Format("mjj_bkg_8TeV_slope2_cat%d",c))->getVal()));
    wAll->factory(TString::Format("CMS_hbb_bkg_8TeV_slope3_cat%d[%g,-100,100]",c, _w->var(TString::Format("mjj_bkg_8TeV_slope3_cat%d",c))->getVal()));

     cout<<"here"<<endl;
    }
  } // close ncat
  // (2) do reparametrization of background
  for (int c = 0; c < _NCAT; ++c){
   if(_sigMass == 0 || (_sigMass != 0 && c==1)){
    wAll->factory(
		  TString::Format("EDIT::CMS_bkg_8TeV_cat%d(BkgPdf_cat%d,",c,c) +
		  TString::Format(" bkg_8TeV_norm_cat%d=CMS_bkg_8TeV_cat%d_norm,", c,c)+
		  TString::Format(" mgg_bkg_8TeV_slope1_cat%d=CMS_hgg_bkg_8TeV_slope1_cat%d,", c,c) +
		  TString::Format(" mjj_bkg_8TeV_slope1_cat%d=CMS_hbb_bkg_8TeV_slope1_cat%d)", c,c)  );
   }
   else if(_sigMass != 0 && c==0){
    wAll->factory(
          TString::Format("EDIT::CMS_bkg_8TeV_cat%d(BkgPdf_cat%d,",c,c) +
		  TString::Format(" bkg_8TeV_norm_cat%d=CMS_bkg_8TeV_cat%d_norm,", c,c)+
		  TString::Format(" mgg_bkg_8TeV_slope2_cat%d=CMS_hgg_bkg_8TeV_slope2_cat%d,", c,c) +
          TString::Format(" mgg_bkg_8TeV_slope3_cat%d=CMS_hgg_bkg_8TeV_slope3_cat%d,", c,c) +
		  TString::Format(" mjj_bkg_8TeV_slope2_cat%d=CMS_hbb_bkg_8TeV_slope2_cat%d,", c,c) +
          TString::Format(" mjj_bkg_8TeV_slope3_cat%d=CMS_hbb_bkg_8TeV_slope3_cat%d)", c,c)  );
   }
  } // close for cat

  TString filename(wsDir+TString(fileBaseName)+".root");
  wAll->writeToFile(filename);
  cout << "Write background workspace in: " << filename << " file" << endl;
  std::cout << "observation ";
  for (int c = 0; c < _NCAT; ++c) {
    std::cout << " " << wAll->data(TString::Format("data_obs_cat%d",c))->sumEntries();
  }
  std::cout << std::endl;
  return;
} // close make BKG workspace

void bbgg2DFitter::MakeDataCard(const char* fileBaseName, const char* fileBkgName , 
                    const char* fileHiggsNameggh, const char* fileHiggsNametth,
                    const char* fileHiggsNamevbf, const char* fileHiggsNamevh, const char* fileHiggsNamebbh,
                    Bool_t useSigTheoryUnc)
{
  TString cardDir = "datacards/";
  RooDataSet* data[_NCAT];
  RooDataSet* sigToFit[_NCAT];
  RooDataSet* higToFitggh[_NCAT];
  RooDataSet* higToFittth[_NCAT];
  RooDataSet* higToFitvbf[_NCAT];
  RooDataSet* higToFitvh[_NCAT];
  RooDataSet* higToFitbbh[_NCAT];
  for (int c = 0; c < _NCAT; ++c) {
    data[c] = (RooDataSet*) _w->data(TString::Format("Data_cat%d",c));
    sigToFit[c] = (RooDataSet*) _w->data(TString::Format("Sig_cat%d",c));
    //
    higToFitggh[c] = (RooDataSet*) _w->data(TString::Format("Hig_%d_cat%d",0,c));
    higToFittth[c] = (RooDataSet*) _w->data(TString::Format("Hig_%d_cat%d",1,c));
    higToFitvbf[c] = (RooDataSet*) _w->data(TString::Format("Hig_%d_cat%d",2,c));
    higToFitvh[c] = (RooDataSet*) _w->data(TString::Format("Hig_%d_cat%d",3,c));
    higToFitbbh[c] = (RooDataSet*) _w->data(TString::Format("Hig_%d_cat%d",4,c));
  } // close cat
  ////////////////////////////////////////////////////////////////////////////////////
  //RooRealVar* lumi = w->var("lumi");
  cout << "======== Expected Events Number =====================" << endl;
  cout << ".........Measured Data for L = " << _lumi << " pb-1 ............................" << endl;
  if(!_doblinding){ cout << "#Events data: " << _w->data("Data")->sumEntries() << endl; }
  else cout << "#Events data: -1 " << endl;

  if(!_doblinding){
     for (int c = 0; c < _NCAT; ++c) {
          cout << TString::Format("#Events data cat%d: ",c) << data[c]->sumEntries() << endl;
     }
  }else{
     for (int c = 0; c < _NCAT; ++c) {
          cout << TString::Format("#Events data cat%d: ",c) << -1 << endl;
     }
  }
  cout << ".........Expected Signal for L = " << _lumi << " pb-1 ............................" << endl;
  if(!_doblinding){ cout << "#Events Signal: " << _w->data("Data")->sumEntries() << endl; }
  else cout << "#Events Signal: -1 " << endl;
  Float_t siglikeErr[_NCAT];
  for (int c = 0; c < _NCAT; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << sigToFit[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*sigToFit[c]->sumEntries();
  }
  cout << "====================================================" << endl;
  TString filename(cardDir+TString(fileBaseName)+".txt");
  ofstream outFile(filename);

  // outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() << " pb-1 " << endl;
  outFile << "#Run with: combine -d hgg.mH350.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0 -b 3500 -i 50000 --optimizeSim=1 --tries 30" << endl;
  outFile << "# Lumi = " << _lumi << " pb-1" << endl;
  outFile << "imax "<<_NCAT << endl;
  outFile << "jmax 6" << endl; // number of BKG
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;
  outFile << "shapes data_obs cat0 " << TString(fileBkgName)+".root" << " w_all:data_obs_cat0" << endl;
  outFile << "shapes data_obs cat1 "<< TString(fileBkgName)+".root" << " w_all:data_obs_cat1" << endl;
  if ( _NCAT > 2 ){
  outFile << "shapes data_obs cat2 " << TString(fileBkgName)+".root" << " w_all:data_obs_cat2" << endl;
  outFile << "shapes data_obs cat3 "<< TString(fileBkgName)+".root" << " w_all:data_obs_cat3" << endl;
  }
  outFile << "############## shape with reparametrization" << endl;
  outFile << "shapes Bkg cat0 " << TString(fileBkgName)+".root" << " w_all:CMS_bkg_8TeV_cat0" << endl;
  outFile << "shapes Bkg cat1 "<< TString(fileBkgName)+".root" << " w_all:CMS_bkg_8TeV_cat1" << endl;
  if ( _NCAT > 2 ){
  outFile << "shapes Bkg cat2 " << TString(fileBkgName)+".root" << " w_all:CMS_bkg_8TeV_cat2" << endl;
  outFile << "shapes Bkg cat3 "<< TString(fileBkgName)+".root" << " w_all:CMS_bkg_8TeV_cat3" << endl;
  }
  outFile << "# signal" << endl;
  outFile << "shapes Sig cat0 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_sig_cat0" << endl;
  outFile << "shapes Sig cat1 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_sig_cat1" << endl;
  if ( _NCAT > 2 ){
  outFile << "shapes Sig cat2 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_sig_cat2" << endl;
  outFile << "shapes Sig cat3 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_sig_cat3" << endl;
  }
  outFile << "# ggh" << endl;
  outFile << "shapes Higggh cat0 " << TString(fileHiggsNameggh)+".inputhig.root" << " w_all:CMS_hig_0_cat0" << endl;
  outFile << "shapes Higggh cat1 " << TString(fileHiggsNameggh)+".inputhig.root" << " w_all:CMS_hig_0_cat1" << endl;
  if ( _NCAT > 2 ){
  outFile << "shapes Higggh cat2 " << TString(fileHiggsNameggh)+".inputhig.root" << " w_all:CMS_hig_0_cat2" << endl;
  outFile << "shapes Higggh cat3 " << TString(fileHiggsNameggh)+".inputhig.root" << " w_all:CMS_hig_0_cat3" << endl;
  }
  outFile << "# tth" << endl;
  outFile << "shapes Higtth cat0 " << TString(fileHiggsNametth)+".inputhig.root" << " w_all:CMS_hig_1_cat0" << endl;
  outFile << "shapes Higtth cat1 " << TString(fileHiggsNametth)+".inputhig.root" << " w_all:CMS_hig_1_cat1" << endl;
  if (_NCAT > 2 ){
  outFile << "shapes Higtth cat2 " << TString(fileHiggsNametth)+".inputhig.root" << " w_all:CMS_hig_1_cat2" << endl;
  outFile << "shapes Higtth cat3 " << TString(fileHiggsNametth)+".inputhig.root" << " w_all:CMS_hig_1_cat3" << endl;
  }
  outFile << "# vbf" << endl;
  outFile << "shapes Higvbf cat0 " << TString(fileHiggsNamevbf)+".inputhig.root" << " w_all:CMS_hig_2_cat0" << endl;
  outFile << "shapes Higvbf cat1 " << TString(fileHiggsNamevbf)+".inputhig.root" << " w_all:CMS_hig_2_cat1" << endl;
  if ( _NCAT > 2 ){
  outFile << "shapes Higvbf cat2 " << TString(fileHiggsNamevbf)+".inputhig.root" << " w_all:CMS_hig_2_cat2" << endl;
  outFile << "shapes Higvbf cat3 " << TString(fileHiggsNamevbf)+".inputhig.root" << " w_all:CMS_hig_2_cat3" << endl;
  }
  outFile << "# vh" << endl;
  outFile << "shapes Higvh cat0 " << TString(fileHiggsNamevh)+".inputhig.root" << " w_all:CMS_hig_3_cat0" << endl;
  outFile << "shapes Higvh cat1 " << TString(fileHiggsNamevh)+".inputhig.root" << " w_all:CMS_hig_3_cat1" << endl;
  if ( _NCAT > 2 ){
  outFile << "shapes Higvh cat2 " << TString(fileHiggsNamevh)+".inputhig.root" << " w_all:CMS_hig_3_cat2" << endl;
  outFile << "shapes Higvh cat3 " << TString(fileHiggsNamevh)+".inputhig.root" << " w_all:CMS_hig_3_cat3" << endl;
  }
  outFile << "# bbh" << endl;
  outFile << "shapes Higbbh cat0 " << TString(fileHiggsNamebbh)+".inputhig.root" << " w_all:CMS_hig_4_cat0" << endl;
  outFile << "shapes Higbbh cat1 " << TString(fileHiggsNamebbh)+".inputhig.root" << " w_all:CMS_hig_4_cat1" << endl;
  if ( _NCAT > 2 ){
  outFile << "shapes Higbbh cat2 " << TString(fileHiggsNamebbh)+".inputhig.root" << " w_all:CMS_hig_4_cat2" << endl;
  outFile << "shapes Higbbh cat3 " << TString(fileHiggsNamebbh)+".inputhig.root" << " w_all:CMS_hig_4_cat3" << endl;
  }
  outFile << "---------------" << endl;
  /////////////////////////////////////
  if(_addHiggs) { //
    outFile << "bin cat0 cat1 ";
    if ( _NCAT > 2 ) outFile << "cat2 cat3 ";
    cout<<"here"<<endl;
    if(!_doblinding){ 
      outFile << "\nobservation "<< data[0]->sumEntries() <<" " << data[1]->sumEntries() <<" "; 
      if ( _NCAT > 2 ) outFile << data[2]->sumEntries() <<" " << data[3]->sumEntries() <<" "; 
    }
    else{
      outFile << "\nobservation -1 -1 ";
      if ( _NCAT > 2 ) outFile << "-1 -1 ";
    }
    outFile << "\n------------------------------" << endl;
    outFile << "bin cat0 cat0 cat0 cat0 cat0 cat0 cat0"
	    <<" cat1 cat1 cat1 cat1 cat1 cat1 cat1";
    if ( _NCAT > 2 ){
    outFile << " cat2 cat2 cat2 cat2 cat2 cat2 cat2"
	    <<" cat3 cat3 cat3 cat3 cat3 cat3 cat3";
    }
    outFile << "\nprocess Sig Bkg Higggh Higtth Higvbf Higvh Higbbh"
	    <<" Sig Bkg Higggh Higtth Higvbf Higvh Higbbh";
    if ( _NCAT > 2 ){
    outFile << " Sig Bkg Higggh Higtth Higvbf Higvh Higbbh"
	    << " Sig Bkg Higggh Higtth Higvbf Higvh Higbbh";
    }
    outFile << "\nprocess 0 1 2 3 4 5 6"
	    <<" 0 1 2 3 4 5 6 ";
    if ( _NCAT > 2 ){
    outFile << "0 1 2 3 4 5 6 0 1 2 3 4 5 6";
    }
    outFile << "\nrate "
	    <<" "<<sigToFit[0]->sumEntries()<<" "<<1<<" "<<higToFitggh[0]->sumEntries()<<" "<<higToFittth[0]->sumEntries()<<" "<<higToFitvbf[0]->sumEntries()<<" "<<higToFitvh[0]->sumEntries()<< " "<<higToFitbbh[0]->sumEntries()
	    <<" "<<sigToFit[1]->sumEntries()<<" "<<1<<" "<<higToFitggh[1]->sumEntries()<<" "<<higToFittth[1]->sumEntries()<<" "<<higToFitvbf[1]->sumEntries()<<" "<<higToFitvh[1]->sumEntries()<<" "<<higToFitbbh[1]->sumEntries()
	    <<" ";
    if ( _NCAT > 2 ){
    outFile <<" "<<sigToFit[2]->sumEntries()<<" "<<1<<" "<<higToFitggh[2]->sumEntries()<<" "<<higToFittth[2]->sumEntries()<<" "<<higToFitvbf[2]->sumEntries()<<" "<<higToFitvh[2]->sumEntries()<< " "<<higToFitbbh[2]->sumEntries()
	    <<" "<<sigToFit[3]->sumEntries()<<" "<<1<<" "<<higToFitggh[3]->sumEntries()<<" "<<higToFittth[3]->sumEntries()<<" "<<higToFitvbf[3]->sumEntries()<<" "<<higToFitvh[3]->sumEntries()<< " "<<higToFitbbh[3]->sumEntries();
    }
    outFile << " " << endl;
    outFile << "############## Total normalisation" << endl;
    outFile << "lumi_8TeV lnN "
	    << "1.026 - 1.026 1.026 1.026 1.026 1.026 "
	    << "1.026 - 1.026 1.026 1.026 1.026 1.026 ";
    if ( _NCAT > 2 ){
    outFile << "1.026 - 1.026 1.026 1.026 1.026 1.026 "
	    << "1.026 - 1.026 1.026 1.026 1.026 1.026 ";
    }
    outFile << " " << endl << endl;
    outFile << "############## Photon selection normalisation uncertainties " << endl;
    outFile << "DiphoTrigger lnN "
	    << "1.01 - 1.010 1.010 1.010 1.010 1.010 "
	    << "1.01 - 1.010 1.010 1.010 1.010 1.010 ";
    if ( _NCAT > 2 ){
    outFile << "1.01 - 1.010 1.010 1.010 1.010 1.010 "
	    << "1.01 - 1.010 1.010 1.010 1.010 1.010 ";
    }
    outFile << "# Trigger efficiency" << endl;
    outFile << "CMS_hgg_eff_g lnN "
	    << "1.010 - 1.010 1.010 1.010 1.010 1.010 "
	    << "1.010 - 1.010 1.010 1.010 1.010 1.010 ";
    if ( _NCAT > 2 ){
      outFile << "1.010 - 1.010 1.010 1.010 1.010 1.010 "
	      << "1.010 - 1.010 1.010 1.010 1.010 1.010 ";
    }
    outFile << "# photon selection accep." << endl;
    outFile << " " << endl;
    outFile << "############## Jet selection and phase space cuts normalisation uncertainties " << endl;
    outFile << "pTj_cut_acceptance lnN "
	    << "1.010 - 1.010 1.010 1.010 1.010 1.010 "
	    << "1.010 - 1.010 1.010 1.010 1.010 1.010 ";
    if ( _NCAT > 2 ){
    outFile << "1.010 - 1.010 1.010 1.010 1.010 1.010 "
	    << "1.010 - 1.010 1.010 1.010 1.010 1.010 ";
    }
    outFile <<"# JER and JES " << endl;
    if ( _NCAT == 2){
    outFile << "btag_eff lnN "
	    << "1.050 - 1.0508 1.050 1.050 1.050 1.050 "
	    << "0.979 - 0.979 0.979 0.979 0.979 0.979 ";
    }
    if ( _NCAT > 2 ){
    outFile << "btag_eff lnN "
	    << "1.050 - 1.050 1.050 1.050 1.050 1.050 "
	    << "0.979 - 0.979 0.979 0.979 0.979 0.979 "
	    << "1.050 - 1.050 1.050 1.050 1.050 1.050 "
	    << "0.979 - 0.979 0.979 0.979 0.979 0.979 ";
    }
    outFile <<"# b tag efficiency uncertainty" << endl;
    if (_NCAT == 2){
    outFile << "maajj_cut_acceptance lnN "
	    << "1.015 - 1.015 1.015 1.015 1.015 1.015 "
	    << "1.015 - 1.015 1.015 1.015 1.015 1.015 ";
    }
    else if (_NCAT > 2){
    outFile << "maajj_cut_acceptance lnN "
	    << "0.995 - 0.995 0.995 0.995 0.995 0.995 "
	    << "0.995 - 0.995 0.995 0.995 0.995 0.995 "
	    << "1.015 - 1.015 1.015 1.015 1.015 1.015 "
	    << "1.015 - 1.015 1.015 1.015 1.015 1.015 ";
    }
    outFile << "# uncertainty on mggjj cut acceptance" << endl;
    outFile << " " << endl << endl;
    outFile << "############## Theory uncertainties on SM Higgs production " << endl;
    outFile << "PDF lnN "
	    << " - - 0.931/1.075 0.919/1.081 0.972/1.026 0.976/1.024 0.976/1.024 "
	    << " - - 0.931/1.075 0.919/1.081 0.972/1.026 0.976/1.024 0.976/1.024 ";
    if ( _NCAT > 2 ){
    outFile << " - - 0.931/1.075 0.919/1.081 0.972/1.026 0.976/1.024 0.976/1.024 "
	    << " - - 0.931/1.075 0.919/1.081 0.972/1.026 0.976/1.024 0.976/1.024 ";
    }
    outFile << "\nQCD_scale lnN "
	    << " - - 0.922/1.072 0.907/1.038 0.998/1.002 0.980/1.020 0.980/1.020 "
	    << " - - 0.922/1.072 0.907/1.038 0.998/1.002 0.980/1.020 0.980/1.020 ";
    if ( _NCAT > 2 ){
    outFile << " - - 0.922/1.072 0.907/1.038 0.998/1.002 0.980/1.020 0.980/1.020 "
	    << " - - 0.922/1.072 0.907/1.038 0.998/1.002 0.980/1.020 0.980/1.020 ";
    }
    outFile << "\ngg_migration lnN "
	    << " - - 1.25 1.25 1.08 1.08 1.08 "
	    << " - - 1.25 1.25 1.08 1.08 1.08 ";
    if ( _NCAT > 2 ){
    outFile << " - - 1.25 1.25 1.08 1.08 1.08 "
	    << " - - 1.25 1.25 1.08 1.08 1.08 ";
    }
    outFile << "# UEPS" << endl;
    outFile << "gluonSplitting lnN "
	    << " - - 1.40 1.40 1.40 1.40 1.40 "
	    << " - - 1.40 1.40 1.40 1.40 1.40 ";
    if ( _NCAT > 2 ){
    outFile << " - - 1.40 1.40 1.40 1.40 1.40 "
	    << " - - 1.40 1.40 1.40 1.40 1.40 ";
    }
    outFile << " " << endl<<endl;
    if(useSigTheoryUnc){
      outFile << "############## Theory uncertainty on SM diHiggs production " << endl;
      outFile << "SM_diHiggs_Theory lnN "
	      << " 0.857/1.136 - - - - - - "
	      << " 0.857/1.136 - - - - - - ";
      if ( _NCAT > 2 ){
	outFile << " 0.857/1.136 - - - - - - "
		<< " 0.857/1.136 - - - - - - ";
      }
      outFile << " # from 9.96 + 1.35 - 1.42 fb " << endl << endl;
    }
    outFile << "############## Signal parametric shape uncertainties " << endl;
    if ( _sigMass > 0)
      outFile << "CMS_hgg_sig_m0_absShift param 1 0.0045 # displacement of the dipho mean error = sqrt(0.4^ 2 + 0.2^ 2) " << endl;
    else
      outFile << "CMS_hgg_sig_m0_absShift param 1 0.0054 # displacement of the dipho mean error = sqrt(0.5^ 2 + 0.2^ 2) " << endl;
    outFile << "CMS_hgg_sig_sigmaScale param 1 0.05 # optimistic estimate of resolution uncertainty " << endl;
    //
    outFile << "CMS_hbb_sig_m0_absShift param 1 0.026 # displacement of the dijet mean error " << endl;
    outFile << "CMS_hbb_sig_sigmaScale param 1 0.10 # optimistic estimate of resolution uncertainty " << endl;
    //
    outFile << "############## for mggxmjj fit - slopes" << endl;
    outFile << "CMS_bkg_8TeV_cat0_norm flatParam # Normalization uncertainty on background slope" << endl;
    outFile << "CMS_bkg_8TeV_cat1_norm flatParam # Normalization uncertainty on background slope" << endl;
    if ( _NCAT > 2 ){
    outFile << "CMS_bkg_8TeV_cat2_norm flatParam # Normalization uncertainty on background slope" << endl;
    outFile << "CMS_bkg_8TeV_cat3_norm flatParam # Normalization uncertainty on background slope" << endl;
    }

    if ( _NCAT > 2 ){
    outFile << "CMS_hgg_bkg_8TeV_slope1_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
    outFile << "CMS_hgg_bkg_8TeV_slope1_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;
    outFile << "CMS_hgg_bkg_8TeV_slope1_cat2 flatParam # Mean and absolute uncertainty on background slope" << endl;
    outFile << "CMS_hgg_bkg_8TeV_slope1_cat3 flatParam # Mean and absolute uncertainty on background slope" << endl;
    }
    else{
    outFile << "CMS_hgg_bkg_8TeV_slope2_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
    outFile << "CMS_hgg_bkg_8TeV_slope3_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
    outFile << "CMS_hgg_bkg_8TeV_slope1_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;
    }

    if ( _NCAT > 2 ){
    outFile << "CMS_hbb_bkg_8TeV_slope1_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
    outFile << "CMS_hbb_bkg_8TeV_slope1_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;
    outFile << "CMS_hbb_bkg_8TeV_slope1_cat2 flatParam # Mean and absolute uncertainty on background slope" << endl;
    outFile << "CMS_hbb_bkg_8TeV_slope1_cat3 flatParam # Mean and absolute uncertainty on background slope" << endl;
    }
    else{
    outFile << "CMS_hbb_bkg_8TeV_slope2_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
    outFile << "CMS_hbb_bkg_8TeV_slope3_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
    outFile << "CMS_hbb_bkg_8TeV_slope1_cat1 flatParam # Mean and absolute uncertainty on background slope" << endl;
    }

  } // if ncat == 2 or 4
  /////////////////////////////////////
  outFile.close();
  cout << "Write data card in: " << filename << " file" << endl;
  return;
} // close write full datacard

void bbgg2DFitter::MakeDataCardonecatnohiggs(const char* fileBaseName, const char* fileBkgName, Bool_t useSigTheoryUnc)
{
  TString cardDir = "datacards/";
  RooDataSet* data[_NCAT];
  RooDataSet* sigToFit[_NCAT];
  RooDataSet* higToFit[_NCAT];
  for (int c = 0; c < _NCAT; ++c) {
    data[c] = (RooDataSet*) _w->data(TString::Format("Data_cat%d",c));
    sigToFit[c] = (RooDataSet*) _w->data(TString::Format("Sig_cat%d",c));
    higToFit[c] = (RooDataSet*) _w->data(TString::Format("Hig_cat%d",c));
  }
  //RooRealVar* lumi = w->var("lumi");
  cout << "======== Expected Events Number =====================" << endl;
  cout << ".........Measured Data for L = " << _lumi << " pb-1 ............................" << endl;
  if(!_doblinding){ cout << "#Events data: " << _w->data("Data")->sumEntries() << endl; }
  else cout << "#Events data: -1 " << endl;
  if(!_doblinding){
     for (int c = 0; c < _NCAT; ++c) cout << TString::Format("#Events data cat%d: ",c) << data[c]->sumEntries() << endl;
  }
  else{
     for (int c = 0; c < _NCAT; ++c) cout << TString::Format("#Events data cat%d: ",c) << -1 << endl;
  }
  // cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;
  cout << ".........Expected Signal for L = " << _lumi << " pb-1 ............................" << endl;
  if(!_doblinding){ cout << "#Events Signal: " << _w->data("Data")->sumEntries() << endl; }
  else cout << "#Events Signal: -1 " << endl;
  Float_t siglikeErr[_NCAT];
  for (int c = 0; c < _NCAT; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << sigToFit[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*sigToFit[c]->sumEntries();
  }
  cout << "====================================================" << endl;
  TString filename(cardDir+TString(fileBaseName)+"onecatnohiggs.txt");
  ofstream outFile(filename);

  //outFile << "#CMS-HGG DataCard for Unbinned Limit Setting, " << lumi->getVal() << " pb-1 " << endl;
  outFile << "#Run with: combine -d hgg.mH350.0.shapes-Unbinned.txt -U -m 130 -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0 -b 3500 -i 50000 --optimizeSim=1 --tries 30" << endl;
  // outFile << "# Lumi = " << lumi->getVal() << " pb-1" << endl;
  outFile << "# Lumi = " << _lumi << " pb-1" << endl;
  outFile << "imax 1" << endl;
  outFile << "jmax 1" << endl; // number of BKG
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;

  cout<<"here"<<endl;
  outFile << "shapes data_obs cat0 " << TString(fileBkgName)+".root" << " w_all:data_obs_cat0" << endl;
  outFile << "############## shape with reparametrization" << endl;
  outFile << "shapes Bkg cat0 " << TString(fileBkgName)+".root" << " w_all:CMS_bkg_8TeV_cat0" << endl;
  outFile << "# signal" << endl;
  outFile << "shapes Sig cat0 " << TString(fileBaseName)+".inputsig.root" << " w_all:CMS_sig_cat0" << endl;


  outFile << "---------------" << endl;
  /////////////////////////////////////
  if(_addHiggs) { //
    outFile << "bin cat0 " << endl;
    if(!_doblinding){ 
    outFile << "observation "
	    << data[0]->sumEntries() << " " << endl; 
    }else{
    outFile << "observation -1 " << endl;   
    }
    outFile << "------------------------------" << endl;
    outFile << "bin cat0 cat0 " << endl;
    outFile << "process Sig Bkg " << endl;
    outFile << "process 0 1 " << endl;
    outFile << "rate "
	    << " " << sigToFit[0]->sumEntries() << " " << 1
	    << " " << endl;
    outFile << "--------------------------------" << endl;
    outFile << "lumi_8TeV lnN "
	    << "1.026 - " << endl;
    outFile << "############## jet" << endl;
    outFile << "pTj_acceptance lnN "
	    << "1.010 - "
	    <<"# JER and JES " << endl;
    outFile << "btag_eff lnN ";
    if (_sigMass==0)
      outFile << "1.050 - ";
    else
      outFile << "1.050 - ";
    outFile <<"# b tag efficiency uncertainty" << endl;
    outFile << "mggjj_acceptance lnN ";
    if (_sigMass==0)
      outFile << "0.995 - " << endl;
    else
      outFile << "1.015 - " << endl;
    outFile << "############## photon " << endl;
    outFile << "CMS_hgg_eff_g lnN "
	    << "1.010 - "
	    << "# photon selection accep." << endl;
    outFile << "DiphoTrigger lnN "
	    << "1.01 - "
	    << "# Trigger efficiency" << endl;
    if(useSigTheoryUnc){
      outFile << "SM_diHiggs_Theory lnN 0.857/1.136 - " << endl;
    }
    outFile << "# Parametric shape uncertainties, entered by hand. they act on signal " << endl;
    outFile << "CMS_hgg_sig_m0_absShift param 1 0.0054 # displacement of the dipho mean" << endl;
    outFile << "CMS_hgg_sig_sigmaScale param 1 0.05 # optimistic estimate of resolution uncertainty " << endl;
    outFile << "CMS_hbb_sig_m0_absShift param 1 0.026 # displacement of the dijet mean" << endl;
    outFile << "CMS_hbb_sig_sigmaScale param 1 0.10 # optimistic estimate of resolution uncertainty " << endl;
    outFile << "############## for mggxmjj fit - slopes" << endl;
    outFile << "CMS_bkg_8TeV_cat0_norm flatParam # Normalization uncertainty on background slope" << endl;
    if ( _NCAT > 2 ){
    outFile << "CMS_hgg_bkg_8TeV_slope1_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
    }
    else{
    outFile << "CMS_hgg_bkg_8TeV_slope2_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
    outFile << "CMS_hgg_bkg_8TeV_slope3_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
    }
    if ( _NCAT > 2 ){
    outFile << "CMS_hbb_bkg_8TeV_slope1_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
    }
    else{
    outFile << "CMS_hbb_bkg_8TeV_slope2_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
    outFile << "CMS_hbb_bkg_8TeV_slope3_cat0 flatParam # Mean and absolute uncertainty on background slope" << endl;
    }
  } // if ncat ==2
  /////////////////////////////////////

  outFile.close();
  cout << "Write data card in: " << filename << " file" << endl;
  return;
} // close write datacard one cat

void bbgg2DFitter::SetConstantParams(const RooArgSet* params)
{
  // set constant parameters for signal fit, ... NO IDEA !!!!
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }
} // close set const parameters

TStyle * bbgg2DFitter::style()
{
  TStyle *defaultStyle = new TStyle("defaultStyle","Default Style");
  defaultStyle->SetOptStat(0000);
  defaultStyle->SetOptFit(000);
  defaultStyle->SetPalette(1);
  /////// pad ////////////
  defaultStyle->SetPadBorderMode(1);
  defaultStyle->SetPadBorderSize(1);
  defaultStyle->SetPadColor(0);
  defaultStyle->SetPadTopMargin(0.05);
  defaultStyle->SetPadBottomMargin(0.13);
  defaultStyle->SetPadLeftMargin(0.13);
  defaultStyle->SetPadRightMargin(0.02);
  /////// canvas /////////
  defaultStyle->SetCanvasBorderMode(0);
  defaultStyle->SetCanvasColor(0);
  defaultStyle->SetCanvasDefH(600);
  defaultStyle->SetCanvasDefW(600);
  /////// frame //////////
  defaultStyle->SetFrameBorderMode(0);
  defaultStyle->SetFrameBorderSize(1);
  defaultStyle->SetFrameFillColor(0);
  defaultStyle->SetFrameLineColor(1);
  /////// label //////////
  defaultStyle->SetLabelOffset(0.005,"XY");
  defaultStyle->SetLabelSize(0.05,"XY");
  defaultStyle->SetLabelFont(42,"XY");
  /////// title //////////
  defaultStyle->SetTitleOffset(1.1,"X");
  defaultStyle->SetTitleSize(0.01,"X");
  defaultStyle->SetTitleOffset(1.25,"Y");
  defaultStyle->SetTitleSize(0.05,"Y");
  defaultStyle->SetTitleFont(42, "XYZ");
  /////// various ////////
  defaultStyle->SetNdivisions(505,"Y");
  defaultStyle->SetLegendBorderSize(0); // For the axis titles:

  defaultStyle->SetTitleColor(1, "XYZ");
  defaultStyle->SetTitleFont(42, "XYZ");
  defaultStyle->SetTitleSize(0.06, "XYZ");
 
  // defaultStyle->SetTitleYSize(Float_t size = 0.02);
  defaultStyle->SetTitleXOffset(0.9);
  defaultStyle->SetTitleYOffset(1.05);
  // defaultStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:
  defaultStyle->SetLabelColor(1, "XYZ");
  defaultStyle->SetLabelFont(42, "XYZ");
  defaultStyle->SetLabelOffset(0.007, "XYZ");
  defaultStyle->SetLabelSize(0.04, "XYZ");

  // For the axis:
  defaultStyle->SetAxisColor(1, "XYZ");
  defaultStyle->SetStripDecimals(kTRUE);
  defaultStyle->SetTickLength(0.03, "XYZ");
  defaultStyle->SetNdivisions(510, "XYZ");
  defaultStyle->SetPadTickX(1); // To get tick marks on the opposite side of the frame
  defaultStyle->SetPadTickY(1);
  defaultStyle->cd();
  return defaultStyle;
}

RooFitResult* bbgg2DFitter::BkgModelFit(Bool_t dobands)
{
  const Int_t ncat = _NCAT;
  std::vector<TString> catdesc;
  if ( _NCAT == 2 ){
    catdesc.push_back(" High Purity");
    catdesc.push_back(" Med. Purity");
  }
  else{
    catdesc.push_back(" #splitline{High Purity}{High m_{#gamma#gammajj}^{kin}}");
    catdesc.push_back(" #splitline{Med. Purity}{High m_{#gamma#gammajj}^{kin}}");
    catdesc.push_back(" #splitline{High Purity}{Low m_{#gamma#gammajj}^{kin}}");
    catdesc.push_back(" #splitline{Med. Purity}{Low m_{#gamma#gammajj}^{kin}}");
  }
  //******************************************//
  // Fit background with model pdfs
  //******************************************//
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[ncat];
  RooDataSet* dataplot[ncat]; // the data
  RooBernstein* mggBkg[ncat];// the polinomial of 4* order
  RooBernstein* mjjBkg[ncat];// the polinomial of 4* order
  RooPlot* plotmggBkg[ncat];
  RooPlot* plotmjjBkg[ncat];
  RooDataSet* sigToFit0[ncat];
  RooDataSet* sigToFit1[ncat];
  RooDataSet* sigToFit2[ncat];
  RooDataSet* sigToFit3[ncat];
  RooDataSet* sigToFit4[ncat];
  RooAbsPdf* mggSig[ncat];
  RooAbsPdf* mggSig0[ncat];
  RooAbsPdf* mggSig1[ncat];
  RooAbsPdf* mggSig2[ncat];
  RooAbsPdf* mggSig3[ncat];
  RooAbsPdf* mggSig4[ncat];
  RooAbsPdf* mjjSig[ncat];
  RooAbsPdf* mjjSig0[ncat];
  RooAbsPdf* mjjSig1[ncat];
  RooAbsPdf* mjjSig2[ncat];
  RooAbsPdf* mjjSig3[ncat];
  RooAbsPdf* mjjSig4[ncat];
  RooProdPdf* BkgPdf = 0;
  RooAbsPdf* mjjBkgTmpPow1 = 0;
  RooAbsPdf* mggBkgTmpPow1 = 0;
  RooExponential* mjjBkgTmpExp1 = 0;
  RooExponential* mggBkgTmpExp1 = 0;
  RooBernstein* mjjBkgTmpBer1 = 0;
  RooBernstein* mggBkgTmpBer1 = 0;
  Float_t minMggMassFit(100),maxMggMassFit(180);
  Float_t minMjjMassFit(60),maxMjjMassFit(180);
  if(_sigMass == 260) maxMggMassFit = 145;
  if(_sigMass == 270) maxMggMassFit = 155;
  // Fit data with background pdf for data limit
  RooRealVar* mgg = _w->var("mgg");
  RooRealVar* mjj = _w->var("mjj");
  RooRealVar* mtot = _w->var("mtot");
  mgg->setUnit("GeV");
  mjj->setUnit("GeV");
  mgg->setRange("BkgFitRange",minMggMassFit,maxMggMassFit);
  mjj->setRange("BkgFitRange",minMjjMassFit,maxMjjMassFit);
  RooFitResult* fitresults = new RooFitResult();
  //
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  //
  cout << "####\t Starting loop for each category" << endl;
  for (int c = 0; c < ncat; ++c) { // to each category
    data[c] = (RooDataSet*) _w->data(TString::Format("Data_cat%d",c));
    cout << "!!!!!!!!!!!!! CATEGORY: " << c <<  endl;
    ////////////////////////////////////
    // these are the parameters for the bkg polinomial
    // one slope by category - float from -10 > 10
    // we first wrap the normalization of mggBkgTmp0, mjjBkgTmp0
    _w->factory(TString::Format("bkg_8TeV_norm_cat%d[1.0,0.,100000]",c));
    if(_sigMass == 0)
    {
       if(c == 0 || c == 2)
       {
          RooFormulaVar *mgg_p1mod = new RooFormulaVar(
					     TString::Format("mgg_p1mod_cat%d",c),
					     "","@0",
					     *_w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c)));
          RooFormulaVar *mjj_p1mod = new RooFormulaVar(
					     TString::Format("mjj_p1mod_cat%d",c),
					     "","@0",
					     *_w->var(TString::Format("mjj_bkg_8TeV_slope1_cat%d",c)));

    
          mggBkgTmpExp1 = new RooExponential(
                                             TString::Format("mggBkgTmpExp1_cat%d",c), 
                                             "",
                                             *mgg,*mgg_p1mod);
          mjjBkgTmpExp1 = new RooExponential(
                                             TString::Format("mjjBkgTmpExp1_cat%d",c), 
                                             "",
                                             *mjj,*mjj_p1mod);

          BkgPdf = new RooProdPdf(TString::Format("BkgPdf_cat%d",c), "", RooArgList(*mggBkgTmpExp1, *mjjBkgTmpExp1));
       }
       else if(c == 1 || c == 3)
       { 
          RooFormulaVar *mgg_p1mod = new RooFormulaVar(
					     TString::Format("mgg_p1mod_cat%d",c),
					     "","@0",
					     *_w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c)));
          RooFormulaVar *mjj_p1mod = new RooFormulaVar(
					     TString::Format("mjj_p1mod_cat%d",c),
					     "","@0",
					     *_w->var(TString::Format("mjj_bkg_8TeV_slope1_cat%d",c)));

          mggBkgTmpPow1 = new RooGenericPdf(
					      TString::Format("mggBkgTmpPow1_cat%d",c),
					      "1./pow(@0,@1)",//"1./exp(@0*@1)",//
					      RooArgList(*mgg, *mgg_p1mod));
          mjjBkgTmpPow1 = new RooGenericPdf(
					      TString::Format("mjjBkgTmpPow1_cat%d",c),
					      "1./pow(@0,@1)",//"1./exp(@0*@1)",//
					      RooArgList(*mjj, *mjj_p1mod));

          BkgPdf = new RooProdPdf(TString::Format("BkgPdf_cat%d",c), "", RooArgList(*mggBkgTmpPow1, *mjjBkgTmpPow1));
       }
    }
    else if(_sigMass != 0)
    {
        if(c == 0)
        {
          RooFormulaVar *mgg_p0amp = new RooFormulaVar(
					     TString::Format("mgg_p0amp_cat%d",c),
					     "","@0*@0",
					     *_w->var(TString::Format("mgg_bkg_8TeV_slope2_cat%d",c)));
          RooFormulaVar *mgg_p1amp = new RooFormulaVar(
					     TString::Format("mgg_p1amp_cat%d",c),
					     "","@0*@0",
					     *_w->var(TString::Format("mgg_bkg_8TeV_slope3_cat%d",c)));
          RooFormulaVar *mjj_p0amp = new RooFormulaVar(
					     TString::Format("mjj_p0amp_cat%d",c),
					     "","@0*@0",
					     *_w->var(TString::Format("mjj_bkg_8TeV_slope2_cat%d",c)));
          RooFormulaVar *mjj_p1amp = new RooFormulaVar(
					     TString::Format("mjj_p1amp_cat%d",c),
					     "","@0*@0",
					     *_w->var(TString::Format("mjj_bkg_8TeV_slope3_cat%d",c)));
 

          mggBkgTmpBer1 = new RooBernstein(
                                                   TString::Format("mggBkgTmpBer1_cat%d",c),
                                                   "",
                                                   *mgg,RooArgList(*mgg_p0amp,*mgg_p1amp));
          mjjBkgTmpBer1 = new RooBernstein(
                                                   TString::Format("mjjBkgTmpBer1_cat%d",c),
                                                   "",
                                                   *mjj,RooArgList(*mjj_p0amp,*mjj_p1amp));

          BkgPdf = new RooProdPdf(TString::Format("BkgPdf_cat%d",c), "", RooArgList(*mggBkgTmpBer1, *mjjBkgTmpBer1));
        }
        else if(c == 1)
        {
          RooFormulaVar *mgg_p1mod = new RooFormulaVar(
					     TString::Format("mgg_p1mod_cat%d",c),
					     "","@0",
					     *_w->var(TString::Format("mgg_bkg_8TeV_slope1_cat%d",c)));
          RooFormulaVar *mjj_p1mod = new RooFormulaVar(
					     TString::Format("mjj_p1mod_cat%d",c),
					     "","@0",
					     *_w->var(TString::Format("mjj_bkg_8TeV_slope1_cat%d",c)));

          mggBkgTmpPow1 = new RooGenericPdf(
					      TString::Format("mggBkgTmpPow1_cat%d",c),
					      "1./pow(@0,@1)",//"1./exp(@0*@1)",//
					      RooArgList(*mgg, *mgg_p1mod));
          mjjBkgTmpPow1 = new RooGenericPdf(
					      TString::Format("mjjBkgTmpPow1_cat%d",c),
					      "1./pow(@0,@1)",//"1./exp(@0*@1)",//
					      RooArgList(*mjj, *mjj_p1mod));

          BkgPdf = new RooProdPdf(TString::Format("BkgPdf_cat%d",c), "", RooArgList(*mggBkgTmpPow1, *mjjBkgTmpPow1));
        } 
    }
    RooExtendPdf BkgPdfExt(TString::Format("BkgPdfExt_cat%d",c),"", *BkgPdf,*_w->var(TString::Format("bkg_8TeV_norm_cat%d",c)));
    fitresults = BkgPdfExt.fitTo(*data[c], Strategy(1),Minos(kFALSE), Range("BkgFitRange"),SumW2Error(kTRUE), Save(kTRUE));
    _w->import(BkgPdfExt);
    cout << "#### FINISHED FIT!" << endl;
    //BkgPdf.fitTo(*data[c], Strategy(1),Minos(kFALSE), Range("BkgFitRange"),SumW2Error(kTRUE), Save(kTRUE));
    //w->import(BkgPdf);

    //************************************************//
    // Plot mgg background fit results per categories
    //************************************************//
    TCanvas* ctmp = new TCanvas(TString::Format("ctmpBkgMgg_cat%d",c),"mgg Background Categories",0,0,500,500);
    Int_t nBinsMass(80);
    plotmggBkg[c] = mgg->frame(nBinsMass);
    cout<<" here 1"<<endl;
    dataplot[c] = (RooDataSet*) _w->data(TString::Format("Dataplot_cat%d",c));
    cout<<" here 1"<<endl;
    if(_sigMass == 260) plotmggBkg[c]->GetXaxis()->SetRangeUser(100.,145.);
    if(_sigMass == 270) plotmggBkg[c]->GetXaxis()->SetRangeUser(100.,155.);
    if(_doblinding) dataplot[c]->plotOn(plotmggBkg[c],Invisible());
    else dataplot[c]->plotOn(plotmggBkg[c]);
    if(_sigMass == 0)
    {
       if(c == 0 || c == 2) mggBkgTmpExp1->plotOn(plotmggBkg[c],LineColor(kBlue),Range("BkgFitRange"),NormRange("BkgFitRange"));
       else if(c == 1 || c == 3) mggBkgTmpPow1->plotOn(plotmggBkg[c],LineColor(kBlue),Range("BkgFitRange"),NormRange("BkgFitRange"));
    }
    else if(_sigMass != 0)
    {
       if(c == 0) mggBkgTmpBer1->plotOn(plotmggBkg[c],LineColor(kBlue),Range("BkgFitRange"),NormRange("BkgFitRange"));
       else if(c == 1) mggBkgTmpPow1->plotOn(plotmggBkg[c],LineColor(kBlue),Range("BkgFitRange"),NormRange("BkgFitRange"));
    }
    if(_doblinding) dataplot[c]->plotOn(plotmggBkg[c], Invisible());
    else dataplot[c]->plotOn(plotmggBkg[c]);
    cout << "!!!!!!!!!!!!!!!!!" << endl;
    cout << "!!!!!!!!!!!!!!!!!" << endl; // now we fit the gaussian on signal
    //plotmggBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
    if(c==0||c==2)plotmggBkg[c]->SetMinimum(0.001); // no error bar in bins with zero events
    if(c==1||c==3)plotmggBkg[c]->SetMinimum(0.001); // no error bar in bins with zero events
    plotmggBkg[c]->Draw();
    //plotmggBkg[c]->SetTitle("CMS preliminary 19.7/fb");
    //plotmggBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
    plotmggBkg[c]->SetMaximum(1.40*plotmggBkg[c]->GetMaximum());
    plotmggBkg[c]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
    //double test = sigToFit[c]->sumEntries();
    //cout<<"number of events on dataset "<<test<<endl;
    TPaveText *pt = new TPaveText(0.1,0.94,0.9,0.99, "brNDC");
    // pt->SetName("title");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetTextSize(0.035);
    pt->AddText("            CMS Preliminary                     L = 19.7 fb^{-1}    #sqrt{s} = 8 TeV   ");
    pt->Draw();
    TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
    TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
    if (dobands) {
      RooAbsPdf *cpdf;
      cpdf = mggBkgTmpExp1;
      if(_sigMass == 0 && (c == 0 || c == 2)) cpdf = mggBkgTmpExp1;
      if(_sigMass == 0 && (c == 1 || c == 3)) cpdf = mggBkgTmpPow1;
      if(_sigMass != 0 && c == 0) cpdf = mggBkgTmpBer1;
      if(_sigMass != 0 && c == 1) cpdf = mggBkgTmpPow1;
//      onesigma = new TGraphAsymmErrors();
//      twosigma = new TGraphAsymmErrors();
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotmggBkg[c]->getObject(1));
      for (int i=1; i<(plotmggBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
        double lowedge = plotmggBkg[c]->GetXaxis()->GetBinLowEdge(i);
        double upedge = plotmggBkg[c]->GetXaxis()->GetBinUpEdge(i);
        double center = plotmggBkg[c]->GetXaxis()->GetBinCenter(i);
        double nombkg = nomcurve->interpolate(center);
        nlim->setVal(nombkg);
        mgg->setRange("errRange",lowedge,upedge);
        RooAbsPdf *epdf = 0;
        epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
        RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
        RooMinimizer minim(*nll);
        minim.setStrategy(0);
        double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
        double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
        minim.migrad();
        minim.minos(*nlim);
        // printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
        onesigma->SetPoint(i-1,center,nombkg);
        onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
        minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2));
        // the 0.5 is because qmu is -2*NLL
        // eventually if cl = 0.95 this is the usual 1.92!
        minim.migrad();
        minim.minos(*nlim);
        twosigma->SetPoint(i-1,center,nombkg);
        twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
        delete nll;
        delete epdf;
      } // close for bin
      mgg->setRange("errRange",minMggMassFit,maxMggMassFit);
      twosigma->SetLineColor(kYellow);
      twosigma->SetFillColor(kYellow);
      twosigma->SetMarkerColor(kYellow);
      twosigma->Draw("L3 SAME");
      onesigma->SetLineColor(kGreen);
      onesigma->SetFillColor(kGreen);
      onesigma->SetMarkerColor(kGreen);
      onesigma->Draw("L3 SAME");
      plotmggBkg[c]->Draw("SAME");
    } else plotmggBkg[c]->Draw("SAME"); // close dobands
    //plotmggBkg[c]->getObject(1)->Draw("SAME");
    //plotmggBkg[c]->getObject(2)->Draw("P SAME");
    ////////////////////////////////////////////////////////// plot higgs
    sigToFit0[c] = (RooDataSet*) _w->data(TString::Format("Hig_0_cat%d",c));
    double norm0; norm0 = 1.0*sigToFit0[c]->sumEntries(); //
    //norm0 = 0.0000001;
    mggSig0[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggHig_0_cat%d",c));
    // we are not constructing signal pdf, this is constructed on sig to fit function...
    mggSig0[c] ->plotOn(
			plotmggBkg[c],
			Normalization(norm0,RooAbsPdf::NumEvent),
			DrawOption("F"),
			LineColor(kRed),FillStyle(1001),FillColor(19));
    mggSig0[c]->plotOn(
		       plotmggBkg[c],
		       Normalization(norm0,RooAbsPdf::NumEvent),LineColor(kRed),LineStyle(1));
    //
    sigToFit1[c] = (RooDataSet*) _w->data(TString::Format("Hig_1_cat%d",c));
    double norm1 = 1.0*sigToFit1[c]->sumEntries(); //
    mggSig1[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggHig_1_cat%d",c));
    // we are not constructing signal pdf, this is constructed on sig to fit function...
    mggSig1[c] ->plotOn(
			plotmggBkg[c],
			Normalization(norm1,RooAbsPdf::NumEvent),
			DrawOption("F"),
			LineColor(kGreen),FillStyle(1001),FillColor(19));
    mggSig1[c]->plotOn(
		       plotmggBkg[c],
		       Normalization(norm1,RooAbsPdf::NumEvent),LineColor(kGreen),LineStyle(1));
    //
    sigToFit2[c] = (RooDataSet*) _w->data(TString::Format("Hig_2_cat%d",c));
    double norm2;
    //if(sigToFit2[c]->sumEntries()>0)
    norm2 = 1.0*sigToFit2[c]->sumEntries(); //else
    //norm2 = 0.0000000000001; //
    mggSig2[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggHig_2_cat%d",c));
    // we are not constructing signal pdf, this is constructed on sig to fit function...
    mggSig2[c] ->plotOn(
			plotmggBkg[c],
			Normalization(norm2,RooAbsPdf::NumEvent),
			DrawOption("F"),
			LineColor(kMagenta),FillStyle(1001),FillColor(19));
    mggSig2[c]->plotOn(
		       plotmggBkg[c],
		       Normalization(norm2,RooAbsPdf::NumEvent),LineColor(kMagenta),LineStyle(1));
    //
    sigToFit3[c] = (RooDataSet*) _w->data(TString::Format("Hig_3_cat%d",c));
    double norm3 = 1.0*sigToFit3[c]->sumEntries(); //
    mggSig3[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggHig_3_cat%d",c));
    // we are not constructing signal pdf, this is constructed on sig to fit function...
    mggSig3[c] ->plotOn(
			plotmggBkg[c],
			Normalization(norm3,RooAbsPdf::NumEvent),
			DrawOption("F"),
			LineColor(kCyan),FillStyle(1001),FillColor(19));
    mggSig3[c]->plotOn(
		       plotmggBkg[c],
		       Normalization(norm3,RooAbsPdf::NumEvent),LineColor(kCyan),LineStyle(1));

    sigToFit4[c] = (RooDataSet*) _w->data(TString::Format("Hig_4_cat%d",c));
    double norm4 = 1.0*sigToFit4[c]->sumEntries(); //
    mggSig4[c] = (RooAbsPdf*) _w->pdf(TString::Format("mggHig_4_cat%d",c));
    // we are not constructing signal pdf, this is constructed on sig to fit function...
    mggSig4[c] ->plotOn(
			plotmggBkg[c],
			Normalization(norm4,RooAbsPdf::NumEvent),
			DrawOption("F"),
			LineColor(kBlue),FillStyle(1001),FillColor(19));
    mggSig4[c]->plotOn(
		       plotmggBkg[c],
		       Normalization(norm4,RooAbsPdf::NumEvent),LineColor(kBlue),LineStyle(1));

    //////////////////////////////////////////////////////////
    plotmggBkg[c]->Draw("SAME");
    if(c==0||c==2)plotmggBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
    if(c==1||c==3)plotmggBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
    if(c==0||c==2)plotmggBkg[c]->SetMaximum(5.3); // no error bar in bins with zero events
    if(c==1||c==3)plotmggBkg[c]->SetMaximum(20); // no error bar in bins with zero events
    // plotmggBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
    //plotmggBkg[c]->SetLogy(0);
    cout << "!!!!!!!!!!!!!!!!!" << endl;
    TLegend *legmc = new TLegend(0.40,0.72,0.62,0.9);
    TLegend *legmcH = new TLegend(0.66,0.72,0.94,0.9);
    if(_doblinding) legmc->AddEntry(plotmggBkg[c]->getObject(2),"Data ","");
    else  legmc->AddEntry(plotmggBkg[c]->getObject(2),"Data ","LPE");
    legmc->AddEntry(plotmggBkg[c]->getObject(1),"Bkg Fit","L");
    if(dobands)legmc->AddEntry(onesigma,"Fit #pm1 #sigma","F");
    if(dobands)legmc->AddEntry(twosigma,"Fit #pm2 #sigma ","F"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(3),"ggH ","LPE"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(5),"ttH ","LPE"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(7),"VBF ","LPE"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(9),"VH ","LPE"); // not...
    legmcH->AddEntry(plotmggBkg[c]->getObject(11),"bbH ","LPE"); // not...
    if(_sigMass==0)
      legmc->SetHeader(" Nonresonant HH");
    else
      legmc->SetHeader(TString::Format(" m_{X} = %d GeV",_sigMass));
    legmcH->SetHeader(" Higgs");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmcH->SetBorderSize(0);
    legmcH->SetFillStyle(0);
    legmc->Draw();
    legmcH->Draw();
    TLatex *lat2 = new TLatex(minMggMassFit+1.5,0.85*plotmggBkg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    //
    ctmp->SaveAs(TString::Format("databkgoversigMgg_cat%d.pdf",c));
    ctmp->SaveAs(TString::Format("databkgoversigMgg_cat%d.png",c));

    if(c==0||c==2)plotmggBkg[c]->SetMaximum(100); // no error bar in bins with zero events
    if(c==1||c==3)plotmggBkg[c]->SetMaximum(1000); // no error bar in bins with zero events
    ctmp->SetLogy(1);
    ctmp->SaveAs(TString::Format("databkgoversigMgg_cat%d_log.pdf",c));
    ctmp->SaveAs(TString::Format("databkgoversigMgg_cat%d_log.png",c));
    // ctmp->SaveAs(TString::Format("databkgoversigMgg_cat%d.C",c));

    //************************************************//
    // Plot mjj background fit results per categories
    //************************************************//
    ctmp = new TCanvas(TString::Format("ctmpBkgMjj_cat%d",c),"mjj Background Categories",0,0,500,500);
    nBinsMass = 60;
    plotmjjBkg[c] = mjj->frame(nBinsMass);
    cout<<" here 1"<<endl;
    dataplot[c] = (RooDataSet*) _w->data(TString::Format("Dataplot_cat%d",c));
    cout<<" here 1"<<endl;
    if(_doblinding) dataplot[c]->plotOn(plotmjjBkg[c],Invisible());
    else dataplot[c]->plotOn(plotmjjBkg[c]);
    if(_sigMass == 0)
    {
       if(c == 0 || c == 2)mjjBkgTmpExp1->plotOn(plotmjjBkg[c],LineColor(kBlue),Range("BkgFitRange"),NormRange("BkgFitRange"));
       else if(c == 1 || c == 3) mjjBkgTmpPow1->plotOn(plotmjjBkg[c],LineColor(kBlue),Range("BkgFitRange"),NormRange("BkgFitRange"));
    }
    else if(_sigMass != 0)
    {
       if(c == 0) mjjBkgTmpBer1->plotOn(plotmjjBkg[c],LineColor(kBlue),Range("BkgFitRange"),NormRange("BkgFitRange"));
       else if(c == 1) mjjBkgTmpPow1->plotOn(plotmjjBkg[c],LineColor(kBlue),Range("BkgFitRange"),NormRange("BkgFitRange"));
    }
    if(_doblinding) dataplot[c]->plotOn(plotmjjBkg[c],Invisible());
    else dataplot[c]->plotOn(plotmjjBkg[c]);

    cout << "!!!!!!!!!!!!!!!!!" << endl;
    cout << "!!!!!!!!!!!!!!!!!" << endl; // now we fit the gaussian on signal
    //plotmjjBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
    if(c==0||c==2)plotmjjBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
    if(c==1||c==3)plotmjjBkg[c]->SetMinimum(0.001); // no error bar in bins with zero events
    plotmjjBkg[c]->Draw();
    //plotmjjBkg[c]->SetTitle("CMS preliminary 19.7/fb");
    //plotmjjBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
    plotmjjBkg[c]->SetMaximum(1.40*plotmjjBkg[c]->GetMaximum());
    plotmjjBkg[c]->GetXaxis()->SetTitle("m_{jj} (GeV)");
    //double test = sigToFit[c]->sumEntries();
    //cout<<"number of events on dataset "<<test<<endl;
    pt = new TPaveText(0.1,0.94,0.9,0.99, "brNDC");
    // pt->SetName("title");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetTextSize(0.035);
    pt->AddText("            CMS Preliminary                     L = 19.7 fb^{-1}    #sqrt{s} = 8 TeV   ");
    pt->Draw();
    if (dobands) {
      RooAbsPdf *cpdf;
      cpdf = mjjBkgTmpExp1;
      if(_sigMass == 0 && (c == 0 || c == 2)) cpdf = mjjBkgTmpExp1;
      if(_sigMass == 0 && (c == 1 || c == 3)) cpdf = mjjBkgTmpPow1;
      if(_sigMass != 0 && c == 0) cpdf = mjjBkgTmpBer1;
      if(_sigMass != 0 && c == 1) cpdf = mjjBkgTmpPow1;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotmjjBkg[c]->getObject(1));
      for (int i=1; i<(plotmjjBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
        double lowedge = plotmjjBkg[c]->GetXaxis()->GetBinLowEdge(i);
        double upedge = plotmjjBkg[c]->GetXaxis()->GetBinUpEdge(i);
        double center = plotmjjBkg[c]->GetXaxis()->GetBinCenter(i);
        double nombkg = nomcurve->interpolate(center);
        nlim->setVal(nombkg);
        mjj->setRange("errRange",lowedge,upedge);
        RooAbsPdf *epdf = 0;
        epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
        RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
        RooMinimizer minim(*nll);
        minim.setStrategy(0);
        double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
        double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
        minim.migrad();
        minim.minos(*nlim);
        // printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
        onesigma->SetPoint(i-1,center,nombkg);
        onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
        minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2));
        // the 0.5 is because qmu is -2*NLL
        // eventually if cl = 0.95 this is the usual 1.92!
        minim.migrad();
        minim.minos(*nlim);
        twosigma->SetPoint(i-1,center,nombkg);
        twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
        delete nll;
        delete epdf;
      } // close for bin
      mjj->setRange("errRange",minMjjMassFit,maxMjjMassFit);
      twosigma->SetLineColor(kYellow);
      twosigma->SetFillColor(kYellow);
      twosigma->SetMarkerColor(kYellow);
      twosigma->Draw("L3 SAME");
      onesigma->SetLineColor(kGreen);
      onesigma->SetFillColor(kGreen);
      onesigma->SetMarkerColor(kGreen);
      onesigma->Draw("L3 SAME");
      plotmjjBkg[c]->Draw("SAME");
    } else plotmjjBkg[c]->Draw("SAME"); // close dobands
    //plotmjjBkg[c]->getObject(1)->Draw("SAME");
    //plotmjjBkg[c]->getObject(2)->Draw("P SAME");
    ////////////////////////////////////////////////////////// plot higgs
    sigToFit0[c] = (RooDataSet*) _w->data(TString::Format("Hig_0_cat%d",c));
    norm0 = 1.0*sigToFit0[c]->sumEntries(); //
    //norm0 = 0.0000001;
    mjjSig0[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_0_cat%d",c));
    // we are not constructing signal pdf, this is constructed on sig to fit function...
    mjjSig0[c] ->plotOn(
			plotmjjBkg[c],
			Normalization(norm0,RooAbsPdf::NumEvent),
			DrawOption("F"),
			LineColor(kRed),FillStyle(1001),FillColor(19));
    mjjSig0[c]->plotOn(
		       plotmjjBkg[c],
		       Normalization(norm0,RooAbsPdf::NumEvent),LineColor(kRed),LineStyle(1));
    //
    sigToFit1[c] = (RooDataSet*) _w->data(TString::Format("Hig_1_cat%d",c));
    norm1 = 1.0*sigToFit1[c]->sumEntries(); //
    mjjSig1[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_1_cat%d",c));
    // we are not constructing signal pdf, this is constructed on sig to fit function...
    mjjSig1[c] ->plotOn(
			plotmjjBkg[c],
			Normalization(norm1,RooAbsPdf::NumEvent),
			DrawOption("F"),
			LineColor(kGreen),FillStyle(1001),FillColor(19));
    mjjSig1[c]->plotOn(
		       plotmjjBkg[c],
		       Normalization(norm1,RooAbsPdf::NumEvent),LineColor(kGreen),LineStyle(1));
    //
    sigToFit2[c] = (RooDataSet*) _w->data(TString::Format("Hig_2_cat%d",c));
    //if(sigToFit2[c]->sumEntries()>0)
    norm2 = 1.0*sigToFit2[c]->sumEntries(); //else
    //norm2 = 0.0000000000001; //
    mjjSig2[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_2_cat%d",c));
    // we are not constructing signal pdf, this is constructed on sig to fit function...
    mjjSig2[c] ->plotOn(
			plotmjjBkg[c],
			Normalization(norm2,RooAbsPdf::NumEvent),
			DrawOption("F"),
			LineColor(kMagenta),FillStyle(1001),FillColor(19));
    mjjSig2[c]->plotOn(
		       plotmjjBkg[c],
		       Normalization(norm2,RooAbsPdf::NumEvent),LineColor(kMagenta),LineStyle(1));
    //
    sigToFit3[c] = (RooDataSet*) _w->data(TString::Format("Hig_3_cat%d",c));
    norm3 = 1.0*sigToFit3[c]->sumEntries(); //
    mjjSig3[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_3_cat%d",c));
    // we are not constructing signal pdf, this is constructed on sig to fit function...
    mjjSig3[c] ->plotOn(
			plotmjjBkg[c],
			Normalization(norm3,RooAbsPdf::NumEvent),
			DrawOption("F"),
			LineColor(kCyan),FillStyle(1001),FillColor(19));
    mjjSig3[c]->plotOn(
		       plotmjjBkg[c],
		       Normalization(norm3,RooAbsPdf::NumEvent),LineColor(kCyan),LineStyle(1));

    sigToFit4[c] = (RooDataSet*) _w->data(TString::Format("Hig_4_cat%d",c));
    norm4 = 1.0*sigToFit4[c]->sumEntries(); //
    mjjSig4[c] = (RooAbsPdf*) _w->pdf(TString::Format("mjjHig_4_cat%d",c));
    // we are not constructing signal pdf, this is constructed on sig to fit function...
    mjjSig4[c] ->plotOn(
			plotmjjBkg[c],
			Normalization(norm4,RooAbsPdf::NumEvent),
			DrawOption("F"),
			LineColor(kBlue),FillStyle(1001),FillColor(19));
    mjjSig4[c]->plotOn(
		       plotmjjBkg[c],
		       Normalization(norm4,RooAbsPdf::NumEvent),LineColor(kBlue),LineStyle(1));

    //////////////////////////////////////////////////////////
    plotmjjBkg[c]->Draw("SAME");
    if(c==0||c==2)plotmjjBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
    if(c==1||c==3)plotmjjBkg[c]->SetMinimum(0.01); // no error bar in bins with zero events
    if(c==0||c==2)plotmjjBkg[c]->SetMaximum(5.3); // no error bar in bins with zero events
    if(c==1||c==3)plotmjjBkg[c]->SetMaximum(20); // no error bar in bins with zero events
    // plotmjjBkg[c]->SetMinimum(0.005); // no error bar in bins with zero events
    //plotmjjBkg[c]->SetLogy(0);
    cout << "!!!!!!!!!!!!!!!!!" << endl;
    legmc = new TLegend(0.40,0.72,0.62,0.9);
    legmcH = new TLegend(0.66,0.72,0.94,0.9);
    if(_doblinding) legmc->AddEntry(plotmjjBkg[c]->getObject(2),"Data ","");
    else legmc->AddEntry(plotmjjBkg[c]->getObject(2),"Data ","LPE");
    legmc->AddEntry(plotmjjBkg[c]->getObject(1),"Fit","L");
    if(dobands)legmc->AddEntry(twosigma,"two sigma ","F"); // not...
    if(dobands)legmc->AddEntry(onesigma,"one sigma","F");
    legmcH->AddEntry(plotmjjBkg[c]->getObject(3),"ggH ","LPE"); // not...
    legmcH->AddEntry(plotmjjBkg[c]->getObject(5),"ttH ","LPE"); // not...
    legmcH->AddEntry(plotmjjBkg[c]->getObject(7),"VBF ","LPE"); // not...
    legmcH->AddEntry(plotmjjBkg[c]->getObject(9),"VH ","LPE"); // not...
    legmcH->AddEntry(plotmjjBkg[c]->getObject(11),"bbH ","LPE"); // not...
    if(_sigMass==0)
      legmc->SetHeader(" Nonresonant HH");
    else
      legmc->SetHeader(TString::Format(" m_{X} = %d GeV",_sigMass));
    legmcH->SetHeader(" Higgs");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmcH->SetBorderSize(0);
    legmcH->SetFillStyle(0);
    legmc->Draw();
    legmcH->Draw();
    lat2 = new TLatex(minMjjMassFit+1.5,0.85*plotmjjBkg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    //
    ctmp->SaveAs(TString::Format("databkgoversigMjj_cat%d.pdf",c));
    ctmp->SaveAs(TString::Format("databkgoversigMjj_cat%d.png",c));

    if(c==0||c==2)plotmjjBkg[c]->SetMaximum(100); // no error bar in bins with zero events
    if(c==1||c==3)plotmjjBkg[c]->SetMaximum(1000); // no error bar in bins with zero events
    ctmp->SetLogy(1);
    ctmp->SaveAs(TString::Format("databkgoversigMjj_cat%d_log.pdf",c));
    ctmp->SaveAs(TString::Format("databkgoversigMjj_cat%d_log.png",c));
    // ctmp->SaveAs(TString::Format("databkgoversigMjj_cat%d.C",c));

  } // close to each category

  /*RooGenericPdf *mggBkgAll = new RooGenericPdf("mggBkgAll", "1./pow(@0,@1)", RooArgList(*mgg,*w->var("mgg_bkg_8TeV_slope1")));
  RooGenericPdf *mjjBkgAll = new RooGenericPdf("mjjBkgAll", "1./pow(@0,@1)", RooArgList(*mjj,*w->var("mjj_bkg_8TeV_slope1")));

  RooProdPdf BkgPdfAll("BkgPdfAll", "Background Pdf", *mggBkgAll, *mjjBkgAll);
  RooFitResult* fitresults = BkgPdfAll.fitTo( // save results to workspace
					     *w->data("Data"),
					     Range("BkgFitRange"),
					     SumW2Error(kTRUE), Save(kTRUE));
  fitresults->Print();
  w->import(BkgPdfAll);*/
  return fitresults;
} // close berestein 3
