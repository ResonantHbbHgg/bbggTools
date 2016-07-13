#define bbggLTMaker_cxx
// The class definition in bbggLTMaker.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("bbggLTMaker.C")
// root> T->Process("bbggLTMaker.C","some options")
// root> T->Process("bbggLTMaker.C+")
//

#include "flashgg/bbggTools/interface/bbggLTMaker.h"
#include <TH2.h>
#include <TStyle.h>
#include <TSelector.h>


void bbggLTMaker::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   TString fileName = "";
   TString cat = "";
   if(!option.Contains(";")) {
       std::cout << "[bbggLTMaker::Begin] Please, make sure your input is in the form \"<fileName>;<Category>\", where <Category> : 0 (no b-tag), 1 (loose b-tag), 2 (tight b-tag)." << endl;
       TSelector::Abort("Option input doesn't contain ;");
   } else {
       TObjArray *options = option.Tokenize(";");
//       std::cout << "[bbggLTMaker::Begin] Number of options passed: " << options->GetEntries() << std::endl;
       if (options->GetEntries() < 5) {
	   std::cout << "[bbggLTMaker::Begin] Error in options parsing from LimitTreeMaker!" << std::endl;
           TSelector::Abort("Array of options divided by ; contains less than 5 options");
       } else {
           fileName = ((TObjString *)(options->At(0)))->String();
	   mtotMin = (((TObjString *)(options->At(1)))->String()).Atof();
	   mtotMax = (((TObjString *)(options->At(2)))->String()).Atof();
	   normalization = (((TObjString *)(options->At(3)))->String()).Atof();
	   photonCR = (((TObjString *)(options->At(4)))->String()).Atoi();
	   doKinFit = (((TObjString *)(options->At(5)))->String()).Atoi();
	   doMX = (((TObjString *)(options->At(6)))->String()).Atoi();
       if(doKinFit == 1 && doMX == 1){
            std::cout << "[bbggLTMaker::Begin] You need to choose either MX or KinFit, not both!" << std::endl;
            TSelector::Abort("Wrong input!");
       }
	   std::cout << "[bbggLTMaker::Begin] Using normalization factor = " << normalization << std::endl;
	   std::cout << "[bbggLTMaker::Begin] Doing photon control region? " << photonCR << std::endl;
	   std::cout << "[bbggLTMaker::Begin] Cutting on Kinematic Fitted M(4body) ? " << doKinFit << std::endl;
	   std::cout << "[bbggLTMaker::Begin] Cutting on MX ? " << doMX << std::endl;
       }
   }
   if(!fileName.Contains(".root")) {
       std::cout << "[bbggLTMaker::Begin] Your input file name is not a root file!" << std::endl;
       TSelector::Abort("Wrong input!");
   } else {
       o_evt = 0;
       o_run = 0;

       o_weight = 0;
       o_bbMass = 0;
       o_ggMass = 0;
       o_bbggMass = 0;
       o_category = -1;
       o_normalization = normalization;
       outFile = new TFile(fileName, "RECREATE");
       outTree = new TTree("TCVARS", "Limit tree for HH->bbgg analyses");
       outTree->Branch("cut_based_ct", &o_category, "o_category/I"); //0: 2btag, 1: 1btag
       outTree->Branch("evWeight", &o_weight, "o_weight/D");
       outTree->Branch("mjj", &o_bbMass, "o_bbMass/D");
       outTree->Branch("mgg", &o_ggMass, "o_ggMass/D");
       outTree->Branch("mtot", &o_bbggMass, "o_bbggMass/D"); //

       outTree->Branch("evt", &o_evt, "o_evt/l");
       outTree->Branch("run", &o_run, "o_run/i");
   }
/*
   std::cout << "[bbggLTMaker::Begin] Category chosen: " << cat << std::endl;
 
  if ( (!cat.Contains("0")) && (!cat.Contains("1")) && (!cat.Contains("2"))) {
       std::cout << "[bbggLTMaker::Begin] Please, make sure your input is in the form \"<fileName>;<Category>\", where <Category> : 0 (no b-tag), 1 (loose b-tag), 2 (tight b-tag)." << endl;
       TSelector::Abort("Category is neither 0, 1 or 2");
   } else {
       if(cat.Contains("0")) o_category = 0;
       if(cat.Contains("1")) o_category = 1;
       if(cat.Contains("2")) o_category = 2;
   }*/
   

}

void bbggLTMaker::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t bbggLTMaker::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either bbggLTMaker::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   bbggLTMaker::GetEntry(entry);
   if( entry%1000 == 0 ) std::cout << "[bbggLTMaker::Process] Reading entry #" << entry << endl;   
   
   o_evt = event;
   o_evt = run;

   o_weight = genTotalWeight*normalization;
   o_bbMass = dijetCandidate->M();
   o_ggMass = diphotonCandidate->M();
   o_bbggMass = diHiggsCandidate->M();
   if(doKinFit)
        o_bbggMass = diHiggsCandidate_KF->M();
   if(doMX)
        o_bbggMass = diHiggsCandidate->M() - dijetCandidate->M() + 125.;

  
//   if( dijetCandidate->Pt() < 50 )
//	return kTRUE;

   //mtot cut
   if(o_bbggMass < mtotMin || o_bbggMass > mtotMax)
	return kTRUE;
   if(photonCR == 1 && isPhotonCR == 0)
	return kTRUE;
   if(photonCR == 0 && isPhotonCR == 1)
	return kTRUE;
   
//   double sumbtag = leadingJet_bDis + subleadingJet_bDis;
//   double upper = 1.83;
//   double lower = 1.11;
//   if ( sumbtag > upper ) o_category = 0;
//   if (sumbtag > lower && sumbtag < upper ) o_category = 1;
//   if (sumbtag < lower ) o_category = -1;
   if ( leadingJet_bDis > 0.8 && subleadingJet_bDis > 0.8 ) {o_category = 0;}
   else if ( leadingJet_bDis > 0.8 && subleadingJet_bDis < 0.8 ) { o_category = 1; }
   else if ( leadingJet_bDis < 0.8 && subleadingJet_bDis > 0.8 ) { o_category = 1; }
   else if ( leadingJet_bDis < 0.8 && subleadingJet_bDis < 0.8 ) { o_category = -1; }
   outTree->Fill();

   return kTRUE;
}

void bbggLTMaker::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void bbggLTMaker::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
    outTree->Write();
    outFile->Close();

}
