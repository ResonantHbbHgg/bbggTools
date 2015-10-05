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
       if (options->GetEntries() < 2) {
           std::cout << "[bbggLTMaker::Begin] Please, make sure your input is in the form \"<fileName>;<Category>\", where <Category> : 0 (no b-tag), 1 (loose b-tag), 2 (tight b-tag)." << endl;
           TSelector::Abort("Array of options divided by ; contains less than 2 options");
       } else {
           fileName = ((TObjString *)(options->At(0)))->String();
           cat = ((TObjString *)(options->At(1)))->String();
       }
   }
   if(!fileName.Contains(".root")) {
       std::cout << "[bbggLTMaker::Begin] Your input file name is not a root file!" << std::endl;
       TSelector::Abort("Wrong input!");
   } else {
       o_weight = 0;
       o_bbMass = 0;
       o_ggMass = 0;
       o_bbggMass = 0;
       o_category = -1;
       outFile = new TFile(fileName, "RECREATE");
       outTree = new TTree("bbggLimitTree", "Limit tree for HH->bbgg analyses");
       outTree->Branch("category", &o_category, "o_category/I");
       outTree->Branch("weight", &o_weight, "o_weight/D");
       outTree->Branch("bbMass", &o_bbMass, "o_bbMass/D");
       outTree->Branch("ggMass", &o_ggMass, "o_ggMass/D");
       outTree->Branch("bbggMass", &o_bbggMass, "o_bbggMass/D");
   }
   std::cout << "[bbggLTMaker::Begin] Category chosen: " << cat << std::endl;
   if ( (!cat.Contains("0")) && (!cat.Contains("1")) && (!cat.Contains("2"))) {
       std::cout << "[bbggLTMaker::Begin] Please, make sure your input is in the form \"<fileName>;<Category>\", where <Category> : 0 (no b-tag), 1 (loose b-tag), 2 (tight b-tag)." << endl;
       TSelector::Abort("Category is neither 0, 1 or 2");
   } else {
       if(cat.Contains("0")) o_category = 0;
       if(cat.Contains("1")) o_category = 1;
       if(cat.Contains("2")) o_category = 2;
   }
   

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
   
   o_weight = genTotalWeight;
   o_bbMass = dijetCandidate->M();
   o_ggMass = diphotonCandidate->M();
   o_bbggMass = diHiggsCandidate->M();
   double sumbtag = leadingJet_bDis + subleadingJet_bDis;
   
   if( o_category == 0 ) outTree->Fill();
   if( o_category == 1 ) {
       if (sumbtag > 0.82 && sumbtag < 1.64 ) outTree->Fill();
   }
   if( o_category == 2 ) {
       if ( sumbtag > 1.64 ) outTree->Fill();
   }

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
