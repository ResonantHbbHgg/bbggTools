//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 15 10:49:49 2015 by ROOT version 6.02/05
// from TTree bbggSelectionTree/Flat tree for HH->bbgg analyses (after pre selection)
// found on file: bbggLTMaker_DoubleEG_37.root
//////////////////////////////////////////////////////////

#ifndef bbggLTMaker_h
#define bbggLTMaker_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <cmath>
#include <Math/LorentzVector.h>
#include <algorithm>
#include <string>
#include <utility>

using namespace std;

class bbggLTMaker : public TSelector {
public :
//   Output file and tree
   TTree *outTree;
   TFile *outFile;
   Int_t           o_category;
   Double_t        o_normalization;
   Double_t        o_weight;
   Double_t        o_bbMass;
   Double_t        o_ggMass;
   Double_t        o_bbggMass;
   double mtotMin;
   double mtotMax;
   double normalization;
   int photonCR;

//   Input tree
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<double>  *genWeights;
   Double_t        genTotalWeight;
   LorentzVector *leadingPhoton;
   vector<int>     *leadingPhotonID;
   vector<int>     *leadingPhotonISO;
   Int_t           leadingPhotonEVeto;
   LorentzVector *subleadingPhoton;
   vector<int>     *subleadingPhotonID;
   vector<int>     *subleadingPhotonISO;
   Int_t           subleadingPhotonEVeto;
   LorentzVector *diphotonCandidate;
   Int_t           nPromptInDiPhoton;
   LorentzVector *leadingJet;
   Float_t         leadingJet_bDis;
   LorentzVector *subleadingJet;
   Float_t         subleadingJet_bDis;
   LorentzVector *dijetCandidate;
   LorentzVector *diHiggsCandidate;
   Int_t	isSignal;
   Int_t	isPhotonCR;

   // List of branches
   TBranch        *b_genWeights;   //!
   TBranch        *b_genTotalWeight;   //!
   TBranch        *b_leadingPhoton;   //!
   TBranch        *b_leadingPhotonID;   //!
   TBranch        *b_leadingPhotonISO;   //!
   TBranch        *b_leadingPhotonEVeto;   //!
   TBranch        *b_subleadingPhoton;   //!
   TBranch        *b_subleadingPhotonID;   //!
   TBranch        *b_subleadingPhotonISO;   //!
   TBranch        *b_subleadingPhotonEVeto;   //!
   TBranch        *b_diphotonCandidate;   //!
   TBranch        *b_nPromptInDiPhoton;   //!
   TBranch        *b_leadingJet;   //!
   TBranch        *b_leadingJet_bDis;   //!
   TBranch        *b_subleadingJet;   //!
   TBranch        *b_subleadingJet_bDis;   //!
   TBranch        *b_dijetCandidate;   //!
   TBranch        *b_diHiggsCandidate;   //!
   TBranch	  *b_isSignal;	//!
   TBranch	  *b_isPhotonCR;  //!


   bbggLTMaker(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~bbggLTMaker() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(bbggLTMaker,0);
};

#endif

#ifdef bbggLTMaker_cxx
void bbggLTMaker::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   leadingPhoton = 0;
   subleadingPhoton = 0;
   diphotonCandidate = 0;
   leadingJet = 0;
   subleadingJet = 0;
   dijetCandidate = 0;
   diHiggsCandidate = 0;
   genWeights = 0;
   leadingPhotonID = 0;
   leadingPhotonISO = 0;
   subleadingPhotonID = 0;
   subleadingPhotonISO = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
//   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("genWeights", &genWeights, &b_genWeights);
   fChain->SetBranchAddress("genTotalWeight", &genTotalWeight, &b_genTotalWeight);
   fChain->SetBranchAddress("leadingPhoton", &leadingPhoton, &b_leadingPhoton);
   fChain->SetBranchAddress("leadingPhotonID", &leadingPhotonID, &b_leadingPhotonID);
   fChain->SetBranchAddress("leadingPhotonISO", &leadingPhotonISO, &b_leadingPhotonISO);
   fChain->SetBranchAddress("leadingPhotonEVeto", &leadingPhotonEVeto, &b_leadingPhotonEVeto);
   fChain->SetBranchAddress("subleadingPhoton", &subleadingPhoton, &b_subleadingPhoton);
   fChain->SetBranchAddress("subleadingPhotonID", &subleadingPhotonID, &b_subleadingPhotonID);
   fChain->SetBranchAddress("subleadingPhotonISO", &subleadingPhotonISO, &b_subleadingPhotonISO);
   fChain->SetBranchAddress("subleadingPhotonEVeto", &subleadingPhotonEVeto, &b_subleadingPhotonEVeto);
   fChain->SetBranchAddress("diphotonCandidate", &diphotonCandidate, &b_diphotonCandidate);
   fChain->SetBranchAddress("nPromptInDiPhoton", &nPromptInDiPhoton, &b_nPromptInDiPhoton);
   fChain->SetBranchAddress("leadingJet", &leadingJet, &b_leadingJet);
   fChain->SetBranchAddress("leadingJet_bDis", &leadingJet_bDis, &b_leadingJet_bDis);
   fChain->SetBranchAddress("subleadingJet", &subleadingJet, &b_subleadingJet);
   fChain->SetBranchAddress("subleadingJet_bDis", &subleadingJet_bDis, &b_subleadingJet_bDis);
   fChain->SetBranchAddress("dijetCandidate", &dijetCandidate, &b_dijetCandidate);
   fChain->SetBranchAddress("diHiggsCandidate", &diHiggsCandidate, &b_diHiggsCandidate);
   fChain->SetBranchAddress("isSignal", &isSignal, &b_isSignal);
   fChain->SetBranchAddress("isPhotonCR", &isPhotonCR, &b_isPhotonCR);
}

Bool_t bbggLTMaker::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef bbggLTMaker_cxx
