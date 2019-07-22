#include <iostream>
#include <cstdint>
#include <sstream>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>
#include <TLeafD.h>

using namespace std;

struct variables
{
    Float_t subleadingPhoton_pt;
    Float_t leadingPhoton_pt;
    Float_t subleadingPhoton_eta;
    Float_t leadingPhoton_eta;
    Float_t subleadingJet_pt;
    Float_t leadingJet_pt;
    Float_t subleadingJet_eta;
    Float_t leadingJet_eta;
    Float_t leadingJet_DeepCSV;
    Float_t PhoJetMinDr;
    Float_t subleadingJet_DeepCSV;
    Float_t dijetCandidatePtOverdiHiggsM;
    Float_t absCosTheta_bb;
    Float_t HHbbggMVA;
    Float_t sigmaMJets;
    Float_t PhoJetotherDr;
};

int main()
{
//     variables var;
    Float_t subleadingPhoton_pt, leadingPhoton_pt,subleadingPhoton_eta,leadingPhoton_eta,subleadingJet_pt, leadingJet_pt, subleadingJet_eta, leadingJet_eta, leadingJet_DeepCSV,  PhoJetMinDr,subleadingJet_DeepCSV, dijetCandidatePtOverdiHiggsM, absCosTheta_bb, HHbbggMVA, absCosThetaStar_CS, sigmaMJets, PhoJetotherDr;
    
 //   TFile *file1= new TFile("output_GluGluHToGG_M-125_13TeV_powheg_pythia8_2016.root", "READ");
 //   TFile *file2= new TFile("output_GluGluHToGG_M-125_13TeV_powheg_pythia8.root", "READ");
    TFile *file1= new TFile("output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph_2016.root", "READ");
    TFile *file2= new TFile("output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root", "READ");
  //  TFile *file1= new TFile("Bkg_merged.root", "READ");
  // TFile *file2= new TFile("Signal_merged.root", "READ");
    TFile *file_res= new TFile("Signal_merged.root", "RECREATE");
    if (file1->IsZombie() || file2->IsZombie()) cout<<"files didn't opened properly"<<endl;
   // TTree *tree1 = (TTree*) file1->Get("tagsDumper/trees/GluGluHToGG_M_125_13TeV_powheg_pythia8_13TeV_DoubleHTag_0");
  //  TTree *tree1 = (TTree*) file1->Get("tagsDumper/trees/GluGluToHHTo2B2G_node_SM_13TeV_madgraph_13TeV_DoubleHTag_0");
//    TTree *tree1 = (TTree*) file1->Get("bbggSelectionTree");
    TTree *tree1 = (TTree*) file1->Get("bbggSelectionTree");
    TTree *tree2 = (TTree*) file2->Get("bbggSelectionTree");
    if (tree1 == nullptr) cout<<"problems with tree1"<<endl;
    if( tree2 == nullptr) cout<<"problems with tree2"<<endl;
    TTree *tree_res = new TTree("bbggSelectionTree","");
    
    tree_res->Branch("subleadingPhoton_pt", &subleadingPhoton_pt, "subleadingPhoton_pt/F");
    tree_res->Branch("leadingPhoton_pt", &leadingPhoton_pt, "leadingPhoton_pt/F");
    tree_res->Branch("subleadingPhoton_eta", &subleadingPhoton_eta, "subleadingPhoton_eta/F");
    tree_res->Branch("leadingPhoton_eta", &leadingPhoton_eta, "leadingPhoton_eta/F");
    tree_res->Branch("subleadingJet_pt", &subleadingJet_pt, "subleadingJet_pt/F");
    tree_res->Branch("leadingJet_pt", &leadingJet_pt, "leadingJet_pt/F");
    tree_res->Branch("subleadingJet_eta", &subleadingJet_eta, "subleadingPhoton_eta/F");
    tree_res->Branch("leadingJet_eta", &leadingJet_eta, "leadingJet_eta/F");
    tree_res->Branch("leadingJet_DeepCSV", &leadingJet_DeepCSV, "leadingJet_DeepCSV/F");
    tree_res->Branch("PhoJetMinDr", &PhoJetMinDr, "PhoJetMinDr/F");
    tree_res->Branch("subleadingJet_DeepCSV", &subleadingJet_DeepCSV, "subleadingJet_DeepCSV/F");
    tree_res->Branch("dijetCandidatePtOverdiHiggsM", &dijetCandidatePtOverdiHiggsM, "dijetCandidatePtOverdiHiggsM/F");
    tree_res->Branch("absCosTheta_bb", &absCosTheta_bb, "absCosTheta_bb/F");
    tree_res->Branch("absCosThetaStar_CS", &absCosThetaStar_CS, "absCosThetaStar_CS/F");
    tree_res->Branch("HHbbggMVA", &HHbbggMVA, "HHbbggMVA/F");
    tree_res->Branch("sigmaMJets", &sigmaMJets, "sigmaMJets/F");
    tree_res->Branch("PhoJetotherDr", &PhoJetotherDr, "PhoJetotherDr/F");
    
    
    long int numpart1 = tree1->GetEntries();
    long int numpart2 = tree2->GetEntries();
    cout<<numpart1<<" "<<numpart2<<endl;
    cout<<"summ: "<<numpart1 + numpart2<<endl;
    for (int i=0; i<numpart1; i++)
    {
        tree1->GetEntry(i);
        subleadingPhoton_pt = ((TLeafD*)(tree1->GetLeaf("subleadingPhoton_pt")))->GetValue();
        leadingPhoton_pt = ((TLeafD*)(tree1->GetLeaf("leadingPhoton_pt")))->GetValue();
        subleadingPhoton_eta = ((TLeafD*)(tree1->GetLeaf("subleadingPhoton_eta")))->GetValue();
        leadingPhoton_eta = ((TLeafD*)(tree1->GetLeaf("leadingPhoton_eta")))->GetValue();
        subleadingJet_pt = ((TLeafD*)(tree1->GetLeaf("subleadingJet_pt")))->GetValue();
        leadingJet_pt = ((TLeafD*)(tree1->GetLeaf("leadingJet_pt")))->GetValue();
        subleadingJet_eta = ((TLeafD*)(tree1->GetLeaf("subleadingJet_eta")))->GetValue();
        leadingJet_eta = ((TLeafD*)(tree1->GetLeaf("leadingJet_eta")))->GetValue();
        leadingJet_DeepCSV =((TLeafD*)(tree1->GetLeaf("leadingJet_DeepCSV")))->GetValue();
        PhoJetMinDr=((TLeafD*)(tree1->GetLeaf("PhoJetMinDr")))->GetValue();
        subleadingJet_DeepCSV=((TLeafD*)(tree1->GetLeaf("subleadingJet_DeepCSV")))->GetValue();
        dijetCandidatePtOverdiHiggsM=((TLeafD*)(tree1->GetLeaf("dijetCandidatePtOverdiHiggsM")))->GetValue();
        absCosTheta_bb=((TLeafD*)(tree1->GetLeaf("absCosTheta_bb")))->GetValue(); 
        absCosThetaStar_CS=((TLeafD*)(tree1->GetLeaf("absCosThetaStar_CS")))->GetValue(); 
        HHbbggMVA=((TLeafD*)(tree1->GetLeaf("HHbbggMVA")))->GetValue(); 
        sigmaMJets=((TLeafD*)(tree1->GetLeaf("sigmaMJets")))->GetValue(); 
        PhoJetotherDr=((TLeafD*)(tree1->GetLeaf("PhoJetotherDr")))->GetValue(); 
        tree_res->Fill();
    }
    
    for (int i=0; i<numpart2; i++)
    {
        tree2->GetEntry(i);
        
//	cout<<((TLeafD*)(tree2->GetLeaf("absCosThetaStar_CS")))->GetValue()<<endl;        
        subleadingPhoton_pt = ((TLeafD*)(tree2->GetLeaf("subleadingPhoton_pt")))->GetValue();
        leadingPhoton_pt = ((TLeafD*)(tree2->GetLeaf("leadingPhoton_pt")))->GetValue();
        subleadingPhoton_eta = ((TLeafD*)(tree2->GetLeaf("subleadingPhoton_eta")))->GetValue();
        leadingPhoton_eta = ((TLeafD*)(tree2->GetLeaf("leadingPhoton_eta")))->GetValue();
        subleadingJet_pt = ((TLeafD*)(tree2->GetLeaf("subleadingJet_pt")))->GetValue();
        leadingJet_pt = ((TLeafD*)(tree2->GetLeaf("leadingJet_pt")))->GetValue();
        subleadingJet_eta = ((TLeafD*)(tree2->GetLeaf("subleadingJet_eta")))->GetValue();
        leadingJet_eta = ((TLeafD*)(tree2->GetLeaf("leadingJet_eta")))->GetValue();
        leadingJet_DeepCSV =((TLeafD*)(tree2->GetLeaf("leadingJet_DeepCSV")))->GetValue();
        PhoJetMinDr=((TLeafD*)(tree2->GetLeaf("PhoJetMinDr")))->GetValue();
        subleadingJet_DeepCSV=((TLeafD*)(tree2->GetLeaf("subleadingJet_DeepCSV")))->GetValue();
        dijetCandidatePtOverdiHiggsM=((TLeafD*)(tree2->GetLeaf("dijetCandidatePtOverdiHiggsM")))->GetValue();
        absCosTheta_bb=((TLeafD*)(tree2->GetLeaf("absCosTheta_bb")))->GetValue(); 
        absCosThetaStar_CS=((TLeafD*)(tree2->GetLeaf("absCosThetaStar_CS")))->GetValue(); 
        HHbbggMVA=((TLeafD*)(tree2->GetLeaf("HHbbggMVA")))->GetValue(); 
        sigmaMJets=((TLeafD*)(tree2->GetLeaf("sigmaMJets")))->GetValue(); 
        PhoJetotherDr=((TLeafD*)(tree2->GetLeaf("PhoJetotherDr")))->GetValue(); 
        tree_res->Fill();
    }
    long int numpart = tree_res->GetEntries();
    cout<<numpart<<" entries"<<endl;
    file_res->Write();
    
    TBranch* branch = tree_res->GetBranch("merged tree");
    if (branch) cout<<"success"<<endl;
    delete tree1;
    delete tree2;
    delete tree_res;
    
    delete file1;
    delete file2;
    delete file_res;
    
}
