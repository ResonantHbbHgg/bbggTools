// -*- C++ -*-
//
// Package:    flashgg/bbggMCMatch
// Class:      bbggMCMatch
// 
/**\class bbggMCMatch bbggMCMatch.cc flashgg/bbggMCMatch/plugins/bbggMCMatch.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rafael Teixeira De Lima
//         Created:  Tue, 16 Jun 2015 17:11:20 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <string>
#include <utility>
#include <cmath>

//root include files
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//FLASHgg files
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//Local
#include "flashgg/bbggTools/interface/bbggTools.h"
#include "flashgg/bbggTools/interface/bbggMC.h"
//
// class declaration
//

const int DEBUG = 0;

class bbggMCMatch : public edm::EDAnalyzer {
   public:
      explicit bbggMCMatch(const edm::ParameterSet&);
      ~bbggMCMatch();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      typedef math::XYZTLorentzVector LorentzVector;


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
	  bbggMC	mcTools_;
	  bbggTools tools_;
      //Parameter tokens
      edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diPhotonToken_;
      edm::EDGetTokenT<edm::View<flashgg::Jet> > thejetToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genToken_;
      edm::InputTag rhoFixedGrid_;
      std::string bTagType;

      //OutFile & Hists
      TFile* outFile;
      std::map<std::string, TH1F> hists;
      std::map<std::string, TH2F> hists2D;
      std::string fileName;

      //Event counter for cout's
      long unsigned int EvtCount;
};

bbggMCMatch::bbggMCMatch(const edm::ParameterSet& iConfig) :
diPhotonToken_( consumes<edm::View<flashgg::DiPhotonCandidate> >( iConfig.getUntrackedParameter<edm::InputTag> ( "DiPhotonTag", edm::InputTag( "flashggDiPhotons" ) ) ) ),
thejetToken_( consumes<edm::View<flashgg::Jet> >( iConfig.getUntrackedParameter<edm::InputTag>( "JetTag", edm::InputTag( "flashggJets" ) ) ) ),
genToken_( consumes<edm::View<reco::GenParticle> >( iConfig.getUntrackedParameter<edm::InputTag>( "GenTag", edm::InputTag( "prunedGenParticles" ) ) ) )
{
   //now do what ever initialization is needed
	  mcTools_ = bbggMC();
	  tools_ = bbggTools();
   	  EvtCount = 0;
      std::string def_fileName =  "out.root";
	  std::string def_bTagType = "pfCombinedInclusiveSecondaryVertexV2BJetTags";
      rhoFixedGrid_  = iConfig.getUntrackedParameter<edm::InputTag>( "rhoFixedGridCollection", edm::InputTag( "fixedGridRhoAll" ) );
      bTagType = iConfig.getUntrackedParameter<std::string>( "bTagType", def_bTagType );
      fileName = iConfig.getUntrackedParameter<std::string>( "OutFileName", def_fileName );	  
	  
      std::cout << "Parameters initialized... out file name: " << fileName << std::endl;

}


bbggMCMatch::~bbggMCMatch()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
bbggMCMatch::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   if( EvtCount%100 == 0 ) std::cout << "[bbggMCMatch::analyze] Analyzing event number: " << EvtCount << std::endl;
   EvtCount++;
   using namespace edm;
	             
   Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
   iEvent.getByToken( diPhotonToken_, diPhotons );
   Handle<View<flashgg::Jet> > theJets;
   iEvent.getByToken( thejetToken_, theJets );
   Handle<View<reco::GenParticle> > genParticles;
   iEvent.getByToken( genToken_, genParticles);
   Handle<double> rhoHandle;        // the old way for now...move to getbytoken?
   iEvent.getByLabel( rhoFixedGrid_, rhoHandle );
   const double rhoFixedGrd = *( rhoHandle.product() );
   tools_.setRho(rhoFixedGrd);
   
	bool passedSelection = mcTools_.MatchTruth(diPhotons, theJets, genParticles);
	if(!passedSelection) return;
		
	edm::Ptr<flashgg::DiPhotonCandidate> diphoCand = mcTools_.GetSelected_diphoCandidate();
	//Photon Variables
 	double dipho_pt = diphoCand->pt();
 	double dipho_eta = diphoCand->eta();
 	double dipho_phi = diphoCand->phi();
 	double dipho_mass = diphoCand->mass();
 	
 	double pho1_pt = diphoCand->leadingPhoton()->pt();		
 	double pho2_pt = diphoCand->subLeadingPhoton()->pt();
 	double pho1_eta = diphoCand->leadingPhoton()->superCluster()->eta();
 	double pho2_eta = diphoCand->subLeadingPhoton()->superCluster()->eta();
 	double pho1_phi = diphoCand->leadingPhoton()->superCluster()->phi();
 	double pho2_phi = diphoCand->subLeadingPhoton()->superCluster()->phi();
 	double pho1_hoe = diphoCand->leadingPhoton()->hadronicOverEm(); 	
 	double pho2_hoe = diphoCand->subLeadingPhoton()->hadronicOverEm();
 	double pho1_sieie = diphoCand->leadingPhoton()->full5x5_sigmaIetaIeta();
 	double pho2_sieie = diphoCand->subLeadingPhoton()->full5x5_sigmaIetaIeta();
 	double pho1_r9 = diphoCand->leadingPhoton()->r9();		
 	double pho2_r9 = diphoCand->subLeadingPhoton()->r9();
 	double pho1_elveto = diphoCand->leadingPhoton()->passElectronVeto();
 	double pho2_elveto = diphoCand->subLeadingPhoton()->passElectronVeto();
 	
     // double pho1_chiso = tools_.getCHisoToCutValue( diphoCand, 0);
     // double pho2_chiso = tools_.getCHisoToCutValue( diphoCand, 1);
     // double pho1_nhiso = tools_.getNHisoToCutValue( diphoCand->leadingPhoton() );
     // double pho2_nhiso = tools_.getNHisoToCutValue( diphoCand->subLeadingPhoton() );
     // double pho1_phiso = tools_.getPHisoToCutValue( diphoCand->leadingPhoton() );
     // double pho2_phiso = tools_.getPHisoToCutValue( diphoCand->subLeadingPhoton() );
 	//END photon Variables
	edm::Ptr<flashgg::Jet> LeadingJet = mcTools_.GetSelected_leadingJetCandidate();
	//Leading Jet Variables
	double jet1_pt = LeadingJet->pt();
	double jet1_eta = LeadingJet->eta();
	double jet1_phi = LeadingJet->phi();
	double jet1_bDis = LeadingJet->bDiscriminator(bTagType);
	double jet1_PUid = LeadingJet->passesPuJetId(diphoCand->vtx());
	double jet1_drPho1 = tools_.DeltaR(LeadingJet->p4(), diphoCand->leadingPhoton()->p4() );
	double jet1_drPho2 = tools_.DeltaR(LeadingJet->p4(), diphoCand->subLeadingPhoton()->p4() );
	//END Leading Jet Variables
	edm::Ptr<flashgg::Jet> SubLeadingJet = mcTools_.GetSelected_subleadingJetCandidate();
	//SubLeading Jet Variables
	double jet2_pt = SubLeadingJet->pt();
	double jet2_eta = SubLeadingJet->eta();
	double jet2_phi = SubLeadingJet->phi();
	double jet2_bDis = SubLeadingJet->bDiscriminator(bTagType);
	double jet2_PUid = SubLeadingJet->passesPuJetId(diphoCand->vtx());
	double jet2_drPho1 = tools_.DeltaR(SubLeadingJet->p4(), diphoCand->leadingPhoton()->p4() );
	double jet2_drPho2 = tools_.DeltaR(SubLeadingJet->p4(), diphoCand->subLeadingPhoton()->p4() );
	//END SubLeading Jet Variables
	//DiJet Variables
	LorentzVector DiJet = LeadingJet->p4() + SubLeadingJet->p4();
    double dijet_pt = DiJet.pt();
 	double dijet_eta = DiJet.eta();
 	double dijet_phi = DiJet.phi();
 	double dijet_mass = DiJet.mass();
    double deltaR = tools_.DeltaR(DiJet, diphoCand->p4());
    double dijet_sumbdis = jet1_bDis + jet2_bDis;
	//END DiJet Variables
	//HH Candidate Variables
	LorentzVector HHCandidate = DiJet + diphoCand->p4();
    double cand4_pt = HHCandidate.pt();
 	double cand4_eta = HHCandidate.eta();
 	double cand4_phi = HHCandidate.phi();
 	double cand4_mass = HHCandidate.mass();
	//END HH Candidate Variables
	//Other variables
	double candmass_minDeltaR_lJP = (jet1_drPho1 < jet1_drPho2) ? jet1_drPho1 : jet1_drPho2;
	double candmass_minDeltaR_sJP = (jet2_drPho1 < jet2_drPho2) ? jet2_drPho1 : jet2_drPho2;

	if(candmass_minDeltaR_lJP < 0 || candmass_minDeltaR_sJP < 0 ) std::cout << "NEGATIVE!!!! " << candmass_minDeltaR_lJP << "\t" << candmass_minDeltaR_sJP << std::endl;

	int nbjets = -1;
	
   if(DEBUG) std::cout << "GOT TO THE END!!" << std::endl;

   hists2D["jet1bdis_jet2bdis"].Fill(jet1_bDis, jet2_bDis);

   hists["n_bjets"].Fill(nbjets);

   hists["dr_dijetdipho"].Fill( deltaR );

   hists2D["candmass_deltar"].Fill( cand4_mass, deltaR );
   hists2D["candmass_dijetmass"].Fill(cand4_mass, dijet_mass);
   hists2D["candmass_diphomass"].Fill(cand4_mass, dipho_mass);
   hists2D["candmass_minDeltaR_lJP"].Fill(cand4_mass, candmass_minDeltaR_lJP);
   hists2D["candmass_minDeltaR_sJP"].Fill(cand4_mass, candmass_minDeltaR_sJP);
	   
   hists["dipho_pt"].Fill(dipho_pt);
   hists["dipho_eta"].Fill(dipho_eta);
   hists["dipho_phi"].Fill(dipho_phi);
   hists["dipho_mass"].Fill(dipho_mass);

   hists2D["diphomass_lpho"].Fill(dipho_mass, pho1_pt);
	
   hists["dijet_pt"].Fill(dijet_pt);
   hists["dijet_eta"].Fill(dijet_eta);
   hists["dijet_phi"].Fill(dijet_phi);
   hists["dijet_mass"].Fill(dijet_mass);
   hists["dijet_sumbdis"].Fill(dijet_sumbdis);
   
   hists2D["dijetmass_deltar"].Fill(dijet_mass, deltaR);
   hists2D["dijetmass_minDeltaR_lJP"].Fill(dijet_mass, candmass_minDeltaR_lJP);
   hists2D["dijetmass_minDeltaR_sJP"].Fill(dijet_mass, candmass_minDeltaR_lJP);
   hists2D["dijetmass_lJpt"].Fill(dijet_mass, jet1_pt);
   hists2D["dijetmass_sJpt"].Fill(dijet_mass, jet2_pt);

   hists2D["dijetmass_sumbdis"].Fill(dijet_mass, dijet_sumbdis);
   hists2D["dijetmass_ratio_l"].Fill(dijet_mass, jet1_pt/dijet_mass);
   hists2D["dijetmass_ratio_s"].Fill(dijet_mass, jet2_pt/dijet_mass);
	
   hists["cand4_pt"].Fill(cand4_pt);
   hists["cand4_eta"].Fill(cand4_eta);
   hists["cand4_phi"].Fill(cand4_phi);
   hists["cand4_mass"].Fill(cand4_mass);
	
   hists["pho1_pt"].Fill(pho1_pt);
   hists["pho1_eta"].Fill(pho1_eta);
   hists["pho1_phi"].Fill(pho1_phi);
	
   hists["pho1_hoe"].Fill(pho1_hoe);
   hists["pho1_sieie"].Fill(pho1_sieie);
   hists["pho1_r9"].Fill(pho1_r9);
   // hists["pho1_chiso"].Fill(pho1_chiso);
   // hists["pho1_nhiso"].Fill(pho1_nhiso);
   // hists["pho1_phiso"].Fill(pho1_phiso);
   hists["pho1_elveto"].Fill(pho1_elveto);
	
   hists["pho2_pt"].Fill(pho2_pt);
   hists["pho2_eta"].Fill(pho2_eta);
   hists["pho2_phi"].Fill(pho2_phi);
	
   hists["pho2_hoe"].Fill(pho2_hoe);
   hists["pho2_sieie"].Fill(pho2_sieie);
   hists["pho2_r9"].Fill(pho2_r9);
   // hists["pho2_chiso"].Fill(pho2_chiso);
   // hists["pho2_nhiso"].Fill(pho2_nhiso);
   // hists["pho2_phiso"].Fill(pho2_phiso);
   hists["pho2_elveto"].Fill(pho2_elveto);
	
   hists["jet1_pt"].Fill(jet1_pt);
   hists["jet1_eta"].Fill(jet1_eta);
   hists["jet1_phi"].Fill(jet1_phi);
   hists["jet1_drPho1"].Fill(jet1_drPho1);
   hists["jet1_drPho2"].Fill(jet1_drPho2);
   hists["jet1_bDis"].Fill(jet1_bDis);
   hists["jet1_PUid"].Fill(jet1_PUid);
	
   hists["jet2_pt"].Fill(jet2_pt);
   hists["jet2_eta"].Fill(jet2_eta);
   hists["jet2_phi"].Fill(jet2_phi);
   hists["jet2_drPho1"].Fill(jet2_drPho1);
   hists["jet2_drPho2"].Fill(jet2_drPho2);
   hists["jet2_bDis"].Fill(jet2_bDis);
   hists["jet2_PUid"].Fill(jet2_PUid);


   if(DEBUG) std::cout << "Histograms filled!" << std::endl;
}

// ------------ method called once each job just before starting event loop  ------------
void 
bbggMCMatch::beginJob()
{
	outFile = new TFile(fileName.c_str(), "RECREATE");

	hists["dipho_pt"] 	= TH1F("dipho_pt", "DiPhoton p_{T}; p_{T}(#gamma#gamma) (GeV); Events", 100, 0, 400);
	hists["dipho_eta"] 	= TH1F("dipho_eta", "DiPhoton #eta; #eta(#gamma#gamma); Events", 100, -5, 5);
	hists["dipho_phi"]	= TH1F("dipho_phi", "DiPhoton #phi; #phi(#gamma#gamma); Events", 100, -3.5, 3.5);
	hists["dipho_mass"] 	= TH1F("dipho_mass", "DiPhoton Mass; M(#gamma#gamma); Events", 100, 100, 150);

	hists["dijet_pt"] 	= TH1F("dijet_pt", "DiJet p_{T}; p_{T}(jj) (GeV); Events", 100, 0, 500);
	hists["dijet_eta"] 	= TH1F("dijet_eta", "DiJet #eta; #eta(jj); Events", 100, -5, 5);
	hists["dijet_phi"]	= TH1F("dijet_phi", "DiJet #phi; #phi(jj); Events", 100, -3.5, 3.5);
	hists["dijet_mass"] 	= TH1F("dijet_mass", "DiJet Mass; M(jj); Events", 100, 30, 300);
	hists["dijet_sumbdis"]  = TH1F("dijet_sumbdis", "Sum of b-Discriminant of DiJet; Sum b-Dis; Events", 100, -0.1, 2.1); 

	hists["dr_dijetdipho"]	= TH1F("dr_dijetdipho", "DeltaR between DiJet and DiPhoton; #DeltaR(#gamma#gamma,jj); Events", 100, -1, 10);

	hists["cand4_pt"] 	= TH1F("cand_pt", "DiHiggs Candidate (jj#gamma#gamma) p_{T}; p_{T}(jj#gamma#gamma) (GeV); Events", 100, 0, 700);
	hists["cand4_eta"] 	= TH1F("cand_eta", "DiHiggs Candidate (jj#gamma#gamma) #eta; #eta(jj#gamma#gamma); Events", 100, -5, 5);
	hists["cand4_phi"]	= TH1F("cand_phi", "DiHiggs Candidate (jj#gamma#gamma) #phi; #phi(jj#gamma#gamma); Events", 100, -3.5, 3.5); 
	hists["cand4_mass"] 	= TH1F("cand_mass", "DiHiggs Candidate (jj#gamma#gamma) Mass; M(jj#gamma#gamma) (GeV); Events", 100, 100, 1000);

	hists["pho1_pt"] 	= TH1F("pho1_pt", "Leading Photon p_{T}; p_{T}(leading #gamma) (GeV); Events", 100, 10, 400);
	hists["pho1_eta"] 	= TH1F("pho1_eta", "Leading Photon #eta; #eta(leading #gamma); Events", 100, -5., 5.);
	hists["pho1_phi"]	= TH1F("pho1_phi", "Leading Photon #phi; #phi(leading #gamma); Events", 100, -3.5, 3.5);

	hists["pho1_hoe"] 	= TH1F("pho1_hoe", "Leading Photon H/E; H/E(leading #gamma); Events", 100, -0.01, 0.1);
	hists["pho1_sieie"] 	= TH1F("pho1_sieie", "Leading Photon #sigma_{i#etai#eta}; #sigma_{i#etai#eta}(leading #gamma); Events", 100, 0.0, 0.04);
	hists["pho1_r9"] 	= TH1F("pho1_r9", "Leading Photon R9; R9(leading #gamma); Events", 100, 0, 1.1);
	hists["pho1_chiso"] 	= TH1F("pho1_chiso", "Leading Photon Corrected Charged Isolation; Corrected Charged Isolation(leading #gamma); Events", 100, -1., 5);
	hists["pho1_nhiso"] 	= TH1F("pho1_nhiso", "Leading Photon Corrected Neutral Isolation; Corrected Neutral Isolation(leading #gamma); Events", 100, -1., 5);
	hists["pho1_phiso"] 	= TH1F("pho1_phiso", "Leading Photon Corrected Photon Isolation; Corrected Photon Isolation(leading #gamma); Events", 100, -1., 5);
	hists["pho1_elveto"] 	= TH1F("pho1_elveto", "Leading Photon Electron Veto; Electron Veto(leading #gamma); events", 8, -1, 3);

	hists["pho2_pt"] 	= TH1F("pho2_pt", "SubLeading Photon p_{T}; p_{T}(subLeading #gamma) (GeV); Events", 100, 10, 200);
	hists["pho2_eta"] 	= TH1F("pho2_eta", "SubLeading Photon #eta; #eta(subLeading #gamma); Events", 100, -5., 5.);
	hists["pho2_phi"]	= TH1F("pho2_phi", "SubLeading Photon #phi; #phi(subleading #gamma); Events", 100, -3.5, 3.5);

	hists["pho2_hoe"] 	= TH1F("pho2_hoe", "SubLeading Photon H/E; H/E(subLeading #gamma); Events", 100, -0.01, 0.1);
	hists["pho2_sieie"] 	= TH1F("pho2_sieie", "SubLeading Photon #sigma_{i#etai#eta}; #sigma_{i#etai#eta}(subLeading #gamma); Events", 100, 0.0, 0.04);
	hists["pho2_r9"] 	= TH1F("pho2_r9", "SubLeading Photon R9; R9(subLeading #gamma); Events", 100, 0, 1.1);
	hists["pho2_chiso"] 	= TH1F("pho2_chiso", "SubLeading Photon Corrected Charged Isolation; Corrected Charged Isolation(subLeading #gamma); Events", 100, -1., 5);
	hists["pho2_nhiso"] 	= TH1F("pho2_nhiso", "SubLeading Photon Corrected Neutral Isolation; Corrected Neutral Isolation(subLeading #gamma); Events", 100, -1., 5);
	hists["pho2_phiso"] 	= TH1F("pho2_phiso", "SubLeading Photon Corrected Photon Isolation; Corrected Photon Isolation(subLeading #gamma); Events", 100, -1., 5);
	hists["pho2_elveto"] 	= TH1F("pho2_elveto", "SubLeading Photon Electron Veto; Electron Veto(subLeading #gamma); events", 8, -1, 3);


	hists["jet1_pt"] 	= TH1F("jet1_pt", "Leading Jet p_{T}; p_{T}(leading jet) (GeV); Events", 100, 0, 300);
	hists["jet1_eta"] 	= TH1F("jet1_eta", "Leading Jet #eta; #eta(leading jet); Events", 100, -5., 5.);
	hists["jet1_phi"]	= TH1F("jet1_phi", "Leading Jet #phi; #phi(leading jet); Events", 100, -3.5, 3.5);
	hists["jet1_drPho1"]	= TH1F("jet1_drPho1", "Leading Jet DeltaR with Leading Photon; DeltaR(leading jet, leading #gamma); Events", 100, -1, 10);
	hists["jet1_drPho2"]	= TH1F("jet1_drPho2", "Leading Jet DeltaR with SubLeading Photon; DeltaR(leading jet, subleading #gamma); Events", 100, -1, 10);
	hists["jet1_bDis"] 	= TH1F("jet1_bDis", "Leading Jet b-Discriminant; b-Discriminant(leading jet); Events", 100, -0.01, 1.01);
	hists["jet1_PUid"] 	= TH1F("jet1_PUid", "Leading Jet PU ID; PU ID(leading jet); Events", 8, -1, 3);

	hists["jet2_pt"] 	= TH1F("jet2_pt", "SubLeading Jet p_{T}; p_{T}(subLeading jet) (GeV); Events", 100, 0, 300);
	hists["jet2_eta"] 	= TH1F("jet2_eta", "SubLeading Jet #eta; #eta(subLeading jet); Events", 100, -5., 5.);
	hists["jet2_phi"]	= TH1F("jet2_phi", "SubLeading Jet #phi; #phi(subleading jet); Events", 100, -3.5, 3.5);
	hists["jet2_drPho1"]	= TH1F("jet2_drPho1", "SubLeading Jet DeltaR with Leading Photon; DeltaR(subleading jet, leading #gamma); Events", 100, -1, 10);
	hists["jet2_drPho2"]	= TH1F("jet2_drPho2", "SubLeading Jet DeltaR with SubLeading Photon; DeltaR(subleadign jet, subleading #gamma); Events", 100, -1, 10);
	hists["jet2_bDis"] 	= TH1F("jet2_bDis", "SubLeading Jet b-Discriminant; b-Discriminant(subLeading jet); Events", 100, -0.01, 1.01);
	hists["jet2_PUid"] 	= TH1F("jet2_PUid", "SubLeading Jet PU ID; PU ID(subleading jet); Events", 8, -1, 3);

	hists["n_bjets"]	= TH1F("n_bjets", "Number of jets passing b-Dis requirement; #bJets; Events", 22, -1, 11);

	hists2D["candmass_dijetmass"] = TH2F("candmass_dijetmass", "DiHiggs Candidate Mass vs DiJet Mass; M(jj#gamma#gamma) (GeV); M(jj)", 100, 100, 800, 100, 50, 300);
	hists2D["candmass_diphomass"] = TH2F("candmass_diphomass", "DiHiggs Candidate Mass vs DiPhoton Mass; M(jj#gamma#gamma) (GeV); M(#gamma#gamma)", 100, 100, 800, 100, 110, 150);
	hists2D["candmass_deltar"]    = TH2F("candmass_deltar", "DiHiggs Candidate Mass vs DeltaR between DiCandidates; M(jj#gamma#gamma) (GeV); #DeltaR(#gamma#gamma,jj)", 100, 100, 800, 100, -1, 10);
	hists2D["candmass_minDeltaR_lJP"]= TH2F("candmass_minDeltaR_lJP", "DiHiggs Candidate Mass vs Minimum DeltaR between LeadingJet and Photons; M(jj#gamma#gamma) (GeV); Min(#DeltaR(Leading Jet,Leading Photon), #DeltaR(Leading Jet,subLeading Photon))", 100, 100, 800, 100, 0, 5 );
	hists2D["candmass_minDeltaR_sJP"]= TH2F("candmass_minDeltaR_sJP", "DiHiggs Candidate Mass vs Minimum DeltaR between subLeadingJet and Photons; M(jj#gamma#gamma) (GeV); Min(#DeltaR(subLeadingJet Jet,Leading Photon), #DeltaR(subLeadingJet Jet,subLeading Photon))", 100, 100, 800, 100, 0, 5 );
	hists2D["jet1bdis_jet2bdis"]    = TH2F("jet1bdis_jet2bdis", "Jet1 bDis vs Jet2 bDis; Jet1 bDis; Jet2 bDis", 100, 0., 1., 100, 0., 1.);

	hists2D["dijetmass_deltar"] 	= TH2F("dijetmass_deltar", "DiJet Mass vs DeltaR between DiCandidates; M(jj); #DeltaR(#gamma#gamma,jj)", 100, 50, 300, 100, -0.1, 10);
	hists2D["dijetmass_minDeltaR_lJP"] = TH2F("dijetmass_minDeltaR_lJP", "DiJet Mass vs Minimum DeltaR between LeadingJet and Photons; M(jj); Min(#DeltaR())", 100, 50, 300, 100, 0, 5 );
	hists2D["dijetmass_minDeltaR_sJP"] = TH2F("dijetmass_minDeltaR_sJP", "DiJet Mass vs Minimum DeltaR between subLeadingJet and Photons; M(jj); Min(#DeltaR())", 100, 50, 300, 100, 0, 5 );
	hists2D["dijetmass_lJpt"] = TH2F("dijetmass_lJpt", "DiJet Mass vs LeadingJet Pt (GeV); M(jj); p_{T}(leading jet) (GeV)", 100, 50, 300, 100, 0, 300);
	hists2D["dijetmass_sJpt"] = TH2F("dijetmass_sJpt", "DiJet Mass vs subLeadingJet Pt (GeV); M(jj); p_{T}(leading jet) (GeV)", 100, 50, 300, 100, 0, 300);

	hists2D["dijetmass_sumbdis"] = TH2F("dijetmass_sumbdis", "DiJet Mass vs Sum of b discriminant; M(jj); Sum of b discriminant", 100, 50, 300, 100, 0, 2.);
	hists2D["dijetmass_ratio_l"] = TH2F("dijetmass_ratio_l", "DiJet Mass vs Ratio of Leading Jet Pt and DiJet Mass; M(jj); Leading Jet Pt/DiJet Mass", 100, 50, 300, 100, 0, 3.);
	hists2D["dijetmass_ratio_s"] = TH2F("dijetmass_ratio_s", "DiJet Mass vs Ratio of subLeading Jet Pt and DiJet Mass; M(jj); subLeading Jet Pt/DiJet Mass", 100, 50, 300, 100, 0, 3.);



	hists2D["diphomass_lpho"] = TH2F("diphomass_lpho", "DiPhoton Mass vs Leading Photon Pt (GeV); M(#gamma#gamma); p_{T}(leading #gamma)", 100, 100, 150, 100, 10, 400);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
bbggMCMatch::endJob() 
{
	outFile->cd();
	for(std::map<std::string, TH1F>::iterator it = hists.begin(); it != hists.end(); ++it)
	{
		std::cout << "Saving histogram... " << it->first << std::endl;
		it->second.Write();
	}
	for(std::map<std::string, TH2F>::iterator it = hists2D.begin(); it != hists2D.end(); ++it)
	{
		std::cout << "Saving histogram... " << it->first << std::endl;
		it->second.Write();
	}
	outFile->Close();
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
bbggMCMatch::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(bbggMCMatch);
