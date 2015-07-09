// -*- C++ -*-
//
// Package:    flashgg/bbggPlotter
// Class:      bbggPlotter
// 
/**\class bbggPlotter bbggPlotter.cc flashgg/bbggPlotter/plugins/bbggPlotter.cc

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

//Local
#include "flashgg/bbggTools/interface/bbggTools.h"

//
// class declaration
//

const int DEBUG = 0;

class bbggPlotter : public edm::EDAnalyzer {
   public:
      explicit bbggPlotter(const edm::ParameterSet&);
      ~bbggPlotter();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      typedef math::XYZTLorentzVector LorentzVector;


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
	  bbggTools tools_;
      //Parameter tokens
      edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diPhotonToken_;
      edm::EDGetTokenT<edm::View<flashgg::Jet> > thejetToken_;
      edm::InputTag rhoFixedGrid_;
      std::string bTagType;

      //Thresholds
      std::vector<double> ph_pt;
      std::vector<double> ph_eta;
      std::vector<double> ph_hoe;
      std::vector<double> ph_sieie;
      std::vector<double> ph_r9;
      std::vector<double> ph_chIso;
      std::vector<double> ph_nhIso;
      std::vector<double> ph_phIso;
      std::vector<int> ph_elVeto;
      std::vector<int> ph_doID;
      std::vector<int> ph_doISO;

      std::vector<double> diph_pt;
      std::vector<double> diph_eta;
      std::vector<double> diph_mass;
      unsigned int diph_onlyfirst;

      std::vector<double> jt_pt;
      std::vector<double> jt_eta;
      std::vector<double> jt_drPho;
      std::vector<double> jt_bDis;
      std::vector<int> jt_doPU;
      unsigned int n_bJets;

      std::vector<double> dijt_pt;
      std::vector<double> dijt_eta;
      std::vector<double> dijt_mass; 

      std::vector<double> cand_pt;
      std::vector<double> cand_eta;
      std::vector<double> cand_mass;
	  
	  std::vector<double> dr_cands;

      //OutFile & Hists
      TFile* outFile;
      std::map<std::string, TH1F> hists;
      std::map<std::string, TH2F> hists2D;
      std::string fileName;

      //Event counter for cout's
      long unsigned int EvtCount;
};

bbggPlotter::bbggPlotter(const edm::ParameterSet& iConfig) :
diPhotonToken_( consumes<edm::View<flashgg::DiPhotonCandidate> >( iConfig.getUntrackedParameter<edm::InputTag> ( "DiPhotonTag", edm::InputTag( "flashggDiPhotons" ) ) ) ),
thejetToken_( consumes<edm::View<flashgg::Jet> >( iConfig.getUntrackedParameter<edm::InputTag>( "JetTag", edm::InputTag( "flashggJets" ) ) ) )
{
   //now do what ever initialization is needed
	  tools_ = bbggTools();
   	  EvtCount = 0;
//Default values for thresholds
      std::vector<double> def_ph_pt;
      std::vector<double> def_ph_eta;
      std::vector<double> def_ph_hoe;
      std::vector<double> def_ph_sieie;
      std::vector<double> def_ph_r9;
      std::vector<double> def_ph_chIso;
      std::vector<double> def_ph_nhIso;
      std::vector<double> def_ph_phIso;
      std::vector<int> def_ph_elVeto;
      std::vector<int> def_ph_doID;
      std::vector<int> def_ph_doISO;

      std::vector<double> def_diph_pt;
      std::vector<double> def_diph_eta;
      std::vector<double> def_diph_mass;
      unsigned int def_diph_onlyfirst;

      std::vector<double> def_jt_pt;
      std::vector<double> def_jt_eta;
      std::vector<double> def_jt_drPho;
      std::vector<double> def_jt_bDis;
      std::vector<int> def_jt_doPU;
      unsigned int def_n_bJets;

      std::vector<double> def_dijt_pt;
      std::vector<double> def_dijt_eta;
      std::vector<double> def_dijt_mass;

      std::vector<double> def_cand_pt;
      std::vector<double> def_cand_eta;
      std::vector<double> def_cand_mass;

	  std::vector<double> def_dr_cands;

      std::string def_bTagType;

      std::string def_fileName;

      def_ph_pt.push_back(10.);         def_ph_pt.push_back(10.);
      def_ph_eta.push_back(20.);         def_ph_eta.push_back(20.);
      def_ph_hoe.push_back(-1.);        def_ph_hoe.push_back(-1.);
      def_ph_sieie.push_back(-1.);      def_ph_sieie.push_back(-1.);
      def_ph_r9.push_back(-1.);         def_ph_r9.push_back(-1.);
      def_ph_chIso.push_back(-1.);      def_ph_chIso.push_back(-1.);
      def_ph_nhIso.push_back(-1.);      def_ph_nhIso.push_back(-1.);
      def_ph_phIso.push_back(-1.);      def_ph_phIso.push_back(-1.);
      def_ph_elVeto.push_back(-1.);     def_ph_elVeto.push_back(-1.);
      def_ph_doID.push_back(0);		def_ph_doID.push_back(0);
      def_ph_doISO.push_back(0);	def_ph_doISO.push_back(0);

      def_diph_pt.push_back(10.);       def_diph_pt.push_back(10.);
      def_diph_eta.push_back(0.);       def_diph_eta.push_back(0.);
      def_diph_mass.push_back(0.);      def_diph_mass.push_back(1000.);
      def_diph_onlyfirst = 0;
			
      def_jt_pt.push_back(10.);         def_jt_pt.push_back(10.);
      def_jt_eta.push_back(20.);         def_jt_eta.push_back(20.);
      def_jt_bDis.push_back(0.);        def_jt_bDis.push_back(0.);
      def_jt_doPU.push_back(0);         def_jt_doPU.push_back(0);

      def_n_bJets = 0;

      def_dijt_pt.push_back(10.);       def_dijt_pt.push_back(10.);
      def_dijt_eta.push_back(20.);       def_dijt_eta.push_back(20.);
      def_dijt_mass.push_back(0.);      def_dijt_mass.push_back(1000.);
      def_jt_drPho.push_back(0.5);

      def_cand_pt.push_back(0.);
      def_cand_eta.push_back(20.);
      def_cand_mass.push_back(0.);		def_cand_mass.push_back(2000.);
	  
	  def_dr_cands.push_back(0.11);


      def_bTagType = "pfCombinedInclusiveSecondaryVertexV2BJetTags";
      def_fileName =  "out.root";

//Get thresholds from config file
      ph_pt     = iConfig.getUntrackedParameter<std::vector<double > >("PhotonPtOverDiPhotonMass", def_ph_pt);
      ph_eta    = iConfig.getUntrackedParameter<std::vector<double > >("PhotonEta", def_ph_eta);
      ph_hoe    = iConfig.getUntrackedParameter<std::vector<double > >("PhotonHoverE", def_ph_hoe);
      ph_sieie  = iConfig.getUntrackedParameter<std::vector<double > >("PhotonSieie", def_ph_sieie);
      ph_r9     = iConfig.getUntrackedParameter<std::vector<double > >("PhotonR9", def_ph_r9);
      ph_chIso  = iConfig.getUntrackedParameter<std::vector<double > >("PhotonChargedIso", def_ph_chIso);
      ph_nhIso  = iConfig.getUntrackedParameter<std::vector<double > >("PhotonNeutralIso", def_ph_nhIso);
      ph_phIso  = iConfig.getUntrackedParameter<std::vector<double > >("PhotonPhotonIso", def_ph_phIso);
      ph_elVeto = iConfig.getUntrackedParameter<std::vector<int > >("PhotonElectronVeto", def_ph_elVeto);
      ph_doID   = iConfig.getUntrackedParameter<std::vector<int > >("PhotonDoID", def_ph_doID);
      ph_doISO  = iConfig.getUntrackedParameter<std::vector<int > >("PhotonDoISO", def_ph_doID);

      diph_pt   = iConfig.getUntrackedParameter<std::vector<double > >("DiPhotonPt", def_diph_pt);
      diph_eta  = iConfig.getUntrackedParameter<std::vector<double > >("DiPhotonEta", def_diph_eta);
      diph_mass = iConfig.getUntrackedParameter<std::vector<double > >("DiPhotonMassWindow", def_diph_mass);
      diph_onlyfirst = iConfig.getUntrackedParameter<unsigned int>("DiPhotonOnlyFirst", def_diph_onlyfirst);

      jt_pt     = iConfig.getUntrackedParameter<std::vector<double > >("JetPtOverDiJetMass", def_jt_pt);
      jt_eta    = iConfig.getUntrackedParameter<std::vector<double > >("JetEta", def_jt_eta);
      jt_drPho  = iConfig.getUntrackedParameter<std::vector<double > >("JetDrPho", def_jt_drPho);
      jt_bDis   = iConfig.getUntrackedParameter<std::vector<double > >("JetBDiscriminant", def_jt_bDis);
      jt_doPU   = iConfig.getUntrackedParameter<std::vector<int > >("JetDoPUID", def_jt_doPU);

      n_bJets = iConfig.getUntrackedParameter<unsigned int>("n_bJets", def_n_bJets);

      dijt_pt   = iConfig.getUntrackedParameter<std::vector<double > >("DiJetPt", def_dijt_pt);
      dijt_eta  = iConfig.getUntrackedParameter<std::vector<double > >("DiJetEta", def_dijt_eta);
      dijt_mass = iConfig.getUntrackedParameter<std::vector<double > >("DiJetMassWindow", def_dijt_mass);

      cand_pt 	= iConfig.getUntrackedParameter<std::vector<double > >("CandidatePt", def_cand_mass);
      cand_eta 	= iConfig.getUntrackedParameter<std::vector<double > >("CandidateEta", def_cand_mass);
      cand_mass = iConfig.getUntrackedParameter<std::vector<double > >("CandidateMassWindow", def_cand_mass);
	  
	  dr_cands  = iConfig.getUntrackedParameter<std::vector<double > >("CandidatesDeltaR", def_dr_cands);

      rhoFixedGrid_  = iConfig.getUntrackedParameter<edm::InputTag>( "rhoFixedGridCollection", edm::InputTag( "fixedGridRhoAll" ) );

      bTagType = iConfig.getUntrackedParameter<std::string>( "bTagType", def_bTagType );

      fileName = iConfig.getUntrackedParameter<std::string>( "OutFileName", def_fileName );
	  
	  tools_.SetCut_PhotonPtOverDiPhotonMass( ph_pt );
	  tools_.SetCut_PhotonEta( ph_eta );
	  tools_.SetCut_PhotonDoID( ph_doID );
	  tools_.SetCut_PhotonDoISO( ph_doISO );
	  tools_.SetCut_PhotonHoverE( ph_hoe );
	  tools_.SetCut_PhotonSieie( ph_sieie );
	  tools_.SetCut_PhotonR9( ph_r9 );
	  tools_.SetCut_PhotonChargedIso( ph_chIso );
	  tools_.SetCut_PhotonNeutralIso( ph_nhIso);
	  tools_.SetCut_PhotonPhotonIso( ph_phIso);
	  tools_.SetCut_PhotonElectronVeto( ph_elVeto );
	  tools_.SetCut_DiPhotonPt( diph_pt );
	  tools_.SetCut_DiPhotonEta( diph_eta );
	  tools_.SetCut_DiPhotonMassWindow( diph_mass );
	  tools_.SetCut_DiPhotonOnlyFirst( diph_onlyfirst );
	  tools_.SetCut_JetPt( jt_pt );
	  tools_.SetCut_JetEta( jt_eta );
	  tools_.SetCut_JetBDiscriminant( jt_bDis );
	  tools_.SetCut_JetDrPho( jt_drPho  );
	  tools_.SetCut_JetDoPUID( jt_doPU );
	  tools_.SetCut_n_bJets( n_bJets );
	  tools_.SetCut_DiJetPt( dijt_pt );
	  tools_.SetCut_DiJetEta( dijt_eta );
	  tools_.SetCut_DiJetMassWindow( dijt_mass );
	  tools_.SetCut_CandidateMassWindow( cand_mass );
	  tools_.SetCut_CandidatePt( cand_pt );
	  tools_.SetCut_CandidateEta( cand_eta );
	  tools_.SetCut_bTagType( bTagType );
	  tools_.SetCut_CandidatesDeltaR( dr_cands );
	  
	  
      std::cout << "Parameters initialized... cand_mass[0]: " << cand_mass[0] << "\t" << "cand_mass[1] " << cand_mass[1] << std::endl;

}


bbggPlotter::~bbggPlotter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
bbggPlotter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   if( EvtCount%100 == 0 ) std::cout << "[bbggPlotter::analyze] Analyzing event number: " << EvtCount << std::endl;
   EvtCount++;
   using namespace edm;
	             
   Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
   iEvent.getByToken( diPhotonToken_, diPhotons );
   Handle<View<flashgg::Jet> > theJets;
   iEvent.getByToken( thejetToken_, theJets );
   Handle<double> rhoHandle;        // the old way for now...move to getbytoken?
   iEvent.getByLabel( rhoFixedGrid_, rhoHandle );
   const double rhoFixedGrd = *( rhoHandle.product() );
   tools_.setRho(rhoFixedGrd);
   
	bool passedSelection = tools_.AnalysisSelection(diPhotons, theJets);
	if(!passedSelection) return;

	edm::Ptr<flashgg::DiPhotonCandidate> diphoCand = tools_.GetSelected_diphoCandidate();
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
 	
 	double pho1_chiso = tools_.getCHisoToCutValue( diphoCand, 0);
 	double pho2_chiso = tools_.getCHisoToCutValue( diphoCand, 1);
 	double pho1_nhiso = tools_.getNHisoToCutValue( diphoCand->leadingPhoton() );
 	double pho2_nhiso = tools_.getNHisoToCutValue( diphoCand->subLeadingPhoton() );
 	double pho1_phiso = tools_.getPHisoToCutValue( diphoCand->leadingPhoton() );
 	double pho2_phiso = tools_.getPHisoToCutValue( diphoCand->subLeadingPhoton() ); 
 	//END photon Variables
	edm::Ptr<flashgg::Jet> LeadingJet = tools_.GetSelected_leadingJetCandidate();
	//Leading Jet Variables
	double jet1_pt = LeadingJet->pt();
	double jet1_eta = LeadingJet->eta();
	double jet1_phi = LeadingJet->phi();
	double jet1_bDis = LeadingJet->bDiscriminator(bTagType);
	double jet1_PUid = LeadingJet->passesPuJetId(diphoCand->vtx());
	double jet1_drPho1 = tools_.DeltaR(LeadingJet->p4(), diphoCand->leadingPhoton()->p4() );
	double jet1_drPho2 = tools_.DeltaR(LeadingJet->p4(), diphoCand->subLeadingPhoton()->p4() );
	//END Leading Jet Variables
	edm::Ptr<flashgg::Jet> SubLeadingJet = tools_.GetSelected_subleadingJetCandidate();
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

	/*
   bool isValidDiPhotonCandidate = false;

   edm::Ptr<reco::Vertex> CandVtx;
   edm::Ptr<flashgg::DiPhotonCandidate> diphoCand;

   //Begin DiPhoton Loop/Selection -----------------------------------------------------------
   for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ )
   {
	if(diph_onlyfirst && diphoIndex > 0 ) break;

	 edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( diphoIndex );
	 
	 dipho_pt = dipho->pt(); dipho_eta = dipho->eta(); dipho_phi = dipho->phi(); dipho_mass = dipho->mass();
		 
	 if(dipho_mass < diph_mass[0] || dipho_mass > diph_mass[1]) continue;
	 if(fabs(dipho_eta) > diph_eta[0] ) continue;
	 if(dipho_pt < diph_pt[0] ) continue;
		 
	 pho1_pt = dipho->leadingPhoton()->pt();			pho2_pt = dipho->subLeadingPhoton()->pt();
	 pho1_eta = dipho->leadingPhoton()->superCluster()->eta();	pho2_eta = dipho->subLeadingPhoton()->superCluster()->eta();
	 pho1_phi = dipho->leadingPhoton()->superCluster()->phi();	pho2_phi = dipho->subLeadingPhoton()->superCluster()->phi();
	 pho1_hoe = dipho->leadingPhoton()->hadronicOverEm(); 		pho2_hoe = dipho->subLeadingPhoton()->hadronicOverEm();
	 pho1_sieie = dipho->leadingPhoton()->full5x5_sigmaIetaIeta();	pho2_sieie = dipho->subLeadingPhoton()->full5x5_sigmaIetaIeta();
	 pho1_r9 = dipho->leadingPhoton()->r9();			pho2_r9 = dipho->subLeadingPhoton()->r9();
	 pho1_elveto = dipho->leadingPhoton()->passElectronVeto();	pho2_elveto = dipho->subLeadingPhoton()->passElectronVeto();

	 pho1_chiso = tools_.getCHisoToCutValue( dipho, 0);
	 pho2_chiso = tools_.getCHisoToCutValue( dipho, 1);
	 pho1_nhiso = tools_.getNHisoToCutValue( dipho->leadingPhoton() );
	 pho2_nhiso = tools_.getNHisoToCutValue( dipho->subLeadingPhoton() );
	 pho1_phiso = tools_.getPHisoToCutValue( dipho->leadingPhoton() );
	 pho2_phiso = tools_.getPHisoToCutValue( dipho->subLeadingPhoton() ); 
	 		 
	 if( pho1_pt < dipho_mass*ph_pt[0] ) continue;
	 if( fabs(pho1_eta) > ph_eta[1] ) continue;
	 if( pho2_pt < dipho->mass()*ph_pt[1] ) continue;
	 if( fabs(pho2_eta) > ph_eta[1] ) continue;
		 
	 bool pho1_id = true, pho2_id = true;
	 if( ph_doID[0] )
	 {
	   int pho1Index = 0;
	   if( fabs(pho1_eta) > ph_eta[0] ) pho1Index = 1;
			 
	   if( pho1_hoe > ph_hoe[pho1Index] ) 	pho1_id = false;
	   if( pho1_sieie > ph_sieie[pho1Index] ) pho1_id = false;
	   if( pho1_elveto != ph_elVeto[0] )    pho1_id = false;
	}
	if(ph_doISO[0])
	{
	   int pho1Index = 0;
           if( fabs(pho1_eta) > ph_eta[0] ) pho1Index = 1;

	   if( pho1_chiso > ph_chIso[pho1Index] ) pho1_id = false;
	   if( pho1_nhiso > ph_nhIso[pho1Index] ) pho1_id = false;
	   if( pho1_phiso > ph_phIso[pho1Index] ) pho1_id = false;
	}
	if( ph_doID[1] )
	{
	   int pho2Index = 2;
	   if( fabs(pho2_eta) > ph_eta[0] ) pho2Index = 3;

           if( pho2_hoe > ph_hoe[pho2Index] )   pho2_id = false;
           if( pho2_sieie > ph_sieie[pho2Index] ) pho2_id = false;
	   if( pho2_elveto != ph_elVeto[1] )    pho2_id = false;
	}
	if(ph_doISO[1])
	{
           int pho2Index = 2;
           if( fabs(pho2_eta) > ph_eta[0] ) pho2Index = 3;	 

	   if( pho2_chiso > ph_chIso[pho2Index] ) pho2_id = false;
	   if( pho2_nhiso > ph_nhIso[pho2Index] ) pho2_id = false;
	   if( pho2_phiso > ph_phIso[pho2Index] ) pho2_id = false;
	}

	if(pho1_id == true && pho2_id == true){
	  isValidDiPhotonCandidate = true;
  	  CandVtx = dipho->vtx();
	  diphoCand = dipho;
	  break;
	}

   }
   if( isValidDiPhotonCandidate == false ) return;
   if(DEBUG) std::cout << "Passed diphoton selection..." << std::endl;
   //End DiPhoton Loop/Selection -----------------------------------------------------------
   
   //Begin Jets Loop/Selection ------------------------------------------------------------
   std::vector<edm::Ptr<flashgg::Jet>> Jets;
   int nJet1 = 0, nJet2 = 0;
   int nbjets = 0;
   for( unsigned int jetIndex = 0; jetIndex < theJets->size(); jetIndex++ )
   {
   	edm::Ptr<flashgg::Jet> jet = theJets->ptrAt( jetIndex );
   	bool isJet1 = true, isJet2 = true;
	
   	if(jet->pt() < jt_pt[0]) isJet1 = false;
   	if(fabs(jet->eta()) > jt_eta[0] ) isJet1 = false;
   	if( jt_doPU[0] && jet->passesPuJetId(CandVtx) == 0 ) isJet1 = false;
	
   	if(jet->pt() < jt_pt[1]) isJet2 = false;
   	if(fabs(jet->eta()) > jt_eta[1] ) isJet2 = false;
   	if( jt_doPU[1] && jet->passesPuJetId(CandVtx) == 0 ) isJet2 = false;
	
   	if(isJet1) nJet1++;
   	if(isJet1 == false && isJet2) nJet2++;

	if(jet->bDiscriminator(bTagType) < jt_bDis[0]) continue;
	if( !isJet1 && !isJet2 ) continue;
	if( tools_.DeltaR(jet->p4(), diphoCand->leadingPhoton()->p4()) < jt_drPho[0] 
            || tools_.DeltaR(jet->p4(), diphoCand->subLeadingPhoton()->p4()) < jt_drPho[0] ) continue;
	
	Jets.push_back(jet);
	if( jet->bDiscriminator(bTagType) > jt_bDis[1] ) nbjets++;
   }

   if(Jets.size() < 2 ) return;

   edm::Ptr<flashgg::Jet> jet1, jet2;
   bbggPlotter::LorentzVector DiJet(0,0,0,0);
   double dijetPt_ref = 0;
   bool hasDiJet = false;

   for(unsigned int iJet = 0; iJet < Jets.size(); iJet++)
   {
	unsigned int isbjet = 0;
	if( Jets[iJet]->bDiscriminator(bTagType) > jt_bDis[1] ) isbjet = 1;
	for(unsigned int jJet = iJet+1; jJet < Jets.size(); jJet++)
	{
	  unsigned int isbjet2 = 0;
	  if( Jets[jJet]->bDiscriminator(bTagType) > jt_bDis[1] ) isbjet2 = 1;
	  
	  unsigned int totalbjet = isbjet + isbjet2;
	  if(n_bJets && totalbjet != n_bJets) continue;
//	  if(totalbjet != n_bJets) continue;

	  bbggPlotter::LorentzVector dijet = Jets[iJet]->p4() + Jets[jJet]->p4();
	  if(dijet.mass() < dijt_mass[0] || dijet.mass() > dijt_mass[1]) continue;
	  if(dijet.pt() > dijetPt_ref && dijet.pt() > dijt_pt[0] && fabs(dijet.Eta()) < dijt_eta[0] )
	  {
	      hasDiJet = true;
              dijetPt_ref = dijet.pt();
              DiJet = dijet;
              if( Jets[iJet]->pt() > Jets[jJet]->pt() ) {
                     jet1 = Jets.at(iJet);
                     jet2 = Jets.at(jJet);
              } else {
                     jet2 = Jets.at(iJet);
                     jet1 = Jets.at(jJet);
              } 
	  }
	}

   }


   if( hasDiJet == false ) return;
   
   dijet_pt = DiJet.pt();
   dijet_eta = DiJet.eta();
   dijet_phi = DiJet.phi();
   dijet_mass = DiJet.mass();
	
   jet1_pt = jet1->pt();
   jet1_eta = jet1->eta();
   jet1_phi = jet1->phi();
   jet1_bDis = jet1->bDiscriminator(bTagType);
   jet1_PUid = jet1->passesPuJetId(CandVtx);

	
   jet2_pt = jet2->pt();
   jet2_eta = jet2->eta();
   jet2_phi = jet2->phi();
   jet2_bDis = jet2->bDiscriminator(bTagType);
   jet2_PUid = jet2->passesPuJetId(CandVtx);

   jet1_drPho1 = tools_.DeltaR(jet1->p4(), diphoCand->leadingPhoton()->p4());
   jet1_drPho2 = tools_.DeltaR(jet1->p4(), diphoCand->subLeadingPhoton()->p4());
   jet2_drPho1 = tools_.DeltaR(jet2->p4(), diphoCand->leadingPhoton()->p4());
   jet2_drPho2 = tools_.DeltaR(jet2->p4(), diphoCand->subLeadingPhoton()->p4());

   //End Jets Loop/Selection -----------------------------------------------------------
   
   //Candidate assignment --------------------------------------------------------------
   double deltaR = tools_.DeltaR(DiJet, diphoCand->p4());
   double dijet_sumbdis = jet1_bDis + jet2_bDis;
   bbggPlotter::LorentzVector HHCandidate = DiJet + diphoCand->p4();
   cand4_pt = HHCandidate.pt();
   cand4_eta = HHCandidate.eta();
   cand4_phi = HHCandidate.phi();
   cand4_mass = HHCandidate.mass();
   if(DEBUG) std::cout << cand4_pt << "\t < \t" << cand_pt[0] << std::endl;
   if(cand4_pt < cand_pt[0] ) return;
   if(DEBUG) std::cout << fabs(cand4_eta) << "\t > \t" << cand_eta[0] << std::endl;
   if(fabs(cand4_eta) > cand_eta[0] ) return;
   if(DEBUG) std::cout << cand4_mass << "\t" << cand_mass[0] << "\t" << cand_mass[1] << std::endl;
   if(cand4_mass < cand_mass[0] || cand4_mass > cand_mass[1] ) return;
   if(DEBUG) std::cout << "Passed 4-candidate selection..." << std::endl;
   //END Candidate assignment ---------------------------------------------------------- 
*/
	
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
   hists["pho1_chiso"].Fill(pho1_chiso);
   hists["pho1_nhiso"].Fill(pho1_nhiso);
   hists["pho1_phiso"].Fill(pho1_phiso);
   hists["pho1_elveto"].Fill(pho1_elveto);
	
   hists["pho2_pt"].Fill(pho2_pt);
   hists["pho2_eta"].Fill(pho2_eta);
   hists["pho2_phi"].Fill(pho2_phi);
	
   hists["pho2_hoe"].Fill(pho2_hoe);
   hists["pho2_sieie"].Fill(pho2_sieie);
   hists["pho2_r9"].Fill(pho2_r9);
   hists["pho2_chiso"].Fill(pho2_chiso);
   hists["pho2_nhiso"].Fill(pho2_nhiso);
   hists["pho2_phiso"].Fill(pho2_phiso);
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
bbggPlotter::beginJob()
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
	hists2D["dijetmass_lJpt"] = TH2F("dijetmass_minDeltaR_lJP", "DiJet Mass vs LeadingJet Pt (GeV); M(jj); p_{T}(leading jet) (GeV)", 100, 50, 300, 100, 0, 300);
	hists2D["dijetmass_sJpt"] = TH2F("dijetmass_minDeltaR_sJP", "DiJet Mass vs subLeadingJet Pt (GeV); M(jj); p_{T}(leading jet) (GeV)", 100, 50, 300, 100, 0, 300);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
bbggPlotter::endJob() 
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
bbggPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(bbggPlotter);
