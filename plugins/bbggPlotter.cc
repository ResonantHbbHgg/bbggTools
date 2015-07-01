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
#include "flashgg/bbggPlotter/interface/bbggTools.h"

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

      // ---------- extra methods... ---------------------
      double getCHisoToCutValue(edm::Ptr<flashgg::DiPhotonCandidate> dipho, int whichPho, float rho);
      double getNHisoToCutValue(const flashgg::Photon* pho, float rho);
      double getPHisoToCutValue(const flashgg::Photon* pho, float rho);
      double getEA( float eta, int whichEA);
      double DeltaR( bbggPlotter::LorentzVector vec1, bbggPlotter::LorentzVector vec2);


      // ----------member data ---------------------------
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

      rhoFixedGrid_  = iConfig.getUntrackedParameter<edm::InputTag>( "rhoFixedGridCollection", edm::InputTag( "fixedGridRhoAll" ) );

      bTagType = iConfig.getUntrackedParameter<std::string>( "bTagType", def_bTagType );

      fileName = iConfig.getUntrackedParameter<std::string>( "OutFileName", def_fileName );

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
	 
   double dipho_pt = -1, dipho_eta = -1, dipho_phi = -1, dipho_mass = -1;
   
   double dijet_pt = -1, dijet_eta = -1, dijet_phi = -1, dijet_mass = -1;
   
   double cand4_pt = -1, cand4_eta = -1, cand4_phi = -1, cand4_mass = -1;
   
   double pho1_pt = -1, pho1_eta = -1, pho1_phi = -1;
   double pho1_hoe = -1, pho1_sieie = -1, pho1_r9 = -1, pho1_chiso = -1, pho1_nhiso = -1, pho1_phiso = -1;
   int pho1_elveto = -1;
   
   double pho2_pt = -1, pho2_eta = -1, pho2_phi = -1;
   double pho2_hoe = -1, pho2_sieie = -1, pho2_r9 = -1, pho2_chiso = -1, pho2_nhiso = -1, pho2_phiso = -1;
   int pho2_elveto = -1;
   
   double jet1_pt = -1, jet1_eta = -1, jet1_phi = -1, jet1_bDis = -1, jet1_PUid = -1, jet1_drPho1 = -1, jet1_drPho2 = -1;
   double jet2_pt = -1, jet2_eta = -1, jet2_phi = -1, jet2_bDis = -1, jet2_PUid = -1, jet2_drPho1 = -1, jet2_drPho2 = -1;

   Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
   iEvent.getByToken( diPhotonToken_, diPhotons );
   Handle<View<flashgg::Jet> > theJets;
   iEvent.getByToken( thejetToken_, theJets );
   Handle<double> rhoHandle;        // the old way for now...move to getbytoken?
   iEvent.getByLabel( rhoFixedGrid_, rhoHandle );
   const double rhoFixedGrd = *( rhoHandle.product() );

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

	 pho1_chiso = bbggPlotter::getCHisoToCutValue( dipho, 0, rhoFixedGrd);
	 pho2_chiso = bbggPlotter::getCHisoToCutValue( dipho, 1, rhoFixedGrd);
	 pho1_nhiso = bbggPlotter::getNHisoToCutValue( dipho->leadingPhoton(), rhoFixedGrd );
	 pho2_nhiso = bbggPlotter::getNHisoToCutValue( dipho->subLeadingPhoton(), rhoFixedGrd );
	 pho1_phiso = bbggPlotter::getPHisoToCutValue( dipho->leadingPhoton(), rhoFixedGrd );
	 pho2_phiso = bbggPlotter::getPHisoToCutValue( dipho->subLeadingPhoton(), rhoFixedGrd );

		 
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
	if( bbggPlotter::DeltaR(jet->p4(), diphoCand->leadingPhoton()->p4()) < jt_drPho[0] 
            || bbggPlotter::DeltaR(jet->p4(), diphoCand->subLeadingPhoton()->p4()) < jt_drPho[0] ) continue;
	
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

   jet1_drPho1 = bbggPlotter::DeltaR(jet1->p4(), diphoCand->leadingPhoton()->p4());
   jet1_drPho2 = bbggPlotter::DeltaR(jet1->p4(), diphoCand->subLeadingPhoton()->p4());
   jet2_drPho1 = bbggPlotter::DeltaR(jet2->p4(), diphoCand->leadingPhoton()->p4());
   jet2_drPho2 = bbggPlotter::DeltaR(jet2->p4(), diphoCand->subLeadingPhoton()->p4());

   //End Jets Loop/Selection -----------------------------------------------------------
   
   //Candidate assignment --------------------------------------------------------------
   double deltaR = bbggPlotter::DeltaR(DiJet, diphoCand->p4());
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

   if(DEBUG) std::cout << "GOT TO THE END!!" << std::endl;

   hists2D["jet1bdis_jet2bdis"].Fill(jet1_bDis, jet2_bDis);

   hists["n_bjets"].Fill(nbjets);

   hists["dr_dijetdipho"].Fill( deltaR );

   hists2D["candmass_deltar"].Fill( cand4_mass, deltaR );
   hists2D["candmass_dijetmass"].Fill(cand4_mass, dijet_mass);
   hists2D["candmass_diphomass"].Fill(cand4_mass, dipho_mass);

   hists["dipho_pt"].Fill(dipho_pt);
   hists["dipho_eta"].Fill(dipho_eta);
   hists["dipho_phi"].Fill(dipho_phi);
   hists["dipho_mass"].Fill(dipho_mass);
	
   hists["dijet_pt"].Fill(dijet_pt);
   hists["dijet_eta"].Fill(dijet_eta);
   hists["dijet_phi"].Fill(dijet_phi);
   hists["dijet_mass"].Fill(dijet_mass);
   hists["dijet_sumbdis"].Fill(dijet_sumbdis);
	
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

// ------------ Extra methods... -------------
//

double bbggPlotter::DeltaR(bbggPlotter::LorentzVector vec1, bbggPlotter::LorentzVector vec2)
{
	double R2 = (vec1.phi() - vec2.phi())*(vec1.phi() - vec2.phi()) + (vec1.eta() - vec2.eta())*(vec1.eta() - vec2.eta());
	return sqrt(R2);
}

double bbggPlotter::getCHisoToCutValue(edm::Ptr<flashgg::DiPhotonCandidate> dipho, int whichPho, float rho)
{
	if( whichPho > 1 ) {
		std::cout << "[bbggPlotter::getCHisoToCutValue] You chose the wrong photon!" << std::endl;
		return -1;
	}
	double PFIso = -1, eta = -99;
	if(whichPho == 0) {
		PFIso = dipho->leadingView().pfChIso03WrtChosenVtx();
		eta = dipho->leadingPhoton()->superCluster()->eta();
	}
	if(whichPho == 1) {
		PFIso = dipho->subLeadingView().pfChIso03WrtChosenVtx();
		eta = dipho->subLeadingPhoton()->superCluster()->eta();
	}
	
	double EA = bbggPlotter::getEA(eta, 0);
	double finalValue = fmax(PFIso - rho*EA, 0.);
	return finalValue;
}

double bbggPlotter::getNHisoToCutValue(const flashgg::Photon* pho, float rho)
{
	double PFIso = pho->pfNeutIso03();
	double eta = pho->superCluster()->eta();
	double EA = bbggPlotter::getEA(eta, 1);
	double extraFactor = 0;
	if(fabs(eta) < 1.479) extraFactor = exp(0.0028*pho->pt()+0.5408);
	if(fabs(eta) > 1.479) extraFactor = 0.01725*pho->pt();
	double finalValue = fmax(PFIso - rho*EA, 0.) - extraFactor;
	return finalValue;
}

double bbggPlotter::getPHisoToCutValue(const flashgg::Photon* pho, float rho)
{
	double PFIso = pho->pfPhoIso03();
	double eta = pho->superCluster()->eta();
	double EA = bbggPlotter::getEA(eta, 2);
	double extraFactor = 0;
	if(fabs(eta) < 1.479) extraFactor = 0.0014*pho->pt();
	if(fabs(eta) > 1.479) extraFactor = 0.0091*pho->pt();
	double finalValue = fmax(PFIso - rho*EA, 0.) - extraFactor;
	return finalValue;
}

double bbggPlotter::getEA( float eta, int whichEA){
        if(whichEA < 0 || whichEA > 2){
                std::cout << "WRONG EA TYPE" << std::endl;
                return -1;
        }

        float EA[7][3];

        EA[0][0] = 0.0234; EA[0][1] = 0.0053; EA[0][2] = 0.078;
        EA[1][0] = 0.0189; EA[1][1] = 0.0103; EA[1][2] = 0.0629;
        EA[2][0] = 0.0171; EA[2][1] = 0.0057; EA[2][2] = 0.0264;
        EA[3][0] = 0.0129; EA[3][1] = 0.0070; EA[3][2] = 0.0462;
        EA[4][0] = 0.0110; EA[4][1] = 0.0152; EA[4][2] = 0.0740;
        EA[5][0] = 0.0074; EA[5][1] = 0.0232; EA[5][2] = 0.0924;
        EA[6][0] = 0.0035; EA[6][1] = 0.1709; EA[6][2] = 0.1484;

        float feta = fabs(eta);

        if(feta > 0.000 && feta < 1.000 ) return EA[0][whichEA];
        if(feta > 1.000 && feta < 1.479 ) return EA[1][whichEA];
        if(feta > 1.479 && feta < 2.000 ) return EA[2][whichEA];
        if(feta > 2.000 && feta < 2.200 ) return EA[3][whichEA];
        if(feta > 2.200 && feta < 2.300 ) return EA[4][whichEA];
        if(feta > 2.300 && feta < 2.400 ) return EA[5][whichEA];
        if(feta > 2.400 && feta < 10.00 ) return EA[6][whichEA];

        return -1;
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
	hists["dijet_mass"] 	= TH1F("dijet_mass", "DiJet Mass; M(jj); Events", 100, 0, 500);
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

	hists2D["candmass_dijetmass"] = TH2F("candmass_dijetmass", "DiHiggs Candidate Mass vs DiJet Mass; M((jj#gamma#gamma) (GeV); M(jj)", 100, 100, 800, 100, 0, 500);
	hists2D["candmass_diphomass"] = TH2F("candmass_diphomass", "DiHiggs Candidate Mass vs DiPhoton Mass; M((jj#gamma#gamma) (GeV); M(#gamma#gamma)", 100, 100, 800, 100, 110, 150);
	hists2D["candmass_deltar"]    = TH2F("candmass_deltar", "DiHiggs Candidate Mass vs DeltaR between DiCandidates; M((jj#gamma#gamma) (GeV); #DeltaR(#gamma#gamma,jj)", 100, 100, 800, 100, -1, 10);
	hists2D["jet1bdis_jet2bdis"]    = TH2F("jet1bdis_jet2bdis", "Jet1 bDis vs Jet2 bDis; Jet1 bDis; Jet2 bDis", 100, 0., 1., 100, 0., 1.);
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
