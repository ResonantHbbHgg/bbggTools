// -*- C++ -*-
//
// Package:    flashgg/bbggTree
// Class:      bbggTree
// 
/**\class bbggTree bbggTree.cc flashgg/bbggTree/plugins/bbggTree.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rafael Teixeira De Lima
//         Created:  Mon, 06 Jul 2015 08:43:28 GMT
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
#include "TTree.h"
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

class bbggTree : public edm::EDAnalyzer {
   public:
      explicit bbggTree(const edm::ParameterSet&);
      ~bbggTree();

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
	  
	  //Tree objects
	  LorentzVector diphotonCandidate;
	  LorentzVector leadingPhoton;
	  LorentzVector subleadingPhoton;
	  LorentzVector dihiggsCandidate;
	  LorentzVector leadingJet;
	  double leadingJet_bDis;
	  LorentzVector subleadingJet;
	  double subleadingJet_bDis;

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
	  TTree* tree;
      std::string fileName;
	  
      //Event counter for cout's
      long unsigned int EvtCount;
};

//
// constructors and destructor
//
bbggTree::bbggTree(const edm::ParameterSet& iConfig) :
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


bbggTree::~bbggTree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
bbggTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    if( EvtCount%100 == 0 ) std::cout << "[bbggTree::analyze] Analyzing event number: " << EvtCount << std::endl;
    EvtCount++;
    using namespace edm;
	 
    double dipho_pt = -1, dipho_eta = -1, /*dipho_phi = -1,*/ dipho_mass = -1;
//    double dijet_pt = -1, dijet_eta = -1, dijet_phi = -1, dijet_mass = -1;
    double cand4_pt = -1, cand4_eta = -1, /*cand4_phi = -1,*/ cand4_mass = -1;
    double pho1_pt = -1, pho1_eta = -1/*, pho1_phi = -1*/;
    double pho1_hoe = -1, pho1_sieie = -1, /*pho1_r9 = -1,*/ pho1_chiso = -1, pho1_nhiso = -1, pho1_phiso = -1;
    int pho1_elveto = -1;
    double pho2_pt = -1, pho2_eta = -1/*, pho2_phi = -1*/;
    double pho2_hoe = -1, pho2_sieie = -1, /*pho2_r9 = -1,*/ pho2_chiso = -1, pho2_nhiso = -1, pho2_phiso = -1;
    int pho2_elveto = -1;
//    double jet1_pt = -1, jet1_eta = -1, jet1_phi = -1, jet1_bDis = -1, jet1_PUid = -1, jet1_drPho1 = -1, jet1_drPho2 = -1;
//    double jet2_pt = -1, jet2_eta = -1, jet2_phi = -1, jet2_bDis = -1, jet2_PUid = -1, jet2_drPho1 = -1, jet2_drPho2 = -1;
	
	bbggTree::LorentzVector _pho1, _pho2, _jet1, _jet2, _dipho, _dijet, _hhcand;

    Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
    iEvent.getByToken( diPhotonToken_, diPhotons );
    Handle<View<flashgg::Jet> > theJets;
    iEvent.getByToken( thejetToken_, theJets );
    Handle<double> rhoHandle;        // the old way for now...move to getbytoken?
    iEvent.getByLabel( rhoFixedGrid_, rhoHandle );
    const double rhoFixedGrd = *( rhoHandle.product() );
    tools_.setRho(rhoFixedGrd);

    bool isValidDiPhotonCandidate = false;

    edm::Ptr<reco::Vertex> CandVtx;
    edm::Ptr<flashgg::DiPhotonCandidate> diphoCand;

    //Begin DiPhoton Loop/Selection -----------------------------------------------------------
    for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ )
    {
 	if(diph_onlyfirst && diphoIndex > 0 ) break;

 	 edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( diphoIndex );
	 
 	 dipho_pt = dipho->pt(); dipho_eta = dipho->eta(); dipho_mass = dipho->mass();
		 
 	 if(dipho_mass < diph_mass[0] || dipho_mass > diph_mass[1]) continue;
 	 if(fabs(dipho_eta) > diph_eta[0] ) continue;
 	 if(dipho_pt < diph_pt[0] ) continue;
		 
 	 pho1_pt = dipho->leadingPhoton()->pt();			pho2_pt = dipho->subLeadingPhoton()->pt();
 	 pho1_eta = dipho->leadingPhoton()->superCluster()->eta();	pho2_eta = dipho->subLeadingPhoton()->superCluster()->eta();
// 	 pho1_phi = dipho->leadingPhoton()->superCluster()->phi();	pho2_phi = dipho->subLeadingPhoton()->superCluster()->phi();
 	 pho1_hoe = dipho->leadingPhoton()->hadronicOverEm(); 		pho2_hoe = dipho->subLeadingPhoton()->hadronicOverEm();
 	 pho1_sieie = dipho->leadingPhoton()->full5x5_sigmaIetaIeta();	pho2_sieie = dipho->subLeadingPhoton()->full5x5_sigmaIetaIeta();
// 	 pho1_r9 = dipho->leadingPhoton()->r9();			pho2_r9 = dipho->subLeadingPhoton()->r9();
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
	  _pho1 = dipho->leadingPhoton()->p4();
	  _pho2 = dipho->subLeadingPhoton()->p4();
	  _dipho = dipho->p4();
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
    bbggTree::LorentzVector DiJet(0,0,0,0);
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

 	  bbggTree::LorentzVector dijet = Jets[iJet]->p4() + Jets[jJet]->p4();
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
   
	_jet1 = jet1->p4();
	_jet2 = jet2->p4();
	_dijet = DiJet;
	
    //End Jets Loop/Selection -----------------------------------------------------------
   
    //Candidate assignment --------------------------------------------------------------
//    double deltaR = tools_.DeltaR(DiJet, diphoCand->p4());
//    double dijet_sumbdis = jet1->bDiscriminator(bTagType) + jet2->bDiscriminator(bTagType);
    bbggTree::LorentzVector HHCandidate = DiJet + diphoCand->p4();
    cand4_pt = HHCandidate.pt();
    cand4_eta = HHCandidate.eta();
    cand4_mass = HHCandidate.mass();
    if(DEBUG) std::cout << cand4_pt << "\t < \t" << cand_pt[0] << std::endl;
    if(cand4_pt < cand_pt[0] ) return;
    if(DEBUG) std::cout << fabs(cand4_eta) << "\t > \t" << cand_eta[0] << std::endl;
    if(fabs(cand4_eta) > cand_eta[0] ) return;
    if(DEBUG) std::cout << cand4_mass << "\t" << cand_mass[0] << "\t" << cand_mass[1] << std::endl;
    if(cand4_mass < cand_mass[0] || cand4_mass > cand_mass[1] ) return;
    if(DEBUG) std::cout << "Passed 4-candidate selection..." << std::endl;
	_hhcand = HHCandidate;
    //END Candidate assignment ---------------------------------------------------------- 

    if(DEBUG) std::cout << "GOT TO THE END!!" << std::endl;
	diphotonCandidate = _dipho;
	leadingPhoton = _pho1;
	subleadingPhoton = _pho2;
	dihiggsCandidate = _hhcand;
	leadingJet = _jet1;
	leadingJet_bDis = jet1->bDiscriminator(bTagType);
	subleadingJet = _jet2;
	subleadingJet_bDis = jet2->bDiscriminator(bTagType);
	tree->Fill();

    if(DEBUG) std::cout << "Histograms filled!" << std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void 
bbggTree::beginJob()
{
	outFile = new TFile(fileName.c_str(), "RECREATE");
    tree = new TTree("bbggTree", "Tree for HH->bbgg analyses");
	tree->Branch("diphotonCandidate", &diphotonCandidate);
	tree->Branch("leadingPhoton", &leadingPhoton);
	tree->Branch("subleadingPhoton", &subleadingPhoton);
	tree->Branch("dihiggsCandidate", &dihiggsCandidate);
	tree->Branch("leadingJet", &leadingJet);
	tree->Branch("leadingJet_bDis", &leadingJet_bDis, "leadingJet_bDis/D");
	tree->Branch("subleadingJet", &subleadingJet);
	tree->Branch("subleadingJet_bDis", &subleadingJet_bDis, "subleadingJet_bDis/D");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
bbggTree::endJob() 
{
	tree->Write();
	outFile->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
bbggTree::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
bbggTree::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
bbggTree::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
bbggTree::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
bbggTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(bbggTree);
