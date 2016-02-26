// -*- C++ -*-
//
// Package:    flashgg/bbggTruthDiff
// Class:      bbggTruthDiff
// 
/**\class bbggTruthDiff bbggTruthDiff.cc flashgg/bbggTruthDiff/plugins/bbggTruthDiff.cc

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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//FLASHgg files
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"

//Local
#include "flashgg/bbggTools/interface/bbggTools.h"
#include "flashgg/bbggTools/interface/bbggMC.h"

//
// class declaration
//

const int DEBUG = 0;

class bbggTruthDiff : public edm::EDAnalyzer {
   public:
      explicit bbggTruthDiff(const edm::ParameterSet&);
      ~bbggTruthDiff();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      typedef math::XYZTLorentzVector LorentzVector;
      typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
	  bbggTools tools_;
	  bbggMC	mcTools_;
	  
      //Parameter tokens
      edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diPhotonToken_;
      std::vector<edm::InputTag> inputTagJets_;
      edm::EDGetTokenT<edm::View<reco::GenJet> > genjetToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genToken_;
      edm::InputTag rhoFixedGrid_;
      std::string bTagType;
      double genTotalWeight;

      //Thresholds
      std::vector<double> phoIDlooseEB;
      std::vector<double> phoIDlooseEE;
      std::vector<double> phoIDmediumEB;
      std::vector<double> phoIDmediumEE;
      std::vector<double> phoIDtightEB;
      std::vector<double> phoIDtightEE;
      std::vector<double> phoISOlooseEB;
      std::vector<double> phoISOlooseEE;
      std::vector<double> phoISOmediumEB;
      std::vector<double> phoISOmediumEE;
      std::vector<double> phoISOtightEB;
      std::vector<double> phoISOtightEE;
      std::vector<double> nhCorrEB;
      std::vector<double> nhCorrEE;
      std::vector<double> phCorrEB;	  
      std::vector<double> phCorrEE;
	  
      std::vector<double> ph_pt;
      std::vector<double> ph_eta;
      std::vector<double> ph_r9;
      std::vector<int> ph_elVeto;
      std::vector<int> ph_doelVeto;
      std::vector<int> ph_doID;
      std::vector<std::string> ph_whichID;
      std::vector<int> ph_doISO;
      std::vector<std::string> ph_whichISO;
      std::vector<double> diph_pt;
      std::vector<double> diph_eta;
      std::vector<double> diph_mass;
      unsigned int diph_onlyfirst;
      std::vector<double> jt_pt;
      std::vector<double> jt_eta;
      std::vector<double> jt_drPho;
      std::vector<double> jt_bDis;
      std::vector<int> jt_doPU;
      std::vector<int> jt_doID;
      unsigned int n_bJets;
      std::vector<double> dijt_pt;
      std::vector<double> dijt_eta;
      std::vector<double> dijt_mass; 
      std::vector<double> cand_pt;
      std::vector<double> cand_eta;
      std::vector<double> cand_mass;
      std::vector<double> dr_cands;
      unsigned int nPromptPhotons;
      unsigned int doDoubleCountingMitigation;
	  
      //OutFile & Hists
      TFile* outFile;
      std::map<std::string, TH1F> hists;
      std::map<std::string, TH2F> hists2D;
      std::string fileName;

      //Event counter for cout's
      long unsigned int EvtCount;
};

bbggTruthDiff::bbggTruthDiff(const edm::ParameterSet& iConfig) :
diPhotonToken_( consumes<edm::View<flashgg::DiPhotonCandidate> >( iConfig.getUntrackedParameter<edm::InputTag> ( "DiPhotonTag", edm::InputTag( "flashggDiPhotons" ) ) ) ),
inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
genjetToken_( consumes<edm::View<reco::GenJet> >( iConfig.getUntrackedParameter<edm::InputTag>( "GenJetTag", edm::InputTag( "slimmedGenJets" ) ) ) ),
genToken_( consumes<edm::View<reco::GenParticle> >( iConfig.getUntrackedParameter<edm::InputTag>( "GenTag", edm::InputTag( "prunedGenParticles" ) ) ) )
{
   //now do what ever initialization is needed
	  tools_ = bbggTools();
	  mcTools_ = bbggMC();
	  
   	  EvtCount = 0;
//Default values for thresholds
      std::string def_bTagType;
      std::vector<double> def_ph_pt;
      std::vector<double> def_ph_eta;
      std::vector<double> def_ph_r9;
      std::vector<int> def_ph_elVeto;
      std::vector<int> def_ph_doelVeto;
      std::vector<int> def_ph_doID;
      std::vector<std::string> def_ph_whichID;
      std::vector<int> def_ph_doISO;
      std::vector<std::string> def_ph_whichISO;
      std::vector<double> def_diph_pt;
      std::vector<double> def_diph_eta;
      std::vector<double> def_diph_mass;
      unsigned int def_diph_onlyfirst;
      std::vector<double> def_jt_pt;
      std::vector<double> def_jt_eta;
      std::vector<double> def_jt_drPho;
      std::vector<double> def_jt_bDis;
      std::vector<int> def_jt_doPU;
      std::vector<int> def_jt_doID;
      unsigned int def_n_bJets;
      std::vector<double> def_dijt_pt;
      std::vector<double> def_dijt_eta;
      std::vector<double> def_dijt_mass;
      std::vector<double> def_cand_pt;
      std::vector<double> def_cand_eta;
      std::vector<double> def_cand_mass;
      std::vector<double> def_dr_cands;
      def_ph_pt.push_back(10.);           def_ph_pt.push_back(10.);
      def_ph_eta.push_back(20.);          def_ph_eta.push_back(20.);
      def_ph_r9.push_back(-1.);           def_ph_r9.push_back(-1.);
      def_ph_elVeto.push_back(-1);       def_ph_elVeto.push_back(-1);
      def_ph_doelVeto.push_back(0);       def_ph_elVeto.push_back(0);
      def_ph_doID.push_back(0);		    def_ph_doID.push_back(0);
      def_ph_whichID.push_back("loose");  def_ph_whichID.push_back("loose");
      def_ph_doISO.push_back(0);	        def_ph_doISO.push_back(0);
      def_ph_whichISO.push_back("loose"); def_ph_whichISO.push_back("loose");
      def_diph_pt.push_back(10.);         def_diph_pt.push_back(10.);
      def_diph_eta.push_back(0.);         def_diph_eta.push_back(0.);
      def_diph_mass.push_back(0.);        def_diph_mass.push_back(1000.);
      def_diph_onlyfirst = 0;
      def_jt_pt.push_back(10.);           def_jt_pt.push_back(10.);
      def_jt_eta.push_back(20.);          def_jt_eta.push_back(20.);
      def_jt_bDis.push_back(0.);          def_jt_bDis.push_back(0.);
      def_jt_doPU.push_back(0);           def_jt_doPU.push_back(0);
      def_jt_doID.push_back(0);           def_jt_doID.push_back(0);
      def_n_bJets = 0;
      def_dijt_pt.push_back(10.);         def_dijt_pt.push_back(10.);
      def_dijt_eta.push_back(20.);        def_dijt_eta.push_back(20.);
      def_dijt_mass.push_back(0.);        def_dijt_mass.push_back(1000.);
      def_jt_drPho.push_back(0.5);
      def_cand_pt.push_back(0.);
      def_cand_eta.push_back(20.);
      def_cand_mass.push_back(0.);		def_cand_mass.push_back(2000.);
	 
      def_dr_cands.push_back(0.11);
	  
      unsigned int def_nPromptPhotons = 0;
      unsigned int def_doDoubleCountingMitigation = 0;
	  

      std::string def_fileName;

      def_bTagType = "pfCombinedInclusiveSecondaryVertexV2BJetTags";
      def_fileName =  "out.root";

      //Get thresholds from config file
      //Get thresholds from config file
      phoIDlooseEB = iConfig.getUntrackedParameter<std::vector<double > >("phoIDlooseEB");
      phoIDlooseEE = iConfig.getUntrackedParameter<std::vector<double > >("phoIDlooseEE");
      phoIDmediumEB = iConfig.getUntrackedParameter<std::vector<double > >("phoIDmediumEB");
      phoIDmediumEE = iConfig.getUntrackedParameter<std::vector<double > >("phoIDmediumEE");
      phoIDtightEB = iConfig.getUntrackedParameter<std::vector<double > >("phoIDtightEB");
      phoIDtightEE = iConfig.getUntrackedParameter<std::vector<double > >("phoIDtightEE");
      phoISOlooseEB = iConfig.getUntrackedParameter<std::vector<double > >("phoISOlooseEB");
      phoISOlooseEE = iConfig.getUntrackedParameter<std::vector<double > >("phoISOlooseEE");
      phoISOmediumEB = iConfig.getUntrackedParameter<std::vector<double > >("phoISOmediumEB");
      phoISOmediumEE = iConfig.getUntrackedParameter<std::vector<double > >("phoISOmediumEE");
      phoISOtightEB = iConfig.getUntrackedParameter<std::vector<double > >("phoISOtightEB");
      phoISOtightEE = iConfig.getUntrackedParameter<std::vector<double > >("phoISOtightEE");
      nhCorrEB = iConfig.getUntrackedParameter<std::vector<double > >("nhCorrEB");
      nhCorrEE = iConfig.getUntrackedParameter<std::vector<double > >("nhCorrEE");
      phCorrEB = iConfig.getUntrackedParameter<std::vector<double > >("phCorrEB");
      phCorrEE = iConfig.getUntrackedParameter<std::vector<double > >("phCorrEE");
    
      ph_pt     = iConfig.getUntrackedParameter<std::vector<double > >("PhotonPtOverDiPhotonMass", def_ph_pt);
      ph_eta    = iConfig.getUntrackedParameter<std::vector<double > >("PhotonEta", def_ph_eta);
      ph_r9     = iConfig.getUntrackedParameter<std::vector<double > >("PhotonR9", def_ph_r9);
      ph_elVeto = iConfig.getUntrackedParameter<std::vector<int > >("PhotonElectronVeto", def_ph_elVeto);
      ph_doelVeto = iConfig.getUntrackedParameter<std::vector<int > >("PhotonDoElectronVeto", def_ph_doelVeto);
      ph_doID   = iConfig.getUntrackedParameter<std::vector<int > >("PhotonDoID", def_ph_doID);
      ph_whichID   = iConfig.getUntrackedParameter<std::vector<std::string> >("PhotonWhichID", def_ph_whichID);
      ph_doISO  = iConfig.getUntrackedParameter<std::vector<int > >("PhotonDoISO", def_ph_doISO);
      ph_whichISO  = iConfig.getUntrackedParameter<std::vector<std::string > >("PhotonWhichISO", def_ph_whichISO);
      diph_pt   = iConfig.getUntrackedParameter<std::vector<double > >("DiPhotonPt", def_diph_pt);
      diph_eta  = iConfig.getUntrackedParameter<std::vector<double > >("DiPhotonEta", def_diph_eta);
      diph_mass = iConfig.getUntrackedParameter<std::vector<double > >("DiPhotonMassWindow", def_diph_mass);
      diph_onlyfirst = iConfig.getUntrackedParameter<unsigned int>("DiPhotonOnlyFirst", def_diph_onlyfirst);
      jt_pt     = iConfig.getUntrackedParameter<std::vector<double > >("JetPtOverDiJetMass", def_jt_pt);
      jt_eta    = iConfig.getUntrackedParameter<std::vector<double > >("JetEta", def_jt_eta);
      jt_drPho  = iConfig.getUntrackedParameter<std::vector<double > >("JetDrPho", def_jt_drPho);
      jt_bDis   = iConfig.getUntrackedParameter<std::vector<double > >("JetBDiscriminant", def_jt_bDis);
      jt_doPU   = iConfig.getUntrackedParameter<std::vector<int > >("JetDoPUID", def_jt_doPU);
      jt_doID   = iConfig.getUntrackedParameter<std::vector<int > >("JetDoID", def_jt_doID);
      n_bJets = iConfig.getUntrackedParameter<unsigned int>("n_bJets", def_n_bJets);
      dijt_pt   = iConfig.getUntrackedParameter<std::vector<double > >("DiJetPt", def_dijt_pt);
      dijt_eta  = iConfig.getUntrackedParameter<std::vector<double > >("DiJetEta", def_dijt_eta);
      dijt_mass = iConfig.getUntrackedParameter<std::vector<double > >("DiJetMassWindow", def_dijt_mass);
      cand_pt 	= iConfig.getUntrackedParameter<std::vector<double > >("CandidatePt", def_cand_mass);
      cand_eta 	= iConfig.getUntrackedParameter<std::vector<double > >("CandidateEta", def_cand_mass);
      cand_mass = iConfig.getUntrackedParameter<std::vector<double > >("CandidateMassWindow", def_cand_mass);
      dr_cands  = iConfig.getUntrackedParameter<std::vector<double > >("CandidatesDeltaR", def_dr_cands);
	  
      nPromptPhotons = iConfig.getUntrackedParameter<unsigned int>("nPromptPhotons", def_nPromptPhotons);
      doDoubleCountingMitigation = iConfig.getUntrackedParameter<unsigned int>("doDoubleCountingMitigation", def_doDoubleCountingMitigation);

      rhoFixedGrid_  = iConfig.getUntrackedParameter<edm::InputTag>( "rhoFixedGridCollection", edm::InputTag( "fixedGridRhoAll" ) );

      bTagType = iConfig.getUntrackedParameter<std::string>( "bTagType", def_bTagType );

      fileName = iConfig.getUntrackedParameter<std::string>( "OutFileName", def_fileName );
	
      tools_.SetCut_PhotonPtOverDiPhotonMass( ph_pt );
      tools_.SetCut_PhotonEta( ph_eta );
      tools_.SetCut_PhotonDoID( ph_doID );
      tools_.SetCut_PhotonDoISO( ph_doISO );
      tools_.SetCut_PhotonR9( ph_r9 );
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
      tools_.SetCut_JetDoID( jt_doID );
      tools_.SetCut_n_bJets( n_bJets );
      tools_.SetCut_DiJetPt( dijt_pt );
      tools_.SetCut_DiJetEta( dijt_eta );
      tools_.SetCut_DiJetMassWindow( dijt_mass );
      tools_.SetCut_CandidateMassWindow( cand_mass );
      tools_.SetCut_CandidatePt( cand_pt );
      tools_.SetCut_CandidateEta( cand_eta );
      tools_.SetCut_bTagType( bTagType );
      tools_.SetCut_CandidatesDeltaR( dr_cands );
    
      std::map<int, vector<double> > phoIDloose;
      phoIDloose[0] = phoIDlooseEB;
      phoIDloose[1] = phoIDlooseEE;
      tools_.SetCut_phoIDloose(phoIDloose);
    
      std::map<int, vector<double> > phoIDmedium;
      phoIDmedium[0] = phoIDmediumEB;
      phoIDmedium[1] = phoIDmediumEE;
      tools_.SetCut_phoIDmedium(phoIDmedium);
    
      std::map<int, vector<double> > phoIDtight;
      phoIDtight[0] = phoIDtightEB;
      phoIDtight[1] = phoIDtightEE;
      tools_.SetCut_phoIDtight(phoIDtight);
    
      std::map<int, vector<double> > phoISOloose;
      phoISOloose[0] = phoISOlooseEB;
      phoISOloose[1] = phoISOlooseEE;
      tools_.SetCut_phoISOloose(phoISOloose);
    
      std::map<int, vector<double> > phoISOmedium;
      phoISOmedium[0] = phoISOmediumEB;
      phoISOmedium[1] = phoISOmediumEE;
      tools_.SetCut_phoISOmedium(phoISOmedium);
    
      std::map<int, vector<double> > phoISOtight;
      phoISOtight[0] = phoISOtightEB;
      phoISOtight[1] = phoISOtightEE;
      tools_.SetCut_phoISOtight(phoISOtight);
    
      std::map<int, vector<double> > nhCorr;
      nhCorr[0] = nhCorrEB;
      nhCorr[1] = nhCorrEE;
      tools_.SetCut_nhCorr(nhCorr);
    
      std::map<int, vector<double> > phCorr;
      phCorr[0] = phCorrEB;
      phCorr[1] = phCorrEE;
      tools_.SetCut_phCorr(phCorr);
 
      tools_.SetCut_phoWhichID(ph_whichID);
      tools_.SetCut_phoWhichISO(ph_whichISO);
	  
	  
      std::cout << "Parameters initialized... cand_mass[0]: " << cand_mass[0] << "\t" << "cand_mass[1] " << cand_mass[1] << std::endl;

}


bbggTruthDiff::~bbggTruthDiff()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
bbggTruthDiff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   if( EvtCount%100 == 0 ) std::cout << "[bbggTruthDiff::analyze] Analyzing event number: " << EvtCount << std::endl;
   EvtCount++;
   using namespace edm;
	             
   Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
   iEvent.getByToken( diPhotonToken_, diPhotons );
   
   //Get Jets collections!
   JetCollectionVector theJetsCols( inputTagJets_.size() );
   for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
       iEvent.getByLabel( inputTagJets_[j], theJetsCols[j] );
   }

   
   Handle<View<reco::GenJet> > genJets;
   iEvent.getByToken( genjetToken_, genJets );
   
   Handle<double> rhoHandle;        // the old way for now...move to getbytoken?
   iEvent.getByLabel( rhoFixedGrid_, rhoHandle );
   const double rhoFixedGrd = *( rhoHandle.product() );
   tools_.setRho(rhoFixedGrd);
   
   Handle<View<reco::GenParticle> > genParticles;
   iEvent.getByToken( genToken_, genParticles);
   
//PreLoop
	vector<edm::Ptr<flashgg::DiPhotonCandidate>> diphoVec;
	for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ )
	{
		edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( diphoIndex );
		if(doDoubleCountingMitigation){
			bbggMC _mcTools = bbggMC();
			unsigned int nPrompt = _mcTools.CheckNumberOfPromptPhotons(dipho, genParticles);
			if (nPrompt == nPromptPhotons) diphoVec.push_back(dipho);
			else continue;
		} else {
			diphoVec.push_back(dipho);
		}
	}
   
   
//	bool passedSelection = tools_.AnalysisSelection(diPhotons, theJets);
	bool passedSelection = tools_.AnalysisSelection(diphoVec, theJetsCols);
	if(!passedSelection) return;
	
	edm::Ptr<flashgg::DiPhotonCandidate> diphoCand = tools_.GetSelected_diphoCandidate();
    unsigned int jetCollectionIndex = diphoCand->jetCollectionIndex();
    
    
	bool passedMCSelection = mcTools_.MatchTruth(diPhotons, theJetsCols[jetCollectionIndex], genParticles);
	if(!passedMCSelection) return;

	edm::Ptr<flashgg::DiPhotonCandidate> diphoCandMC = mcTools_.GetSelected_diphoCandidate();

	//Photon Variables
	        double dipho_pt         = diphoCand->pt()      ;
        double dipho_eta        = diphoCand->eta()     ;
        double dipho_phi        = diphoCand->phi()     ;
        double dipho_mass       = diphoCand->mass();
        double pho1_pt 					= diphoCand->leadingPhoton()->pt();
        double pho2_pt 					=                diphoCand->subLeadingPhoton()->pt();
        double pho1_eta 				=               diphoCand->leadingPhoton()->superCluster()->eta();
        double pho2_eta 				=               diphoCand->subLeadingPhoton()->superCluster()->eta();
        double pho1_phi 				=               diphoCand->leadingPhoton()->superCluster()->phi();
        double pho2_phi 				=               diphoCand->subLeadingPhoton()->superCluster()->phi();
        double pho1_hoe 				=               diphoCand->leadingPhoton()->hadronicOverEm();
        double pho2_hoe 				=               diphoCand->subLeadingPhoton()->hadronicOverEm();
        double pho1_sieie 			=     diphoCand->leadingPhoton()->full5x5_sigmaIetaIeta();
        double pho2_sieie 			=     diphoCand->subLeadingPhoton()->full5x5_sigmaIetaIeta();
        double pho1_r9 					=                diphoCand->leadingPhoton()->r9();
        double pho2_r9 					=                diphoCand->subLeadingPhoton()->r9();
        double pho1_elveto 			=    diphoCand->leadingPhoton()->passElectronVeto();
        double pho2_elveto 			=    diphoCand->subLeadingPhoton()->passElectronVeto();
        // double pho1_chiso             = tools_.getCHisoToCutValue( diphoCand, 0);
        // double pho2_chiso             = tools_.getCHisoToCutValue( diphoCand, 1);
        // double pho1_nhiso             = tools_.getNHisoToCutValue( diphoCand->leadingPhoton() );
        // double pho2_nhiso             = tools_.getNHisoToCutValue( diphoCand->subLeadingPhoton() );
        // double pho1_phiso             = tools_.getPHisoToCutValue( diphoCand->leadingPhoton() );
        // double pho2_phiso             = tools_.getPHisoToCutValue( diphoCand->subLeadingPhoton() );

 double MCdipho_pt         =  diphoCandMC->pt();
 double MCdipho_eta        =  diphoCandMC->eta();
 double MCdipho_phi        =  diphoCandMC->phi();
 double MCdipho_mass       =  diphoCandMC->mass();
 double MCpho1_pt 					=  diphoCandMC->leadingPhoton()->pt();
 double MCpho2_pt 					=  diphoCandMC->subLeadingPhoton()->pt();
 double MCpho1_eta 				=  diphoCandMC->leadingPhoton()->superCluster()->eta();
 double MCpho2_eta 				=  diphoCandMC->subLeadingPhoton()->superCluster()->eta();
 double MCpho1_phi 				=  diphoCandMC->leadingPhoton()->superCluster()->phi();
 double MCpho2_phi 				=  diphoCandMC->subLeadingPhoton()->superCluster()->phi();
 double MCpho1_hoe 				=  diphoCandMC->leadingPhoton()->hadronicOverEm();
 double MCpho2_hoe 				=  diphoCandMC->subLeadingPhoton()->hadronicOverEm();
 double MCpho1_sieie 			=  diphoCandMC->leadingPhoton()->full5x5_sigmaIetaIeta();
 double MCpho2_sieie 			=  diphoCandMC->subLeadingPhoton()->full5x5_sigmaIetaIeta();
 double MCpho1_r9 					=  diphoCandMC->leadingPhoton()->r9();
 double MCpho2_r9 					=  diphoCandMC->subLeadingPhoton()->r9();
 double MCpho1_elveto 			=  diphoCandMC->leadingPhoton()->passElectronVeto();
 double MCpho2_elveto 			=  diphoCandMC->subLeadingPhoton()->passElectronVeto();
 // double MCpho1_chiso             =  tools_.getCHisoToCutValue( diphoCandMC, 0);
 // double MCpho2_chiso             =  tools_.getCHisoToCutValue( diphoCandMC, 1);
 // double MCpho1_nhiso             =  tools_.getNHisoToCutValue( diphoCandMC->leadingPhoton() );
 // double MCpho2_nhiso             =  tools_.getNHisoToCutValue( diphoCandMC->subLeadingPhoton() );
 // double MCpho1_phiso             =  tools_.getPHisoToCutValue( diphoCandMC->leadingPhoton() );
 // double MCpho2_phiso             =  tools_.getPHisoToCutValue( diphoCandMC->subLeadingPhoton() );

 	//END photon Variables
	
	edm::Ptr<flashgg::Jet> LeadingJet = tools_.GetSelected_leadingJetCandidate();
	edm::Ptr<flashgg::Jet> LeadingJetMC = mcTools_.GetSelected_leadingJetCandidate();
double jet1_pt =                LeadingJet->pt();
double jet1_eta =               LeadingJet->eta();
double jet1_phi =               LeadingJet->phi();
double jet1_bDis =              LeadingJet->bDiscriminator(bTagType);
double jet1_PUid =              LeadingJet->passesPuJetId(diphoCand->vtx());
double jet1_drPho1 =    tools_.DeltaR(LeadingJet->p4(), diphoCand->leadingPhoton()->p4() );
double jet1_drPho2 =    tools_.DeltaR(LeadingJet->p4(), diphoCand->subLeadingPhoton()->p4() );
double MCjet1_pt =     LeadingJetMC->pt();
double MCjet1_eta =    LeadingJetMC->eta();
double MCjet1_phi =    LeadingJetMC->phi();
double MCjet1_bDis =   LeadingJetMC->bDiscriminator(bTagType);
double MCjet1_PUid =   LeadingJetMC->passesPuJetId(diphoCand->vtx());
double MCjet1_drPho1 = tools_.DeltaR(LeadingJetMC->p4(), diphoCand->leadingPhoton()->p4() );
double MCjet1_drPho2 = tools_.DeltaR(LeadingJetMC->p4(), diphoCand->subLeadingPhoton()->p4() );

	//Leading Jet Variables
//	double jet1_pt = 		LeadingJet->pt() - LeadingJetMC->pt();
//	double jet1_eta = 		LeadingJet->eta() - LeadingJetMC->eta();
//	double jet1_phi = 		LeadingJet->phi() - LeadingJetMC->phi();
//	double jet1_bDis = 		LeadingJet->bDiscriminator(bTagType) - LeadingJetMC->bDiscriminator(bTagType);
//	double jet1_PUid = 		LeadingJet->passesPuJetId(diphoCand->vtx()) - LeadingJetMC->passesPuJetId(diphoCand->vtx());
//	double jet1_drPho1 = 	tools_.DeltaR(LeadingJet->p4(), diphoCand->leadingPhoton()->p4() ) - tools_.DeltaR(LeadingJetMC->p4(), diphoCand->leadingPhoton()->p4() );
//	double jet1_drPho2 = 	tools_.DeltaR(LeadingJet->p4(), diphoCand->subLeadingPhoton()->p4() ) - tools_.DeltaR(LeadingJetMC->p4(), diphoCand->subLeadingPhoton()->p4() );
	//END Leading Jet Variables
	
	edm::Ptr<flashgg::Jet> SubLeadingJet = tools_.GetSelected_subleadingJetCandidate();
	edm::Ptr<flashgg::Jet> SubLeadingJetMC = mcTools_.GetSelected_subleadingJetCandidate();
	//SubLeading Jet Variables
double jet2_pt =                SubLeadingJet->pt();
double jet2_eta =               SubLeadingJet->eta();
double jet2_phi =               SubLeadingJet->phi();
double jet2_bDis =              SubLeadingJet->bDiscriminator(bTagType);
double jet2_PUid =              SubLeadingJet->passesPuJetId(diphoCand->vtx());
double jet2_drPho1 =    tools_.DeltaR(SubLeadingJet->p4(), diphoCand->leadingPhoton()->p4() );
double jet2_drPho2 =    tools_.DeltaR(SubLeadingJet->p4(), diphoCand->subLeadingPhoton()->p4() );
double MCjet2_pt =     SubLeadingJetMC->pt();
double MCjet2_eta =    SubLeadingJetMC->eta();
double MCjet2_phi =    SubLeadingJetMC->phi();
double MCjet2_bDis =   SubLeadingJetMC->bDiscriminator(bTagType);
double MCjet2_PUid =   SubLeadingJetMC->passesPuJetId(diphoCand->vtx());
double MCjet2_drPho1 = tools_.DeltaR(SubLeadingJetMC->p4(), diphoCand->leadingPhoton()->p4() );
double MCjet2_drPho2 = tools_.DeltaR(SubLeadingJetMC->p4(), diphoCand->subLeadingPhoton()->p4() );

//	double jet2_pt = 		SubLeadingJet->pt() - SubLeadingJetMC->pt();
//	double jet2_eta = 		SubLeadingJet->eta() - SubLeadingJetMC->eta();
//	double jet2_phi = 		SubLeadingJet->phi() - SubLeadingJetMC->phi();
//	double jet2_bDis = 		SubLeadingJet->bDiscriminator(bTagType) - SubLeadingJetMC->bDiscriminator(bTagType);
//	double jet2_PUid = 		SubLeadingJet->passesPuJetId(diphoCand->vtx()) - SubLeadingJetMC->passesPuJetId(diphoCand->vtx());
//	double jet2_drPho1 = 	tools_.DeltaR(SubLeadingJet->p4(), diphoCand->leadingPhoton()->p4() ) - tools_.DeltaR(SubLeadingJetMC->p4(), diphoCand->leadingPhoton()->p4() );
//	double jet2_drPho2 = 	tools_.DeltaR(SubLeadingJet->p4(), diphoCand->subLeadingPhoton()->p4() ) - tools_.DeltaR(SubLeadingJetMC->p4(), diphoCand->subLeadingPhoton()->p4() );
	//END SubLeading Jet Variables
	
	//DiJet Variables
	LorentzVector DiJet = LeadingJet->p4() + SubLeadingJet->p4();
	LorentzVector DiJetMC = LeadingJetMC->p4() + SubLeadingJetMC->p4();
double dijet_pt =               DiJet.pt();
double dijet_eta =              DiJet.eta();
double dijet_phi =              DiJet.phi();
double dijet_mass =     DiJet.mass();
double deltaR =                 tools_.DeltaR(DiJet, diphoCand->p4());
double dijet_sumbdis =  LeadingJet->bDiscriminator(bTagType) + SubLeadingJet->bDiscriminator(bTagType);
double MCdijet_sumbdis = LeadingJetMC->bDiscriminator(bTagType) + SubLeadingJetMC->bDiscriminator(bTagType);
double MCdijet_pt =   DiJetMC.pt();
double MCdijet_eta =  DiJetMC.eta();
double MCdijet_phi =  DiJetMC.phi();
double MCdijet_mass = DiJetMC.mass();
double MCdeltaR =     tools_.DeltaR(DiJetMC, diphoCand->p4());

//    	double dijet_pt = 		DiJet.pt() - DiJetMC.pt();
 //	double dijet_eta = 		DiJet.eta() - DiJetMC.eta();
 //	double dijet_phi = 		DiJet.phi() - DiJetMC.phi();
 //	double dijet_mass = 	DiJet.mass() - DiJetMC.mass();
   // 	double deltaR = 		tools_.DeltaR(DiJet, diphoCand->p4()) - tools_.DeltaR(DiJetMC, diphoCand->p4());
    //	double dijet_sumbdis = 	LeadingJet->bDiscriminator(bTagType) + SubLeadingJet->bDiscriminator(bTagType) - LeadingJetMC->bDiscriminator(bTagType) - SubLeadingJetMC->bDiscriminator(bTagType);
	//END DiJet Variables
	//HH Candidate Variables
	LorentzVector HHCandidate = DiJet + diphoCand->p4();
	LorentzVector HHCandidateMC = DiJetMC + diphoCandMC->p4();
double cand4_pt =       HHCandidate.pt();
double cand4_eta =      HHCandidate.eta();
double cand4_phi =      HHCandidate.phi();
double cand4_mass = HHCandidate.mass();
double MCcand4_pt =   HHCandidateMC.pt();
double MCcand4_eta =  HHCandidateMC.eta();
double MCcand4_phi =  HHCandidateMC.phi();
double MCcand4_mass = HHCandidateMC.mass();

//    	double cand4_pt = 	HHCandidate.pt() - HHCandidateMC.pt();
 //	double cand4_eta = 	HHCandidate.eta() - HHCandidateMC.eta();
 //	double cand4_phi = 	HHCandidate.phi() - HHCandidateMC.phi();
 //	double cand4_mass = HHCandidate.mass() - HHCandidateMC.mass();
	//END HH Candidate Variables
	//Other variables
	double candmass_minDeltaR_lJP = (jet1_drPho1 < jet1_drPho2) ? jet1_drPho1 : jet1_drPho2;
	double candmass_minDeltaR_sJP = (jet2_drPho1 < jet2_drPho2) ? jet2_drPho1 : jet2_drPho2;

	int nbjets = -1;
	
   if(DEBUG) std::cout << "GOT TO THE END!!" << std::endl;

   hists2D["jet1bdis_jet2bdis"].Fill(jet1_bDis, jet2_bDis);
   hists2D["candmass_deltar"].Fill( cand4_mass, deltaR );
   hists2D["candmass_dijetmass"].Fill(cand4_mass, dijet_mass);
   hists2D["candmass_diphomass"].Fill(cand4_mass, dipho_mass);
   hists2D["candmass_minDeltaR_lJP"].Fill(cand4_mass, candmass_minDeltaR_lJP);
   hists2D["candmass_minDeltaR_sJP"].Fill(cand4_mass, candmass_minDeltaR_sJP);
   hists2D["diphomass_lpho"].Fill(dipho_mass, pho1_pt);
   hists2D["dijetmass_deltar"].Fill(dijet_mass, deltaR);
   hists2D["dijetmass_minDeltaR_lJP"].Fill(dijet_mass, candmass_minDeltaR_lJP);
   hists2D["dijetmass_minDeltaR_sJP"].Fill(dijet_mass, candmass_minDeltaR_lJP);
   hists2D["dijetmass_lJpt"].Fill(dijet_mass, jet1_pt);
   hists2D["dijetmass_sJpt"].Fill(dijet_mass, jet2_pt);
   hists2D["dijetmass_sumbdis"].Fill(dijet_mass, dijet_sumbdis);
   hists2D["dijetmass_ratio_l"].Fill(dijet_mass, jet1_pt/dijet_mass);
   hists2D["dijetmass_ratio_s"].Fill(dijet_mass, jet2_pt/dijet_mass);
   hists["n_bjets"].Fill(nbjets);
   hists2D["dr_dijetdipho"].Fill( deltaR , MCdeltaR );
   hists2D["dipho_pt"].Fill(dipho_pt, MCdipho_pt);
   hists2D["dipho_eta"].Fill(dipho_eta, MCdipho_eta);
   hists2D["dipho_phi"].Fill(dipho_phi, MCdipho_phi);
   hists2D["dipho_mass"].Fill(dipho_mass, MCdipho_mass);
   hists2D["dijet_pt"].Fill(dijet_pt, MCdijet_pt);
   hists2D["dijet_eta"].Fill(dijet_eta, MCdijet_eta);
   hists2D["dijet_phi"].Fill(dijet_phi, MCdijet_phi);
   hists2D["dijet_mass"].Fill(dijet_mass, MCdijet_mass);
   hists2D["dijet_sumbdis"].Fill(dijet_sumbdis, MCdijet_sumbdis);
   hists2D["cand4_pt"].Fill(cand4_pt, MCcand4_pt);
   hists2D["cand4_eta"].Fill(cand4_eta, MCcand4_eta);
   hists2D["cand4_phi"].Fill(cand4_phi, MCcand4_phi);
   hists2D["cand4_mass"].Fill(cand4_mass, MCcand4_mass);
   hists2D["pho1_pt"].Fill(pho1_pt, MCpho1_pt);
   hists2D["pho1_eta"].Fill(pho1_eta, MCpho1_eta);
   hists2D["pho1_phi"].Fill(pho1_phi, MCpho1_phi);
   hists2D["pho1_hoe"].Fill(pho1_hoe, MCpho1_hoe);
   hists2D["pho1_sieie"].Fill(pho1_sieie, MCpho1_sieie);
   hists2D["pho1_r9"].Fill(pho1_r9, MCpho1_r9);
   // hists2D["pho1_chiso"].Fill(pho1_chiso, MCpho1_chiso);
   // hists2D["pho1_nhiso"].Fill(pho1_nhiso, MCpho1_nhiso);
   // hists2D["pho1_phiso"].Fill(pho1_phiso, MCpho1_phiso);
   hists2D["pho1_elveto"].Fill(pho1_elveto, MCpho1_elveto);
   hists2D["pho2_pt"].Fill(pho2_pt, MCpho2_pt);
   hists2D["pho2_eta"].Fill(pho2_eta, MCpho2_eta);
   hists2D["pho2_phi"].Fill(pho2_phi, MCpho2_phi);
   hists2D["pho2_hoe"].Fill(pho2_hoe, MCpho2_hoe);
   hists2D["pho2_sieie"].Fill(pho2_sieie, MCpho2_sieie);
   hists2D["pho2_r9"].Fill(pho2_r9, MCpho2_r9);
   // hists2D["pho2_chiso"].Fill(pho2_chiso, MCpho2_chiso);
   // hists2D["pho2_nhiso"].Fill(pho2_nhiso, MCpho2_nhiso);
   // hists2D["pho2_phiso"].Fill(pho2_phiso, MCpho2_phiso);
   hists2D["pho2_elveto"].Fill(pho2_elveto, MCpho2_elveto);
   hists2D["jet1_pt"].Fill(jet1_pt, MCjet1_pt);
   hists2D["jet1_eta"].Fill(jet1_eta, MCjet1_eta);
   hists2D["jet1_phi"].Fill(jet1_phi, MCjet1_phi);
   hists2D["jet1_drPho1"].Fill(jet1_drPho1, MCjet1_drPho1);
   hists2D["jet1_drPho2"].Fill(jet1_drPho2, MCjet1_drPho2);
   hists2D["jet1_bDis"].Fill(jet1_bDis, MCjet1_bDis);
   hists2D["jet1_PUid"].Fill(jet1_PUid, MCjet1_PUid);
   hists2D["jet2_pt"].Fill(jet2_pt, MCjet2_pt);
   hists2D["jet2_eta"].Fill(jet2_eta, MCjet2_eta);
   hists2D["jet2_phi"].Fill(jet2_phi, MCjet2_phi);
   hists2D["jet2_drPho1"].Fill(jet2_drPho1, MCjet2_drPho1);
   hists2D["jet2_drPho2"].Fill(jet2_drPho2, MCjet2_drPho2);
   hists2D["jet2_bDis"].Fill(jet2_bDis, MCjet2_bDis);
   hists2D["jet2_PUid"].Fill(jet2_PUid, MCjet2_PUid);

/*


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
*/

   if(DEBUG) std::cout << "Histograms filled!" << std::endl;
}

// ------------ method called once each job just before starting event loop  ------------
void 
bbggTruthDiff::beginJob()
{
	outFile = new TFile(fileName.c_str(), "RECREATE");

         hists2D["dipho_pt"]    = TH2F("dipho_pt", "DiPhoton p_{T}; p_{T}(#gamma#gamma) (GeV); MCMacthed", 100, 0, 400, 100, 0, 400);
         hists2D["dipho_eta"]   = TH2F("dipho_eta", "DiPhoton #eta; #eta(#gamma#gamma); MCMacthed", 100, -5, 5, 100, -5, 5);
         hists2D["dipho_phi"]   = TH2F("dipho_phi", "DiPhoton #phi; #phi(#gamma#gamma); MCMacthed", 100, -3.5, 3.5, 100, -3.5, 3.5);
         hists2D["dipho_mass"]  = TH2F("dipho_mass", "DiPhoton Mass; M(#gamma#gamma); MCMacthed", 	100, 100, 150, 100, 100, 150);
         hists2D["dijet_pt"]    = TH2F("dijet_pt", "DiJet p_{T}; p_{T}(jj) (GeV); MCMacthed", 			100, 0, 500, 100, 0, 500);
         hists2D["dijet_eta"]   = TH2F("dijet_eta", "DiJet #eta; #eta(jj); MCMacthed", 		100, -5, 5, 100, -5, 5);
         hists2D["dijet_phi"]   = TH2F("dijet_phi", "DiJet #phi; #phi(jj); MCMacthed", 		100, -3.5, 3.5, 100, -3.5, 3.5);
         hists2D["dijet_mass"]  = TH2F("dijet_mass", "DiJet Mass; M(jj); MCMacthed", 			100, 30, 300, 100, 30, 300);
         hists2D["dijet_sumbdis"]  = TH2F("dijet_sumbdis", "Sum of b-Discriminant of DiJet; Sum b-Dis; MCMacthed", 100, -0.1, 2.1, 100, -0.1, 2.1);
         hists2D["dr_dijetdipho"]  = TH2F("dr_dijetdipho", "DeltaR between DiJet and DiPhoton; #DeltaR(#gamma#gamma,jj); MCMacthed", 			100, -1, 10, 100, -1, 10);
         hists2D["cand4_pt"]    = TH2F("cand_pt", "DiHiggs Candidate (jj#gamma#gamma) p_{T}; p_{T}(jj#gamma#gamma) (GeV); MCMacthed", 			100, 0, 700, 100, 0, 700);
         hists2D["cand4_eta"]   = TH2F("cand_eta", "DiHiggs Candidate (jj#gamma#gamma) #eta; #eta(jj#gamma#gamma); MCMacthed", 		100, -5, 5, 100, -5, 5);
         hists2D["cand4_phi"]   = TH2F("cand_phi", "DiHiggs Candidate (jj#gamma#gamma) #phi; #phi(jj#gamma#gamma); MCMacthed", 		100, -3.5, 3.5, 100, -3.5, 3.5);
         hists2D["cand4_mass"]  = TH2F("cand_mass", "DiHiggs Candidate (jj#gamma#gamma) Mass; M(jj#gamma#gamma) (GeV); MCMacthed", 100, 100, 1000, 100, 100, 1000);
         hists2D["pho1_pt"]     = TH2F("pho1_pt", "Leading Photon p_{T}; p_{T}(leading #gamma) (GeV); MCMacthed", 	100, 10, 400, 100, 10, 400);
         hists2D["pho1_eta"]    = TH2F("pho1_eta", "Leading Photon #eta; #eta(leading #gamma); MCMacthed", 100, -5., 5., 100, -5., 5.);
         hists2D["pho1_phi"]    = TH2F("pho1_phi", "Leading Photon #phi; #phi(leading #gamma); MCMacthed", 100, -3.5, 3.5, 100, -3.5, 3.5);
         hists2D["pho1_hoe"]    = TH2F("pho1_hoe", "Leading Photon H/E; H/E(leading #gamma); MCMacthed", 	100, -0.01, 0.1, 100, -0.01, 0.1);
         hists2D["pho1_sieie"]  = TH2F("pho1_sieie", "Leading Photon #sigma_{i#etai#eta}; #sigma_{i#etai#eta}(leading #gamma); MCMacthed", 100, 0.0, 0.04, 100, 0.0, 0.04);
         hists2D["pho1_r9"]     = TH2F("pho1_r9", "Leading Photon R9; R9(leading #gamma); MCMacthed", 			100, 0, 1.1, 100, 0, 1.1);
         hists2D["pho1_chiso"]  = TH2F("pho1_chiso", "Leading Photon Corrected Charged Isolation; Corrected Charged Isolation(leading #gamma); MCMacthed", 100, -1., 5, 100, -1., 5);
         hists2D["pho1_nhiso"]  = TH2F("pho1_nhiso", "Leading Photon Corrected Neutral Isolation; Corrected Neutral Isolation(leading #gamma); MCMacthed", 100, -1., 5, 100, -1., 5);
         hists2D["pho1_phiso"]  = TH2F("pho1_phiso", "Leading Photon Corrected Photon Isolation; Corrected Photon Isolation(leading #gamma); MCMacthed", 	100, -1., 5, 100, -1., 5);
         hists2D["pho1_elveto"] = TH2F("pho1_elveto", "Leading Photon Electron Veto; Electron Veto(leading #gamma); events", 		8, -1, 3, 8, -1, 3);
         hists2D["pho2_pt"]     = TH2F("pho2_pt", "SubLeading Photon p_{T}; p_{T}(subLeading #gamma) (GeV); MCMacthed", 		100, 10, 200, 100, 10, 200);
         hists2D["pho2_eta"]    = TH2F("pho2_eta", "SubLeading Photon #eta; #eta(subLeading #gamma); MCMacthed", 	100, -5., 5., 100, -5., 5.);
         hists2D["pho2_phi"]    = TH2F("pho2_phi", "SubLeading Photon #phi; #phi(subleading #gamma); MCMacthed", 	100, -3.5, 3.5, 100, -3.5, 3.5);
         hists2D["pho2_hoe"]    = TH2F("pho2_hoe", "SubLeading Photon H/E; H/E(subLeading #gamma); MCMacthed", 		100, -0.01, 0.1, 100, -0.01, 0.1);
         hists2D["pho2_sieie"]  = TH2F("pho2_sieie", "SubLeading Photon #sigma_{i#etai#eta}; #sigma_{i#etai#eta}(subLeading #gamma); MCMacthed", 	100, 0.0, 0.04, 100, 0.0, 0.04);
         hists2D["pho2_r9"]     = TH2F("pho2_r9", "SubLeading Photon R9; R9(subLeading #gamma); MCMacthed", 100, 0, 1.1, 100, 0, 1.1);
         hists2D["pho2_chiso"]  = TH2F("pho2_chiso", "SubLeading Photon Corrected Charged Isolation; Corrected Charged Isolation(subLeading #gamma); MCMacthed", 	100, -1., 5, 100, -1., 5);
         hists2D["pho2_nhiso"]  = TH2F("pho2_nhiso", "SubLeading Photon Corrected Neutral Isolation; Corrected Neutral Isolation(subLeading #gamma); MCMacthed", 	100, -1., 5, 100, -1., 5);
         hists2D["pho2_phiso"]  = TH2F("pho2_phiso", "SubLeading Photon Corrected Photon Isolation; Corrected Photon Isolation(subLeading #gamma); MCMacthed", 		100, -1., 5, 100, -1., 5);
         hists2D["pho2_elveto"] = TH2F("pho2_elveto", "SubLeading Photon Electron Veto; Electron Veto(subLeading #gamma); events", 			8, -1, 3, 8, -1, 3);
         hists2D["jet1_pt"]     = TH2F("jet1_pt", "Leading Jet p_{T}; p_{T}(leading jet) (GeV); MCMacthed", 100, 0, 300, 100, 0, 300);
         hists2D["jet1_eta"]    = TH2F("jet1_eta", "Leading Jet #eta; #eta(leading jet); MCMacthed", 			100, -5., 5., 100, -5., 5.);
         hists2D["jet1_phi"]    = TH2F("jet1_phi", "Leading Jet #phi; #phi(leading jet); MCMacthed", 			100, -3.5, 3.5, 100, -3.5, 3.5);
         hists2D["jet1_drPho1"] = TH2F("jet1_drPho1", "Leading Jet DeltaR with Leading Photon; DeltaR(leading jet, leading #gamma); MCMacthed", 		100, -1, 10, 100, -1, 10);
         hists2D["jet1_drPho2"] = TH2F("jet1_drPho2", "Leading Jet DeltaR with SubLeading Photon; DeltaR(leading jet, subleading #gamma); MCMacthed", 			100, -1, 10, 100, -1, 10);
         hists2D["jet1_bDis"]   = TH2F("jet1_bDis", "Leading Jet b-Discriminant; b-Discriminant(leading jet); MCMacthed", 	100, -0.01, 1.01, 100, -0.01, 1.01);
         hists2D["jet1_PUid"]   = TH2F("jet1_PUid", "Leading Jet PU ID; PU ID(leading jet); MCMacthed", 		8, -1, 3, 8, -1, 3);
         hists2D["jet2_pt"]     = TH2F("jet2_pt", "SubLeading Jet p_{T}; p_{T}(subLeading jet) (GeV); MCMacthed", 	100, 0, 300, 100, 0, 300);
         hists2D["jet2_eta"]    = TH2F("jet2_eta", "SubLeading Jet #eta; #eta(subLeading jet); MCMacthed", 100, -5., 5., 100, -5., 5.);
         hists2D["jet2_phi"]    = TH2F("jet2_phi", "SubLeading Jet #phi; #phi(subleading jet); MCMacthed", 100, -3.5, 3.5, 100, -3.5, 3.5);
         hists2D["jet2_drPho1"] = TH2F("jet2_drPho1", "SubLeading Jet DeltaR with Leading Photon; DeltaR(subleading jet, leading #gamma); MCMacthed", 			100, -1, 10, 100, -1, 10);
         hists2D["jet2_drPho2"] = TH2F("jet2_drPho2", "SubLeading Jet DeltaR with SubLeading Photon; DeltaR(subleadign jet, subleading #gamma); MCMacthed", 100, -1, 10, 100, -1, 10);
         hists2D["jet2_bDis"]   = TH2F("jet2_bDis", "SubLeading Jet b-Discriminant; b-Discriminant(subLeading jet); MCMacthed", 		100, -0.01, 1.01, 100, -0.01, 1.01);
         hists2D["jet2_PUid"]   = TH2F("jet2_PUid", "SubLeading Jet PU ID; PU ID(subleading jet); MCMacthed", 			8, -1, 3, 8, -1, 3);

//	hists["dipho_pt"] 	= TH1F("dipho_pt", "DiPhoton p_{T}; p_{T}(#gamma#gamma) (GeV); Events", 100, -50, 50); // 0, 400);
//	hists["dipho_eta"] 	= TH1F("dipho_eta", "DiPhoton #eta; #eta(#gamma#gamma); Events", 100, -50, 50); // -5, 5);
//	hists["dipho_phi"]	= TH1F("dipho_phi", "DiPhoton #phi; #phi(#gamma#gamma); Events", 100, -50, 50); // -3.5, 3.5);
//	hists["dipho_mass"] 	= TH1F("dipho_mass", "DiPhoton Mass; M(#gamma#gamma); Events", 100, -50, 50); // 100, 150);

//	hists["dijet_pt"] 	= TH1F("dijet_pt", "DiJet p_{T}; p_{T}(jj) (GeV); Events", 100, -50, 50); // 0, 500);
//	hists["dijet_eta"] 	= TH1F("dijet_eta", "DiJet #eta; #eta(jj); Events", 100, -50, 50); // -5, 5);
//	hists["dijet_phi"]	= TH1F("dijet_phi", "DiJet #phi; #phi(jj); Events", 100, -50, 50); // -3.5, 3.5);
//	hists["dijet_mass"] 	= TH1F("dijet_mass", "DiJet Mass; M(jj); Events", 100, -500, 500); // 30, 300);
//	hists["dijet_sumbdis"]  = TH1F("dijet_sumbdis", "Sum of b-Discriminant of DiJet; Sum b-Dis; Events", 100, -50, 50); // -0.1, 2.1); 

//	hists["dr_dijetdipho"]	= TH1F("dr_dijetdipho", "DeltaR between DiJet and DiPhoton; #DeltaR(#gamma#gamma,jj); Events", 100, -50, 50); // -1, 10);

//	hists["cand4_pt"] 	= TH1F("cand_pt", "DiHiggs Candidate (jj#gamma#gamma) p_{T}; p_{T}(jj#gamma#gamma) (GeV); Events", 100, -50, 50); // 0, 700);
//	hists["cand4_eta"] 	= TH1F("cand_eta", "DiHiggs Candidate (jj#gamma#gamma) #eta; #eta(jj#gamma#gamma); Events", 100, -10, 10); // -5, 5);
//	hists["cand4_phi"]	= TH1F("cand_phi", "DiHiggs Candidate (jj#gamma#gamma) #phi; #phi(jj#gamma#gamma); Events", 100, -50, 50); // -3.5, 3.5); 
//	hists["cand4_mass"] 	= TH1F("cand_mass", "DiHiggs Candidate (jj#gamma#gamma) Mass; M(jj#gamma#gamma) (GeV); Events", 100, -500, 500); // 100, 1000);

//	hists["pho1_pt"] 	= TH1F("pho1_pt", "Leading Photon p_{T}; p_{T}(leading #gamma) (GeV); Events", 100, -50, 50); // 10, 400);
//	hists["pho1_eta"] 	= TH1F("pho1_eta", "Leading Photon #eta; #eta(leading #gamma); Events", 100, -50, 50); // -5., 5.);
//	hists["pho1_phi"]	= TH1F("pho1_phi", "Leading Photon #phi; #phi(leading #gamma); Events", 100, -50, 50); // -3.5, 3.5);

//	hists["pho1_hoe"] 	= TH1F("pho1_hoe", "Leading Photon H/E; H/E(leading #gamma); Events", 100, -50, 50); // -0.01, 0.1);
//	hists["pho1_sieie"] 	= TH1F("pho1_sieie", "Leading Photon #sigma_{i#etai#eta}; #sigma_{i#etai#eta}(leading #gamma); Events", 100, -50, 50); // 0.0, 0.04);
//	hists["pho1_r9"] 	= TH1F("pho1_r9", "Leading Photon R9; R9(leading #gamma); Events", 100, -50, 50); // 0, 1.1);
//	hists["pho1_chiso"] 	= TH1F("pho1_chiso", "Leading Photon Corrected Charged Isolation; Corrected Charged Isolation(leading #gamma); Events", 100, -50, 50); // -1., 5);
//	hists["pho1_nhiso"] 	= TH1F("pho1_nhiso", "Leading Photon Corrected Neutral Isolation; Corrected Neutral Isolation(leading #gamma); Events", 100, -50, 50); // -1., 5);
//	hists["pho1_phiso"] 	= TH1F("pho1_phiso", "Leading Photon Corrected Photon Isolation; Corrected Photon Isolation(leading #gamma); Events", 100, -50, 50); // -1., 5);
//	hists["pho1_elveto"] 	= TH1F("pho1_elveto", "Leading Photon Electron Veto; Electron Veto(leading #gamma); events", 8, -1, 3);

//	hists["pho2_pt"] 	= TH1F("pho2_pt", "SubLeading Photon p_{T}; p_{T}(subLeading #gamma) (GeV); Events", 100, -50, 50); // 10, 200);
//	hists["pho2_eta"] 	= TH1F("pho2_eta", "SubLeading Photon #eta; #eta(subLeading #gamma); Events", 100, -50, 50); // -5., 5.);
//	hists["pho2_phi"]	= TH1F("pho2_phi", "SubLeading Photon #phi; #phi(subleading #gamma); Events", 100, -50, 50); // -3.5, 3.5);

//	hists["pho2_hoe"] 	= TH1F("pho2_hoe", "SubLeading Photon H/E; H/E(subLeading #gamma); Events", 100, -50, 50); // -0.01, 0.1);
//	hists["pho2_sieie"] 	= TH1F("pho2_sieie", "SubLeading Photon #sigma_{i#etai#eta}; #sigma_{i#etai#eta}(subLeading #gamma); Events", 100, -50, 50); // 0.0, 0.04);
//	hists["pho2_r9"] 	= TH1F("pho2_r9", "SubLeading Photon R9; R9(subLeading #gamma); Events", 100, -50, 50); // 0, 1.1);
//	hists["pho2_chiso"] 	= TH1F("pho2_chiso", "SubLeading Photon Corrected Charged Isolation; Corrected Charged Isolation(subLeading #gamma); Events", 100, -50, 50); // -1., 5);
//	hists["pho2_nhiso"] 	= TH1F("pho2_nhiso", "SubLeading Photon Corrected Neutral Isolation; Corrected Neutral Isolation(subLeading #gamma); Events", 100, -50, 50); // -1., 5);
//	hists["pho2_phiso"] 	= TH1F("pho2_phiso", "SubLeading Photon Corrected Photon Isolation; Corrected Photon Isolation(subLeading #gamma); Events", 100, -50, 50); // -1., 5);
//	hists["pho2_elveto"] 	= TH1F("pho2_elveto", "SubLeading Photon Electron Veto; Electron Veto(subLeading #gamma); events", 8, -1, 3);


//	hists["jet1_pt"] 	= TH1F("jet1_pt", "Leading Jet p_{T}; p_{T}(leading jet) (GeV); Events", 100, -50, 50); // 0, 300);
//	hists["jet1_eta"] 	= TH1F("jet1_eta", "Leading Jet #eta; #eta(leading jet); Events", 100, -50, 50); // -5., 5.);
//	hists["jet1_phi"]	= TH1F("jet1_phi", "Leading Jet #phi; #phi(leading jet); Events", 100, -50, 50); // -3.5, 3.5);
//	hists["jet1_drPho1"]	= TH1F("jet1_drPho1", "Leading Jet DeltaR with Leading Photon; DeltaR(leading jet, leading #gamma); Events", 100, -50, 50); // -1, 10);
//	hists["jet1_drPho2"]	= TH1F("jet1_drPho2", "Leading Jet DeltaR with SubLeading Photon; DeltaR(leading jet, subleading #gamma); Events", 100, -50, 50); // -1, 10);
//	hists["jet1_bDis"] 	= TH1F("jet1_bDis", "Leading Jet b-Discriminant; b-Discriminant(leading jet); Events", 100, -50, 50); // -0.01, 1.01);
//	hists["jet1_PUid"] 	= TH1F("jet1_PUid", "Leading Jet PU ID; PU ID(leading jet); Events", 8, -1, 3);

//	hists["jet2_pt"] 	= TH1F("jet2_pt", "SubLeading Jet p_{T}; p_{T}(subLeading jet) (GeV); Events", 100, -50, 50); // 0, 300);
//	hists["jet2_eta"] 	= TH1F("jet2_eta", "SubLeading Jet #eta; #eta(subLeading jet); Events", 100, -50, 50); // -5., 5.);
//	hists["jet2_phi"]	= TH1F("jet2_phi", "SubLeading Jet #phi; #phi(subleading jet); Events", 100, -50, 50); // -3.5, 3.5);
//	hists["jet2_drPho1"]	= TH1F("jet2_drPho1", "SubLeading Jet DeltaR with Leading Photon; DeltaR(subleading jet, leading #gamma); Events", 100, -50, 50); // -1, 10);
//	hists["jet2_drPho2"]	= TH1F("jet2_drPho2", "SubLeading Jet DeltaR with SubLeading Photon; DeltaR(subleadign jet, subleading #gamma); Events", 100, -50, 50); // -1, 10);
//	hists["jet2_bDis"] 	= TH1F("jet2_bDis", "SubLeading Jet b-Discriminant; b-Discriminant(subLeading jet); Events", 100, -50, 50); // -0.01, 1.01);
//	hists["jet2_PUid"] 	= TH1F("jet2_PUid", "SubLeading Jet PU ID; PU ID(subleading jet); Events", 8, -1, 3);


	// hists["dipho_pt"] 	= TH1F("dipho_pt", "DiPhoton p_{T}; p_{T}(#gamma#gamma) (GeV); Events", 100, 0, 400);
	// hists["dipho_eta"] 	= TH1F("dipho_eta", "DiPhoton #eta; #eta(#gamma#gamma); Events", 100, -5, 5);
	// hists["dipho_phi"]	= TH1F("dipho_phi", "DiPhoton #phi; #phi(#gamma#gamma); Events", 100, -3.5, 3.5);
	// hists["dipho_mass"] 	= TH1F("dipho_mass", "DiPhoton Mass; M(#gamma#gamma); Events", 100, 100, 150);
	//
	// hists["dijet_pt"] 	= TH1F("dijet_pt", "DiJet p_{T}; p_{T}(jj) (GeV); Events", 100, 0, 500);
	// hists["dijet_eta"] 	= TH1F("dijet_eta", "DiJet #eta; #eta(jj); Events", 100, -5, 5);
	// hists["dijet_phi"]	= TH1F("dijet_phi", "DiJet #phi; #phi(jj); Events", 100, -3.5, 3.5);
	// hists["dijet_mass"] 	= TH1F("dijet_mass", "DiJet Mass; M(jj); Events", 100, 30, 300);
	// hists["dijet_sumbdis"]  = TH1F("dijet_sumbdis", "Sum of b-Discriminant of DiJet; Sum b-Dis; Events", 100, -0.1, 2.1);
	//
	// hists["dr_dijetdipho"]	= TH1F("dr_dijetdipho", "DeltaR between DiJet and DiPhoton; #DeltaR(#gamma#gamma,jj); Events", 100, -1, 10);
	//
	// hists["cand4_pt"] 	= TH1F("cand_pt", "DiHiggs Candidate (jj#gamma#gamma) p_{T}; p_{T}(jj#gamma#gamma) (GeV); Events", 100, 0, 700);
	// hists["cand4_eta"] 	= TH1F("cand_eta", "DiHiggs Candidate (jj#gamma#gamma) #eta; #eta(jj#gamma#gamma); Events", 100, -5, 5);
	// hists["cand4_phi"]	= TH1F("cand_phi", "DiHiggs Candidate (jj#gamma#gamma) #phi; #phi(jj#gamma#gamma); Events", 100, -3.5, 3.5);
	// hists["cand4_mass"] 	= TH1F("cand_mass", "DiHiggs Candidate (jj#gamma#gamma) Mass; M(jj#gamma#gamma) (GeV); Events", 100, 100, 1000);
	//
	// hists["pho1_pt"] 	= TH1F("pho1_pt", "Leading Photon p_{T}; p_{T}(leading #gamma) (GeV); Events", 100, 10, 400);
	// hists["pho1_eta"] 	= TH1F("pho1_eta", "Leading Photon #eta; #eta(leading #gamma); Events", 100, -5., 5.);
	// hists["pho1_phi"]	= TH1F("pho1_phi", "Leading Photon #phi; #phi(leading #gamma); Events", 100, -3.5, 3.5);
	//
	// hists["pho1_hoe"] 	= TH1F("pho1_hoe", "Leading Photon H/E; H/E(leading #gamma); Events", 100, -0.01, 0.1);
	// hists["pho1_sieie"] 	= TH1F("pho1_sieie", "Leading Photon #sigma_{i#etai#eta}; #sigma_{i#etai#eta}(leading #gamma); Events", 100, 0.0, 0.04);
	// hists["pho1_r9"] 	= TH1F("pho1_r9", "Leading Photon R9; R9(leading #gamma); Events", 100, 0, 1.1);
	// hists["pho1_chiso"] 	= TH1F("pho1_chiso", "Leading Photon Corrected Charged Isolation; Corrected Charged Isolation(leading #gamma); Events", 100, -1., 5);
	// hists["pho1_nhiso"] 	= TH1F("pho1_nhiso", "Leading Photon Corrected Neutral Isolation; Corrected Neutral Isolation(leading #gamma); Events", 100, -1., 5);
	// hists["pho1_phiso"] 	= TH1F("pho1_phiso", "Leading Photon Corrected Photon Isolation; Corrected Photon Isolation(leading #gamma); Events", 100, -1., 5);
	// hists["pho1_elveto"] 	= TH1F("pho1_elveto", "Leading Photon Electron Veto; Electron Veto(leading #gamma); events", 8, -1, 3);
	//
	// hists["pho2_pt"] 	= TH1F("pho2_pt", "SubLeading Photon p_{T}; p_{T}(subLeading #gamma) (GeV); Events", 100, 10, 200);
	// hists["pho2_eta"] 	= TH1F("pho2_eta", "SubLeading Photon #eta; #eta(subLeading #gamma); Events", 100, -5., 5.);
	// hists["pho2_phi"]	= TH1F("pho2_phi", "SubLeading Photon #phi; #phi(subleading #gamma); Events", 100, -3.5, 3.5);
	//
	// hists["pho2_hoe"] 	= TH1F("pho2_hoe", "SubLeading Photon H/E; H/E(subLeading #gamma); Events", 100, -0.01, 0.1);
	// hists["pho2_sieie"] 	= TH1F("pho2_sieie", "SubLeading Photon #sigma_{i#etai#eta}; #sigma_{i#etai#eta}(subLeading #gamma); Events", 100, 0.0, 0.04);
	// hists["pho2_r9"] 	= TH1F("pho2_r9", "SubLeading Photon R9; R9(subLeading #gamma); Events", 100, 0, 1.1);
	// hists["pho2_chiso"] 	= TH1F("pho2_chiso", "SubLeading Photon Corrected Charged Isolation; Corrected Charged Isolation(subLeading #gamma); Events", 100, -1., 5);
	// hists["pho2_nhiso"] 	= TH1F("pho2_nhiso", "SubLeading Photon Corrected Neutral Isolation; Corrected Neutral Isolation(subLeading #gamma); Events", 100, -1., 5);
	// hists["pho2_phiso"] 	= TH1F("pho2_phiso", "SubLeading Photon Corrected Photon Isolation; Corrected Photon Isolation(subLeading #gamma); Events", 100, -1., 5);
	// hists["pho2_elveto"] 	= TH1F("pho2_elveto", "SubLeading Photon Electron Veto; Electron Veto(subLeading #gamma); events", 8, -1, 3);
	//
	//
	// hists["jet1_pt"] 	= TH1F("jet1_pt", "Leading Jet p_{T}; p_{T}(leading jet) (GeV); Events", 100, 0, 300);
	// hists["jet1_eta"] 	= TH1F("jet1_eta", "Leading Jet #eta; #eta(leading jet); Events", 100, -5., 5.);
	// hists["jet1_phi"]	= TH1F("jet1_phi", "Leading Jet #phi; #phi(leading jet); Events", 100, -3.5, 3.5);
	// hists["jet1_drPho1"]	= TH1F("jet1_drPho1", "Leading Jet DeltaR with Leading Photon; DeltaR(leading jet, leading #gamma); Events", 100, -1, 10);
	// hists["jet1_drPho2"]	= TH1F("jet1_drPho2", "Leading Jet DeltaR with SubLeading Photon; DeltaR(leading jet, subleading #gamma); Events", 100, -1, 10);
	// hists["jet1_bDis"] 	= TH1F("jet1_bDis", "Leading Jet b-Discriminant; b-Discriminant(leading jet); Events", 100, -0.01, 1.01);
	// hists["jet1_PUid"] 	= TH1F("jet1_PUid", "Leading Jet PU ID; PU ID(leading jet); Events", 8, -1, 3);
	//
	// hists["jet2_pt"] 	= TH1F("jet2_pt", "SubLeading Jet p_{T}; p_{T}(subLeading jet) (GeV); Events", 100, 0, 300);
	// hists["jet2_eta"] 	= TH1F("jet2_eta", "SubLeading Jet #eta; #eta(subLeading jet); Events", 100, -5., 5.);
	// hists["jet2_phi"]	= TH1F("jet2_phi", "SubLeading Jet #phi; #phi(subleading jet); Events", 100, -3.5, 3.5);
	// hists["jet2_drPho1"]	= TH1F("jet2_drPho1", "SubLeading Jet DeltaR with Leading Photon; DeltaR(subleading jet, leading #gamma); Events", 100, -1, 10);
	// hists["jet2_drPho2"]	= TH1F("jet2_drPho2", "SubLeading Jet DeltaR with SubLeading Photon; DeltaR(subleadign jet, subleading #gamma); Events", 100, -1, 10);
	// hists["jet2_bDis"] 	= TH1F("jet2_bDis", "SubLeading Jet b-Discriminant; b-Discriminant(subLeading jet); Events", 100, -0.01, 1.01);
	// hists["jet2_PUid"] 	= TH1F("jet2_PUid", "SubLeading Jet PU ID; PU ID(subleading jet); Events", 8, -1, 3);

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
bbggTruthDiff::endJob() 
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
bbggTruthDiff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(bbggTruthDiff);
