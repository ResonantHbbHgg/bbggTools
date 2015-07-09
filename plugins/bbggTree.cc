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
	  LorentzVector dijetCandidate;

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
	edm::Ptr<flashgg::Jet> LeadingJet = tools_.GetSelected_leadingJetCandidate();
	edm::Ptr<flashgg::Jet> SubLeadingJet = tools_.GetSelected_subleadingJetCandidate();
	
    if(DEBUG) std::cout << "GOT TO THE END!!" << std::endl;
	diphotonCandidate = diphoCand->p4();
	leadingPhoton = diphoCand->leadingPhoton()->p4();
	subleadingPhoton = diphoCand->subLeadingPhoton()->p4();
	leadingJet = LeadingJet->p4();
	subleadingJet = SubLeadingJet->p4();
	leadingJet_bDis = LeadingJet->bDiscriminator(bTagType);
	subleadingJet_bDis = SubLeadingJet->bDiscriminator(bTagType);
	dijetCandidate = leadingJet + subleadingJet;
	dihiggsCandidate = diphotonCandidate + dijetCandidate;
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
	tree->Branch("dijetCandidate", &dijetCandidate);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
bbggTree::endJob() 
{
	outFile->cd();
	tree->Write();
	outFile->Close();
}


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
