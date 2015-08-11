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
	  LorentzVector leadingPhoton;
	  vector<int> leadingPhotonID; 
	  vector<int> leadingPhotonISO;
	  LorentzVector subleadingPhoton;
	  vector<int> subleadingPhotonID;
	  vector<int> subleadingPhotonISO;
	  vector<LorentzVector> jets;
	  vector<float> jets_bDiscriminant;
	  vector<int> jets_PUid;

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
      std::string def_bTagType;

      std::string def_fileName;

      def_bTagType = "pfCombinedInclusiveSecondaryVertexV2BJetTags";
      def_fileName =  "out.root";

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

      rhoFixedGrid_  = iConfig.getUntrackedParameter<edm::InputTag>( "rhoFixedGridCollection", edm::InputTag( "fixedGridRhoAll" ) );

      bTagType = iConfig.getUntrackedParameter<std::string>( "bTagType", def_bTagType );

      fileName = iConfig.getUntrackedParameter<std::string>( "OutFileName", def_fileName );

      std::cout << "Parameters initialized... fileName: " << fileName <<  std::endl;

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
	leadingPhotonID.clear();
	leadingPhotonISO.clear();
	subleadingPhotonID.clear();
	subleadingPhotonISO.clear();
	jets.clear();
	jets_bDiscriminant.clear();
	jets_PUid.clear();
	 
    Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
    iEvent.getByToken( diPhotonToken_, diPhotons );
    Handle<View<flashgg::Jet> > theJets;
    iEvent.getByToken( thejetToken_, theJets );
    Handle<double> rhoHandle;        // the old way for now...move to getbytoken?
    iEvent.getByLabel( rhoFixedGrid_, rhoHandle );
    const double rhoFixedGrd = *( rhoHandle.product() );
    tools_.setRho(rhoFixedGrd);

	if(diPhotons->size() < 1) return;
	if(theJets->size() < 2) return;
	double dipho_pt_ref = 0;
	edm::Ptr<flashgg::DiPhotonCandidate> diphoCand;
    for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ )
	{
    	edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( diphoIndex );
		if(dipho->pt() > dipho_pt_ref){
			diphoCand = dipho;
			dipho_pt_ref = dipho->pt();
		}
	}

	LorentzVector lPho = diphoCand->leadingPhoton()->p4();
	leadingPhoton = lPho;//LorentzVector(lPho.px(), lPho.py(), lPho.pz(), lPho.energy() );	
	LorentzVector sPho = diphoCand->subLeadingPhoton()->p4();
	subleadingPhoton = sPho;//LorentzVector(sPho.px(), sPho.py(), sPho.pz(), sPho.energy() );

	int lphoIDloose = 0;
	int lphoIDmedium = 0;
	int lphoIDtight = 0;
	int sphoIDloose = 0;
	int sphoIDmedium = 0;
	int sphoIDtight = 0;
	int lphoISOloose = 0;
	int lphoISOmedium = 0;
	int lphoISOtight = 0;
	int sphoISOloose = 0;
	int sphoISOmedium = 0;
	int sphoISOtight = 0;

	float lpho_eta_abs = fabs( diphoCand->leadingPhoton()->superCluster()->eta() );
	float spho_eta_abs = fabs( diphoCand->subLeadingPhoton()->superCluster()->eta() );

	if(lpho_eta_abs < 1.44){
		lphoIDloose 	= tools_.isPhoID(diphoCand->leadingPhoton(), phoIDlooseEB);
		lphoIDmedium 	= tools_.isPhoID(diphoCand->leadingPhoton(), phoIDmediumEB);
		lphoIDtight 	= tools_.isPhoID(diphoCand->leadingPhoton(), phoIDtightEB);
		lphoISOloose 	= tools_.isPhoISO(diphoCand, 0, phoISOlooseEB,  nhCorrEB, phCorrEB);
		lphoISOmedium 	= tools_.isPhoISO(diphoCand, 0, phoISOmediumEB, nhCorrEB, phCorrEB);
		lphoISOtight 	= tools_.isPhoISO(diphoCand, 0, phoISOtightEB,  nhCorrEB, phCorrEB);
	}
	if(lpho_eta_abs > 1.44){
		lphoIDloose 	= tools_.isPhoID(diphoCand->leadingPhoton(), phoIDlooseEE);
		lphoIDmedium 	= tools_.isPhoID(diphoCand->leadingPhoton(), phoIDmediumEE);
		lphoIDtight 	= tools_.isPhoID(diphoCand->leadingPhoton(), phoIDtightEE);
		lphoISOloose 	= tools_.isPhoISO(diphoCand, 0, phoISOlooseEE,  nhCorrEE, phCorrEE);
		lphoISOmedium 	= tools_.isPhoISO(diphoCand, 0, phoISOmediumEE, nhCorrEE, phCorrEE);
		lphoISOtight 	= tools_.isPhoISO(diphoCand, 0, phoISOtightEE,  nhCorrEE, phCorrEE);
	}
	
	if(spho_eta_abs < 1.44){
		sphoIDloose 	= tools_.isPhoID(diphoCand->subLeadingPhoton(), phoIDlooseEB);
		sphoIDmedium 	= tools_.isPhoID(diphoCand->subLeadingPhoton(), phoIDmediumEB);
		sphoIDtight 	= tools_.isPhoID(diphoCand->subLeadingPhoton(), phoIDtightEB);
		sphoISOloose 	= tools_.isPhoISO(diphoCand, 1, phoISOlooseEB,  nhCorrEB, phCorrEB);
		sphoISOmedium 	= tools_.isPhoISO(diphoCand, 1, phoISOmediumEB, nhCorrEB, phCorrEB);
		sphoISOtight 	= tools_.isPhoISO(diphoCand, 1, phoISOtightEB,  nhCorrEB, phCorrEB);
	}
	if(spho_eta_abs > 1.44){
		sphoIDloose 	= tools_.isPhoID(diphoCand->subLeadingPhoton(), phoIDlooseEE);
		sphoIDmedium 	= tools_.isPhoID(diphoCand->subLeadingPhoton(), phoIDmediumEE);
		sphoIDtight 	= tools_.isPhoID(diphoCand->subLeadingPhoton(), phoIDtightEE);
		sphoISOloose 	= tools_.isPhoISO(diphoCand, 1, phoISOlooseEE,  nhCorrEE, phCorrEE);
		sphoISOmedium 	= tools_.isPhoISO(diphoCand, 1, phoISOmediumEE, nhCorrEE, phCorrEE);
		sphoISOtight 	= tools_.isPhoISO(diphoCand, 1, phoISOtightEE,  nhCorrEE, phCorrEE);
	}
	
	leadingPhotonID.push_back(lphoIDloose);
	leadingPhotonID.push_back(lphoIDmedium);
	leadingPhotonID.push_back(lphoIDtight);
	leadingPhotonISO.push_back(lphoISOloose);
	leadingPhotonISO.push_back(lphoISOmedium);
	leadingPhotonISO.push_back(lphoISOtight);
	
	subleadingPhotonID.push_back(sphoIDloose);
	subleadingPhotonID.push_back(sphoIDmedium);
	subleadingPhotonID.push_back(sphoIDtight);
	subleadingPhotonISO.push_back(sphoISOloose);
	subleadingPhotonISO.push_back(sphoISOmedium);
	subleadingPhotonISO.push_back(sphoISOtight);

    for( unsigned int jetIndex = 0; jetIndex < theJets->size(); jetIndex++ )
	{
		edm::Ptr<flashgg::Jet> jet = theJets->ptrAt( jetIndex );
		LorentzVector tjet = jet->p4(); //LorentzVector(jet->px(), jet->py(), jet->pz(), jet->energy());
		jets.push_back(tjet);
		int jetpuid = jet->passesPuJetId(diphoCand->vtx());
		jets_PUid.push_back(jetpuid);
		float jetbdis = jet->bDiscriminator(bTagType);
		jets_bDiscriminant.push_back(jetbdis);
	}

    if(DEBUG) std::cout << "GOT TO THE END!!" << std::endl;
	tree->Fill();

    if(DEBUG) std::cout << "Histograms filled!" << std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void 
bbggTree::beginJob()
{
	outFile = new TFile(fileName.c_str(), "RECREATE");
        tree = new TTree("bbggTree", "Flat tree for HH->bbgg analyses");
	tree->Branch("leadingPhoton", &leadingPhoton);
	tree->Branch("leadingPhotonID", &leadingPhotonID);
	tree->Branch("leadingPhotonISO", &leadingPhotonISO);
	tree->Branch("subleadingPhoton", &subleadingPhoton);
	tree->Branch("subleadingPhotonID", &subleadingPhotonID);
	tree->Branch("subleadingPhotonISO", &subleadingPhotonISO);
	tree->Branch("Jets", &jets);
	tree->Branch("Jets_bDiscriminant", &jets_bDiscriminant);
	tree->Branch("Jets_PUid", &jets_PUid);


}

// ------------ method called once each job just after ending the event loop  ------------
void 
bbggTree::endJob() 
{
	outFile->cd();
//	tree->Write();
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
