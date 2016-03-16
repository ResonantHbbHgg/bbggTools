// -*- C++ -*-
//
// Package:    flashgg/bbggBTagEff
// Class:      bbggBTagEff
// 
/**\class bbggBTagEff bbggBTagEff.cc flashgg/bbggBTagEff/plugins/bbggBTagEff.cc

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
#include <map>
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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//FLASHgg files
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/Taggers/interface/GlobalVariablesDumper.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

//Local
#include "flashgg/bbggTools/interface/bbggTools.h"
#include "flashgg/bbggTools/interface/bbggMC.h"

//
// class declaration
//

const int DEBUG = 0;

class bbggBTagEff : public edm::EDAnalyzer {
public:
    explicit bbggBTagEff(const edm::ParameterSet&);
    ~bbggBTagEff();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    typedef math::XYZTLorentzVector LorentzVector;
    typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    bbggTools tools_;
    flashgg::GlobalVariablesDumper* globVar_;
    //Parameter tokens
    edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diPhotonToken_;
    edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > BeforeSelectionDiPhotonToken_;
    std::vector<edm::InputTag> inputTagJets_;
    std::vector<edm::EDGetTokenT<edm::View<flashgg::Jet> > > tokenJets_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > genToken_;
    edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
    

//    edm::InputTag rhoFixedGrid_;
	  
    //Tree objects
    LorentzVector leadingPhoton, subleadingPhoton, diphotonCandidate;
    LorentzVector leadingJet, subleadingJet, dijetCandidate;
    LorentzVector diHiggsCandidate;
    vector<int> leadingPhotonID, leadingPhotonISO, subleadingPhotonID, subleadingPhotonISO;
    vector<double> genWeights;
    float leadingJet_bDis, subleadingJet_bDis;
    double genTotalWeight;
    unsigned int nPromptInDiPhoton;
    int leadingPhotonEVeto, subleadingPhotonEVeto;
    int leadingJet_flavour, subleadingJet_flavour;
    int isSignal, isPhotonCR;
//    int nvtx;
    int passed0, passed1, passed2, passed3, passed4, passed5, passed6, passed7, passed8, passed9, passed10, passed11, passed12, passed13, passed14;

    vector<LorentzVector> leadingjets, subleadingjets, dijets;
    vector<double> leadingjets_bDiscriminant, subleadingjets_bDiscriminant;
    vector<int> leadingjets_partonID, subleadingjets_partonID;
    vector<double> leadingJets_DRleadingPho, leadingJets_DRsubleadingPho;
    vector<double> subleadingJets_DRleadingPho, subleadingJets_DRsubleadingPho;
    vector<double> Jets_minDRPho, Jets_DR, DiJets_DRDiPho;

    //Parameters
    std::vector<double> phoIDlooseEB, phoIDlooseEE;
    std::vector<double> phoIDmediumEB, phoIDmediumEE;
    std::vector<double> phoIDtightEB, phoIDtightEE;
    std::vector<double> phoISOlooseEB, phoISOlooseEE;
    std::vector<double> phoISOmediumEB, phoISOmediumEE;
    std::vector<double> phoISOtightEB,phoISOtightEE;
    std::vector<double> nhCorrEB, nhCorrEE;
    std::vector<double> phCorrEB, phCorrEE;
    std::vector<double> ph_pt, ph_eta, ph_r9;
    std::vector<double> diph_pt, diph_eta, diph_mass, MVAPhotonID;
    std::vector<int> ph_elVeto, ph_doelVeto, ph_doID, ph_doISO;
    std::vector<std::string> ph_whichID, ph_whichISO;
    unsigned int diph_onlyfirst;
    std::vector<double> jt_pt, jt_eta, jt_drPho, jt_bDis;
    std::vector<int> jt_doPU, jt_doID;
    unsigned int n_bJets;
    std::vector<double> dijt_pt, dijt_eta, dijt_mass;
    std::vector<double> cand_pt, cand_eta, cand_mass, dr_cands;
    unsigned int nPromptPhotons, doDoubleCountingMitigation, doPhotonCR, doJetRegression;
    std::vector<std::string> myTriggers;
    std::string bTagType, PhotonMVAEstimator;
    unsigned int doSelection, DoMVAPhotonID; 
    

    //OutFile & Hists
    TFile* outFile;
//    TTree* tree;
    std::string fileName;
	  
    //Event counter for cout's
    long unsigned int EvtCount;

   //Histogram
   TH2F* pass2d;
   TH2F* fail2d;
   TH2F* all2d;
   TH2F* passEff;
   TH2F* failEff;
   TH1F* Efficiencies;
   TH1F* Mjj;
   TH1F* Mgg;
   TH1F* Mggjj;
};

//
// constructors and destructor
//
bbggBTagEff::bbggBTagEff(const edm::ParameterSet& iConfig) :
diPhotonToken_( consumes<edm::View<flashgg::DiPhotonCandidate> >( iConfig.getUntrackedParameter<edm::InputTag> ( "DiPhotonTag", edm::InputTag( "flashggDiPhotons" ) ) ) ),
inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
genToken_( consumes<edm::View<pat::PackedGenParticle> >( iConfig.getUntrackedParameter<edm::InputTag>( "GenTag", edm::InputTag( "prunedGenParticles" ) ) ) )
{

   BeforeSelectionDiPhotonToken_ = consumes<edm::View<flashgg::DiPhotonCandidate> >( iConfig.getUntrackedParameter<edm::InputTag> ( "SimpleDiPhotonTag", edm::InputTag( "flashggDiPhotons" ) ) ); 



    //now do what ever initialization is needed
    tools_ = bbggTools();
    globVar_ = new flashgg::GlobalVariablesDumper(iConfig);
    //Lumi weight
    double lumiWeight_ = ( iConfig.getParameter<double>( "lumiWeight" ) );
    globVar_->dumpLumiFactor(lumiWeight_);
    EvtCount = 0;
    //Default values for thresholds
    std::string def_bTagType, def_PhotonMVAEstimator;
    unsigned int def_doSelection = 0;
    std::vector<double> def_ph_pt, def_ph_eta, def_ph_r9;
    std::vector<double> def_diph_pt, def_diph_eta, def_diph_mass, def_MVAPhotonID;
    std::vector<double> def_jt_pt, def_jt_eta, def_jt_drPho, def_jt_bDis;
    std::vector<double> def_dijt_pt, def_dijt_eta, def_dijt_mass, def_cand_pt, def_cand_eta, def_cand_mass, def_dr_cands;
    std::vector<int> def_ph_elVeto, def_ph_doelVeto, def_ph_doID;
    std::vector<int> def_ph_doISO;
    std::vector<int> def_jt_doPU, def_jt_doID;
    std::vector<std::string> def_ph_whichID, def_ph_whichISO;
    unsigned int def_diph_onlyfirst;
    unsigned int def_n_bJets;

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
    unsigned int def_doPhotonCR = 0;
    unsigned int def_DoMVAPhotonID = 0;

    std::vector<std::string> def_myTriggers;


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
    
    doSelection = iConfig.getUntrackedParameter<unsigned int>("doSelectionTree", def_doSelection);
    ph_pt     = iConfig.getUntrackedParameter<std::vector<double > >("PhotonPtOverDiPhotonMass", def_ph_pt);
    ph_eta    = iConfig.getUntrackedParameter<std::vector<double > >("PhotonEta", def_ph_eta);
    ph_r9     = iConfig.getUntrackedParameter<std::vector<double > >("PhotonR9", def_ph_r9);
    ph_elVeto = iConfig.getUntrackedParameter<std::vector<int > >("PhotonElectronVeto", def_ph_elVeto);
    ph_doelVeto = iConfig.getUntrackedParameter<std::vector<int > >("PhotonDoElectronVeto", def_ph_doelVeto);
    ph_doID   = iConfig.getUntrackedParameter<std::vector<int > >("PhotonDoID", def_ph_doID);
    ph_whichID   = iConfig.getUntrackedParameter<std::vector<std::string > >("PhotonWhichID", def_ph_whichID);
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
    doPhotonCR = iConfig.getUntrackedParameter<unsigned int>("doPhotonCR", def_doPhotonCR);
//    rhoFixedGrid_  = iConfig.getUntrackedParameter<edm::InputTag>( "rhoFixedGridCollection", edm::InputTag( "fixedGridRhoAll" ) );
    bTagType = iConfig.getUntrackedParameter<std::string>( "bTagType", def_bTagType );
    fileName = iConfig.getUntrackedParameter<std::string>( "OutFileName", def_fileName );
    
    myTriggers = iConfig.getUntrackedParameter<std::vector<std::string> >("myTriggers", def_myTriggers);
    triggerToken_ = consumes<edm::TriggerResults>( iConfig.getParameter<edm::InputTag>( "triggerTag" ) );
    
    //MVA Photon ID
    DoMVAPhotonID = iConfig.getUntrackedParameter<unsigned int>("DoMVAPhotonID", def_DoMVAPhotonID);
    MVAPhotonID = iConfig.getUntrackedParameter<std::vector<double>>("MVAPhotonID", def_MVAPhotonID);
    PhotonMVAEstimator = iConfig.getUntrackedParameter<std::string>("PhotonMVAEstimator", def_PhotonMVAEstimator);
    
    tools_.SetCut_DoMVAPhotonID(DoMVAPhotonID);
    tools_.SetCut_MVAPhotonID(MVAPhotonID);
    tools_.SetCut_PhotonMVAEstimator(PhotonMVAEstimator);
	
    tools_.SetPhotonCR(doPhotonCR);
    tools_.SetCut_PhotonPtOverDiPhotonMass( ph_pt );
    tools_.SetCut_PhotonEta( ph_eta );
    tools_.SetCut_PhotonDoID( ph_doID );
    tools_.SetCut_PhotonDoISO( ph_doISO );
    tools_.SetCut_PhotonR9( ph_r9 );
    tools_.SetCut_PhotonElectronVeto( ph_elVeto );
    tools_.SetCut_PhotonDoElectronVeto( ph_doelVeto );
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
        
    /*
    tools_.SetCut_phoIDlooseEB(phoIDlooseEB);
    tools_.SetCut_phoIDlooseEE(phoIDlooseEE);
    tools_.SetCut_phoIDmediumEB(phoIDmediumEB);
    tools_.SetCut_phoIDmediumEE(phoIDmediumEE);
    tools_.SetCut_phoIDtightEB(phoIDtightEB);
    tools_.SetCut_phoIDtightEE(phoIDtightEE);
    */
    
    /*        
    tools_.SetCut_phoISOlooseEB(phoISOlooseEB);
    tools_.SetCut_phoISOlooseEE(phoISOlooseEE);
    tools_.SetCut_phoISOmediumEB(phoISOmediumEB);
    tools_.SetCut_phoISOmediumEE(phoISOmediumEE);
    tools_.SetCut_phoISOtightEB(phoISOtightEB);
    tools_.SetCut_phoISOtightEE(phoISOtightEE);
    */
    
    std::map<int, vector<double> > nhCorr;
    nhCorr[0] = nhCorrEB;
    nhCorr[1] = nhCorrEE;
    tools_.SetCut_nhCorr(nhCorr);
    
    std::map<int, vector<double> > phCorr;
    phCorr[0] = phCorrEB;
    phCorr[1] = phCorrEE;
    tools_.SetCut_phCorr(phCorr);
    
    /*
    tools_.SetCut_nhCorrEB(nhCorrEB);
    tools_.SetCut_nhCorrEE(nhCorrEE);
    tools_.SetCut_phCorrEB(phCorrEB);
    tools_.SetCut_phCorrEE(phCorrEE);
    */
    
    tools_.SetCut_phoWhichID(ph_whichID);
    tools_.SetCut_phoWhichISO(ph_whichISO);

    for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
	auto token = consumes<edm::View<flashgg::Jet> >(inputTagJets_[i]);
	tokenJets_.push_back(token);
    }
    
	  
    std::cout << "Parameters initialized... \n ############ Doing selection tree or before selection tree? : " << (doSelection ? "Selection!":"Before selection!") <<  std::endl;

}


bbggBTagEff::~bbggBTagEff()
{
 
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
    bbggBTagEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    using namespace edm;

    //Get Jets collections!
    int totalNumberofJets = 0;
    JetCollectionVector theJetsCols( inputTagJets_.size() );
    for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
        iEvent.getByToken( tokenJets_[j], theJetsCols[j] );
	totalNumberofJets += theJetsCols[j]->size();
    }

    if (DEBUG) std::cout << "Number of jet collections!!!! " << theJetsCols.size() << std::endl;
	 
    Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
    iEvent.getByToken( diPhotonToken_, diPhotons );

    Handle<View<flashgg::DiPhotonCandidate> > BeforeSeldiPhotons;
    iEvent.getByToken( BeforeSelectionDiPhotonToken_, BeforeSeldiPhotons );

    int totalNumberofPhotons = BeforeSeldiPhotons->size();

/*
    Handle<double> rhoHandle;        // the old way for now...move to getbytoken?
    iEvent.getByLabel( rhoFixedGrid_, rhoHandle );
    const double rhoFixedGrd = *( rhoHandle.product() );
    tools_.setRho(rhoFixedGrd);
*/
    const double rhoFixedGrd = globVar_->valueOf(globVar_->indexOf("rho"));
    tools_.setRho(rhoFixedGrd);

    Handle<View<pat::PackedGenParticle> > genParticles;

    //MC Weights
    Handle<GenEventInfoProduct> genInfo;
    try {
        iEvent.getByLabel( "generator", genInfo );
    } catch (...) {;}
    if( genInfo.isValid() ) {
        genTotalWeight = genInfo->weight();
        iEvent.getByToken( genToken_, genParticles);
    } else {
        genTotalWeight = 1;
    }

    //PreLoop
    vector<edm::Ptr<flashgg::DiPhotonCandidate>> diphoVec;
    for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ )
    {
        edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( diphoIndex );
        if(doDoubleCountingMitigation){
            if ( !genInfo.isValid() ) {
                std::cout << "[bbggTree::analyze] Oh, man! You're trying to get MC information from data. Pay attention to what you're doing!" << std::endl;
                return;
            }
            bbggMC _mcTools = bbggMC();
            unsigned int nPrompt = _mcTools.CheckNumberOfPromptPhotons(dipho, genParticles);
            if (nPromptPhotons > 0 && nPrompt == nPromptPhotons) diphoVec.push_back(dipho);
            if (nPromptPhotons == 0) {
                if (nPrompt == 0) diphoVec.push_back(dipho);
//                else if (nPrompt == 1) diphoVec.push_back(dipho);
            }
        } else {
            diphoVec.push_back(dipho);
        }
    }


passed0 = 0;
passed1 = 0;
passed2 = 0;
passed3 = 0;
passed4 = 0;
passed5 = 0;
passed6 = 0;
passed7 = 0;
passed8 = 0;
passed9 = 0;
passed10 = 0;
passed11 = 0;
passed12 = 0;
passed13 = 0;
passed14 = 0;


//Start efficiencies filling
Efficiencies->Fill(0);
passed0 = 1;

//Check number of diphotons and jets before any sel
if( totalNumberofJets < 2 || totalNumberofPhotons < 1) { return;}
Efficiencies->Fill(1);
passed1 = 1;

//Check number of diphotons after presel:
//76x: Check trigger results
//Trigger
std::vector<int> myTriggerResults;
bool triggerAccepted = false;
if(myTriggers.size() > 0){
    Handle<edm::TriggerResults> trigResults;
    iEvent.getByToken(triggerToken_, trigResults);
    const edm::TriggerNames &names = iEvent.triggerNames(*trigResults);
    myTriggerResults = tools_.TriggerSelection(myTriggers, names, trigResults);
    for( unsigned int tr = 0; tr < myTriggerResults.size(); tr++){
        if(myTriggerResults[tr] == 1){
            triggerAccepted = true;
            break;
        }
    }
}
if( !triggerAccepted ) {return;}
Efficiencies->Fill(2);
passed2 = 1;

//Check photon kinematic selection:
std::vector<edm::Ptr<flashgg::DiPhotonCandidate>> PreSelDiPhos = tools_.DiPhotonKinematicSelection(diphoVec);
if(PreSelDiPhos.size() < 1 ) {return;}
Efficiencies->Fill(3);
passed3 = 1;

//Check photon ID:
std::vector<edm::Ptr<flashgg::DiPhotonCandidate>> SelDiPhos = tools_.DiPhotonIDSelection(PreSelDiPhos);
if(SelDiPhos.size() < 1) {return;}
Efficiencies->Fill(4);
passed4 = 1;

//Check jet presel:
edm::Ptr<flashgg::DiPhotonCandidate> diphoCand = SelDiPhos[0];
std::vector<edm::Ptr<flashgg::Jet>> PreSelJets = tools_.JetPreSelection(theJetsCols, diphoCand);
//Loop over jets
for(unsigned int j = 0; j < PreSelJets.size(); j++)
{
  if(j>4) break;
  edm::Ptr<flashgg::Jet> thisJet = PreSelJets[j];
  float thisEta = fabs(thisJet->eta());
  float thisPt = thisJet->pt();
  float bdis = thisJet->bDiscriminator(bTagType);
  all2d->Fill(thisEta, thisPt, 1.);
  if(bdis > 0.8)
     pass2d->Fill(thisEta, thisPt, 1.);
  else
     fail2d->Fill(thisEta, thisPt, 1.);
}

return;
}


// ------------ method called once each job just before starting event loop  ------------
void 
    bbggBTagEff::beginJob()
{
    if(DEBUG) std::cout << "[bbggBTagEff::beginJob] Setting up output tree..." << std::endl;

    outFile = new TFile(fileName.c_str(), "RECREATE");
    Float_t ptbins[] = {25, 40, 60, 100, 150, 200, 300, 400, 500};
    int nptbins = 8;
    Float_t etabins[] = {0, 0.5, 1, 1.5, 2, 2.5};
    int netabins = 5;

    all2d = new TH2F("all2d", "All jets in PtxEta;Eta;Pt", netabins, etabins, nptbins, ptbins);
    pass2d = new TH2F("pass2d", "Btagged jets in PtxEta;Eta;Pt", netabins, etabins, nptbins, ptbins);
    fail2d = new TH2F("fail2d", "Bvetoed jets in PtxEta;Eta;Pt", netabins, etabins, nptbins, ptbins);

    Efficiencies = new TH1F("Efficiencies", "Efficiencies", 10, 0, 10);
    Mjj = new TH1F("Mjj", "Mjj", 100, 60, 200);
    Mgg = new TH1F("Mgg", "Mgg", 100, 100, 150);
    Mggjj = new TH1F("Mggjj", "Mggjj", 100, 0, 1000);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
    bbggBTagEff::endJob() 
{
    passEff = (TH2F*) pass2d->Clone("passEff");
    failEff = (TH2F*) fail2d->Clone("failEff");
    passEff->Sumw2();
    failEff->Sumw2();
    TH2F* all = (TH2F*) all2d->Clone("all");
    all->Sumw2();
    passEff->Divide(all);
    failEff->Divide(all);

    outFile->cd();
//    tree->Write();
    Efficiencies->Write();
    Mgg->Write();
    Mjj->Write();
    Mggjj->Write();
    all2d->Write();
    pass2d->Write();
    fail2d->Write();
    passEff->Write();
    failEff->Write();
    outFile->Write();
    outFile->Close();
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
bbggBTagEff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(bbggBTagEff);
