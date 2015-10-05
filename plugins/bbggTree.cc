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

class bbggTree : public edm::EDAnalyzer {
public:
    explicit bbggTree(const edm::ParameterSet&);
    ~bbggTree();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    typedef math::XYZTLorentzVector LorentzVector;
    typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    bbggTools tools_;
    //Parameter tokens
    edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diPhotonToken_;
    std::vector<edm::InputTag> inputTagJets_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > genToken_;

    edm::InputTag rhoFixedGrid_;
    std::string bTagType;
    unsigned int doSelection; 
	  
    //Tree objects
    LorentzVector leadingPhoton;
    LorentzVector subleadingPhoton;
    LorentzVector diphotonCandidate;
    vector<int> leadingPhotonID; 
    vector<int> leadingPhotonISO;
    vector<int> subleadingPhotonID;
    vector<int> subleadingPhotonISO;
    vector<LorentzVector> jets;
    vector<float> jets_bDiscriminant;
    vector<int> jets_PUid;
    LorentzVector leadingJet;
    float leadingJet_bDis;
    LorentzVector subleadingJet;
    float subleadingJet_bDis;
    LorentzVector dijetCandidate;
    LorentzVector diHiggsCandidate;
    vector<double> genWeights;
    double genTotalWeight;
    unsigned int nPromptInDiPhoton;
    int leadingPhotonEVeto;
    int subleadingPhotonEVeto;

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
inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
genToken_( consumes<edm::View<pat::PackedGenParticle> >( iConfig.getUntrackedParameter<edm::InputTag>( "GenTag", edm::InputTag( "prunedGenParticles" ) ) ) )
{
    //now do what ever initialization is needed
    tools_ = bbggTools();
    EvtCount = 0;
    //Default values for thresholds
    std::string def_bTagType;
    unsigned int def_doSelection = 0;
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

    rhoFixedGrid_  = iConfig.getUntrackedParameter<edm::InputTag>( "rhoFixedGridCollection", edm::InputTag( "fixedGridRhoAll" ) );

    bTagType = iConfig.getUntrackedParameter<std::string>( "bTagType", def_bTagType );

    fileName = iConfig.getUntrackedParameter<std::string>( "OutFileName", def_fileName );
	
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
    
	  
    std::cout << "Parameters initialized... \n ############ Doing selection tree or before selection tree? : " << (doSelection ? "Selection!":"Before selection!") <<  std::endl;

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
    if( EvtCount%100 == 0 && !DEBUG) std::cout << "[bbggTree::analyze] Analyzing event number: " << EvtCount << std::endl;
    if (DEBUG) std::cout << "[bbggTree::analyze] Analyzing event number: " << EvtCount << std::endl;
	
    EvtCount++;
    using namespace edm;
    leadingPhotonID.clear();
    leadingPhotonISO.clear();
    subleadingPhotonID.clear();
    subleadingPhotonISO.clear();
    jets.clear();
    jets_bDiscriminant.clear();
    jets_PUid.clear();

    //Get Jets collections!
    JetCollectionVector theJetsCols( inputTagJets_.size() );
    for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
        iEvent.getByLabel( inputTagJets_[j], theJetsCols[j] );
    }

    if (DEBUG) std::cout << "Number of jet collections!!!! " << theJetsCols.size() << std::endl;
	 
    Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
    iEvent.getByToken( diPhotonToken_, diPhotons );

    Handle<double> rhoHandle;        // the old way for now...move to getbytoken?
    iEvent.getByLabel( rhoFixedGrid_, rhoHandle );
    const double rhoFixedGrd = *( rhoHandle.product() );
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
                else if (nPrompt == 1) diphoVec.push_back(dipho);
            }
        } else {
            diphoVec.push_back(dipho);
        }
    }
    
    if(DEBUG) std::cout << "Number of diphoton candidates: " << diphoVec.size() << std::endl;
	
    if (diphoVec.size() < 1) return;
    if (theJetsCols.size() < 1) return;
	
    if(doSelection) {
        bool cutsChecked = tools_.CheckCuts();
        if(!cutsChecked) {
            std::cout << "You haven't filled all the cuts correctly!" << std::endl;
            return;
        }

        if(DEBUG) std::cout << "[bbggTree::analyze] About to do event selection! " << std::endl;

        bool passedSelection = tools_.AnalysisSelection(diphoVec, theJetsCols);
        if(!passedSelection) return;
		
        if(DEBUG) std::cout << "[bbggTree::analyze] tools_.AnalysisSelection returned " << passedSelection << std::endl;
		
        edm::Ptr<flashgg::DiPhotonCandidate> diphoCand = tools_.GetSelected_diphoCandidate();
        edm::Ptr<flashgg::Jet> LeadingJet = tools_.GetSelected_leadingJetCandidate();
        edm::Ptr<flashgg::Jet> SubLeadingJet = tools_.GetSelected_subleadingJetCandidate();
        
        if (DEBUG) std::cout << "Jet collection picked: " << diphoCand->jetCollectionIndex() << std::endl;
		
        nPromptInDiPhoton = 999;
        if ( genInfo.isValid() ){
            bbggMC _mcTools = bbggMC();
            nPromptInDiPhoton = _mcTools.CheckNumberOfPromptPhotons(diphoCand, genParticles);
        }
		
        diphotonCandidate = diphoCand->p4();
        leadingPhoton = diphoCand->leadingPhoton()->p4();
        subleadingPhoton = diphoCand->subLeadingPhoton()->p4();
        leadingJet = LeadingJet->p4();
        leadingJet_bDis = LeadingJet->bDiscriminator(bTagType);
        subleadingJet = SubLeadingJet->p4();
        subleadingJet_bDis = SubLeadingJet->bDiscriminator(bTagType);
        dijetCandidate = leadingJet + subleadingJet;
        diHiggsCandidate = diphotonCandidate + dijetCandidate;
		
        if(DEBUG) std::cout << "[bbggTree::analyze] After filling candidates" << std::endl;
				
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
		
        leadingPhotonEVeto = diphoCand->leadingPhoton()->passElectronVeto();
        subleadingPhotonEVeto = diphoCand->subLeadingPhoton()->passElectronVeto();
				
        if(DEBUG) std::cout << "GOT TO THE END!!" << std::endl;
        tree->Fill();

        if(DEBUG) std::cout << "Histograms filled!" << std::endl;
		
    } else {
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
        if(dipho_pt_ref < 1) return;
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
		
        leadingPhotonEVeto = diphoCand->leadingPhoton()->passElectronVeto();
        subleadingPhotonEVeto = diphoCand->subLeadingPhoton()->passElectronVeto();

        unsigned int jetCollectionIndex = diphoCand->jetCollectionIndex();
        for( unsigned int jetIndex = 0; jetIndex < theJetsCols[jetCollectionIndex]->size(); jetIndex++ )
        {
            edm::Ptr<flashgg::Jet> jet = theJetsCols[jetCollectionIndex]->ptrAt( jetIndex );
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
}


// ------------ method called once each job just before starting event loop  ------------
void 
    bbggTree::beginJob()
{
    if(DEBUG) std::cout << "[bbggTree::beginJob] Setting up output tree..." << std::endl;

    outFile = new TFile(fileName.c_str(), "RECREATE");
    if(!doSelection) {
        tree = new TTree("bbggTree", "Flat tree for HH->bbgg analyses (no selection)");
        tree->Branch("genWeights", &genWeights);
        tree->Branch("genTotalWeight", &genTotalWeight, "genTotalWeight/D");
        tree->Branch("leadingPhoton", &leadingPhoton);
        tree->Branch("leadingPhotonID", &leadingPhotonID);
        tree->Branch("leadingPhotonISO", &leadingPhotonISO);
        tree->Branch("leadingPhotonEVeto", &leadingPhotonEVeto, "leadingPhotonEVeto/I");
        tree->Branch("subleadingPhoton", &subleadingPhoton);
        tree->Branch("subleadingPhotonID", &subleadingPhotonID);
        tree->Branch("subleadingPhotonISO", &subleadingPhotonISO);
        tree->Branch("subleadingPhotonEVeto", &subleadingPhotonEVeto, "subleadingPhotonEVeto/I");
        tree->Branch("nPromptInDiPhoton", &nPromptInDiPhoton, "nPromptInDiPhoton/I");
        tree->Branch("Jets", &jets);
        tree->Branch("Jets_bDiscriminant", &jets_bDiscriminant);
        tree->Branch("Jets_PUid", &jets_PUid);
    }
    if(doSelection) {
        tree = new TTree("bbggSelectionTree", "Flat tree for HH->bbgg analyses (after pre selection)");
        tree->Branch("genWeights", &genWeights);
        tree->Branch("genTotalWeight", &genTotalWeight, "genTotalWeight/D");
        tree->Branch("leadingPhoton", &leadingPhoton);
        tree->Branch("leadingPhotonID", &leadingPhotonID);
        tree->Branch("leadingPhotonISO", &leadingPhotonISO);
        tree->Branch("leadingPhotonEVeto", &leadingPhotonEVeto, "leadingPhotonEVeto/I");
        tree->Branch("subleadingPhoton", &subleadingPhoton);
        tree->Branch("subleadingPhotonID", &subleadingPhotonID);
        tree->Branch("subleadingPhotonISO", &subleadingPhotonISO);
        tree->Branch("subleadingPhotonEVeto", &subleadingPhotonEVeto, "subleadingPhotonEVeto/I");
        tree->Branch("diphotonCandidate", &diphotonCandidate);
        tree->Branch("nPromptInDiPhoton", &nPromptInDiPhoton, "nPromptInDiPhoton/I");
        tree->Branch("leadingJet", &leadingJet);
        tree->Branch("leadingJet_bDis", &leadingJet_bDis, "leadingJet_bDis/F");
        tree->Branch("subleadingJet", &subleadingJet);
        tree->Branch("subleadingJet_bDis", &subleadingJet_bDis, "subleadingJet_bDis/F");
        tree->Branch("dijetCandidate", &dijetCandidate);
        tree->Branch("diHiggsCandidate", &diHiggsCandidate);
		
    }

    if(DEBUG) std::cout << "[bbggTree::beginJob] Output tree set!" << std::endl;

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
