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
#include "TLorentzVector.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"


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
#include "flashgg/bbggTools/interface/bbggKinFit.h"
#include "flashgg/bbggTools/interface/bbggJetRegression.h"

//
// class declaration
//

const int DEBUG = 0;

class bbggTree : public edm::EDAnalyzer {
public:
//    explicit bbggTree(const edm::ParameterSet&, edm::ConsumesCollector && cc);
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
    bbggKinFit kinFit_;
    bbggJetRegression jetReg_;
    flashgg::GlobalVariablesDumper* globVar_;
    //Parameter tokens
    edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diPhotonToken_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > genToken_;
//    edm::EDGetTokenT<double> rhoFixedGrid_;
    std::vector<edm::InputTag> inputTagJets_;
    std::vector<edm::EDGetTokenT<edm::View<flashgg::Jet> > > tokenJets_;
    
    std::string bTagType;
    unsigned int doSelection; 
	  
    //Tree objects
    LorentzVector leadingPhoton, subleadingPhoton, diphotonCandidate;
    LorentzVector leadingJet, subleadingJet, dijetCandidate;
    LorentzVector leadingJet_KF, subleadingJet_KF, dijetCandidate_KF;
    LorentzVector leadingJet_Reg, subleadingJet_Reg, dijetCandidate_Reg;
    LorentzVector leadingJet_RegKF, subleadingJet_RegKF, dijetCandidate_RegKF;
    LorentzVector diHiggsCandidate, diHiggsCandidate_KF, diHiggsCandidate_Reg,diHiggsCandidate_RegKF;
    vector<int> leadingPhotonID, leadingPhotonISO, subleadingPhotonID, subleadingPhotonISO;
    vector<double> genWeights;
    float leadingJet_bDis, subleadingJet_bDis, jet1PtRes, jet1EtaRes, jet1PhiRes, jet2PtRes, jet2EtaRes, jet2PhiRes;
    float CosThetaStar;
    
    double genTotalWeight;
    unsigned int nPromptInDiPhoton;
    int leadingPhotonEVeto, subleadingPhotonEVeto;
    int leadingJet_flavour, subleadingJet_flavour;
    int isSignal, isPhotonCR;
//    int nvtx;

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
    std::vector<double> diph_pt, diph_eta, diph_mass;
    std::vector<int> ph_elVeto, ph_doelVeto, ph_doID, ph_doISO;
    std::vector<std::string> ph_whichID, ph_whichISO;
    unsigned int diph_onlyfirst;
    std::vector<double> jt_pt, jt_eta, jt_drPho, jt_bDis;
    std::vector<int> jt_doPU, jt_doID;
    unsigned int n_bJets;
    std::vector<double> dijt_pt, dijt_eta, dijt_mass;
    std::vector<double> cand_pt, cand_eta, cand_mass, dr_cands;
    unsigned int nPromptPhotons, doDoubleCountingMitigation, doPhotonCR;

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
genToken_( consumes<edm::View<pat::PackedGenParticle> >( iConfig.getUntrackedParameter<edm::InputTag>( "GenTag", edm::InputTag( "prunedGenParticles" ) ) ) ),
//rhoFixedGrid_(consumes<double>(iConfig.getParameter<edm::InputTag>( "rhoFixedGridCollection" ) ) ),
inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) )
{
    //now do what ever initialization is needed
    tools_ = bbggTools();
    jetReg_ = bbggJetRegression();
    jetReg_.SetupRegression("BDTG method", "/afs/cern.ch/work/r/rateixei/work/DiHiggs/flashggJets/CMSSW_7_4_15/src/flashgg/bbggTools/Weights/BRegression/TMVARegression_BDTG.weights.xml");
//    globVar_ = new flashgg::GlobalVariablesDumper(iConfig,std::forward<edm::ConsumesCollector>(cc));
    globVar_ = new flashgg::GlobalVariablesDumper(iConfig, consumesCollector() );
    //Lumi weight
    double lumiWeight_ = ( iConfig.getParameter<double>( "lumiWeight" ) );
    globVar_->dumpLumiFactor(lumiWeight_);
    EvtCount = 0;
    //Default values for thresholds
    std::string def_bTagType;
    unsigned int def_doSelection = 0;
    std::vector<double> def_ph_pt, def_ph_eta, def_ph_r9;
    std::vector<double> def_diph_pt, def_diph_eta, def_diph_mass;
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
    
    std::vector<double> ptRes = iConfig.getUntrackedParameter<std::vector<double>>( "ptRes");
    std::vector<double> etaRes = iConfig.getUntrackedParameter<std::vector<double>>( "etaRes");
    std::vector<double> phiRes = iConfig.getUntrackedParameter<std::vector<double>>( "phiRes");   
    std::vector<double> etaBins = iConfig.getUntrackedParameter<std::vector<double>>("etaBins");
 
    kinFit_ = bbggKinFit();
    kinFit_.SetJetResolutionParameters(etaBins, ptRes, etaRes, phiRes);

    for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
	auto token = consumes<edm::View<flashgg::Jet> >(inputTagJets_[i]);
	tokenJets_.push_back(token);
    }

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

    globVar_->fill(iEvent);

    if( EvtCount%100 == 0 && !DEBUG) std::cout << "[bbggTree::analyze] Analyzing event number: " << EvtCount << std::endl;
    if (DEBUG) std::cout << "[bbggTree::analyze] Analyzing event number: " << EvtCount << std::endl;
	
    EvtCount++;
    using namespace edm;
    leadingPhotonID.clear();
    leadingPhotonISO.clear();
    subleadingPhotonID.clear();
    subleadingPhotonISO.clear();

    leadingjets.clear();
    subleadingjets.clear();
    dijets.clear();
    leadingjets_bDiscriminant.clear();
    subleadingjets_bDiscriminant.clear();
    leadingjets_partonID.clear();
    subleadingjets_partonID.clear();
    leadingJets_DRleadingPho.clear();
    leadingJets_DRsubleadingPho.clear();
    leadingJets_DRleadingPho.clear();
    subleadingJets_DRsubleadingPho.clear();
    Jets_minDRPho.clear();
    Jets_DR.clear();
    DiJets_DRDiPho.clear();
    isSignal = 0;
    isPhotonCR = 0;
    CosThetaStar = -999;
//    nvtx = 0;

    diphotonCandidate.SetPxPyPzE(0,0,0,0);// = diphoCand->p4();
    leadingPhoton.SetPxPyPzE(0,0,0,0);// = diphoCand->leadingPhoton()->p4();
    subleadingPhoton.SetPxPyPzE(0,0,0,0);// = diphoCand->subLeadingPhoton()->p4();
    leadingJet.SetPxPyPzE(0,0,0,0);// = LeadingJet->p4();
    leadingJet_bDis = 0;// = LeadingJet->bDiscriminator(bTagType);
    leadingJet_flavour = 0;// = LeadingJet->partonFlavour();
    subleadingJet.SetPxPyPzE(0,0,0,0);// = SubLeadingJet->p4();
    subleadingJet_bDis = 0;//SubLeadingJet->bDiscriminator(bTagType);
    subleadingJet_flavour = 0;//SubLeadingJet->partonFlavour();
    dijetCandidate.SetPxPyPzE(0,0,0,0);// = leadingJet + subleadingJet;
    diHiggsCandidate.SetPxPyPzE(0,0,0,0);// = diphotonCandidate + dijetCandidate;

    //Get Jets collections!
    JetCollectionVector theJetsCols( inputTagJets_.size() );
    for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
        iEvent.getByToken( tokenJets_[j], theJetsCols[j] );
    }

    if (DEBUG) std::cout << "Number of jet collections!!!! " << theJetsCols.size() << std::endl;
	 
    Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
    iEvent.getByToken( diPhotonToken_, diPhotons );

    // Handle<double> rhoHandle;        // the old way for now...move to getbytoken?
    // iEvent.getByToken( rhoFixedGrid_, rhoHandle );
    //
    // const double rhoFixedGrd = *( rhoHandle.product() );
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
	isSignal = tools_.IsSignal();
	isPhotonCR = tools_.IsPhotonCR();
	bool hasLJet = tools_.HasLeadJet();
	bool hasSJet = tools_.HasSubJet();

    if(!isSignal && !isPhotonCR) return; //if event is not signal and is not photon control region, skip
	if(!isSignal && !doPhotonCR) return; //if event is not signal and you don't want to save photon control region, skip
	//if(!passedSelection) return;
	if(!hasLJet || !hasSJet) return;
		
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
	leadingJet_flavour = LeadingJet->partonFlavour();
        subleadingJet = SubLeadingJet->p4();
        subleadingJet_bDis = SubLeadingJet->bDiscriminator(bTagType);
	subleadingJet_flavour = SubLeadingJet->partonFlavour();
        dijetCandidate = leadingJet + subleadingJet;
        diHiggsCandidate = diphotonCandidate + dijetCandidate;
        
        //Kinematic fit:
        if(DEBUG) std::cout << "[bbggTree::analyze] Doing kinematic fit!" << std::endl;
        kinFit_.KinematicFit(leadingJet, subleadingJet);
        leadingJet_KF = kinFit_.GetJet1();
        subleadingJet_KF = kinFit_.GetJet2();
        dijetCandidate_KF = leadingJet_KF + subleadingJet_KF;
        diHiggsCandidate_KF = dijetCandidate_KF + diphotonCandidate;
        jet1PtRes = kinFit_.GetPtResolution(leadingJet);
        jet1EtaRes = kinFit_.GetEtaResolution(leadingJet);
        jet1PhiRes= kinFit_.GetPhiResolution(leadingJet);

	//Regression:
        if(DEBUG) std::cout << "[bbggTree::analyze] Doing regression!" << std::endl;
	leadingJet_Reg = jetReg_.GetRegression(LeadingJet, rhoFixedGrd);
	subleadingJet_Reg = jetReg_.GetRegression(SubLeadingJet, rhoFixedGrd);
	dijetCandidate_Reg = leadingJet_Reg + subleadingJet_Reg;
	diHiggsCandidate_Reg = dijetCandidate_Reg + diphotonCandidate;

	//Regression + Kinematic Fit:
        kinFit_.KinematicFit(leadingJet_Reg, subleadingJet_Reg);
        leadingJet_RegKF = kinFit_.GetJet1();
        subleadingJet_RegKF = kinFit_.GetJet2();
        dijetCandidate_RegKF = leadingJet_RegKF + subleadingJet_RegKF;
        diHiggsCandidate_RegKF = dijetCandidate_RegKF + diphotonCandidate;

        //Calculating costheta star
        TLorentzVector BoostedHgg(0,0,0,0);
        BoostedHgg.SetPtEtaPhiE( diphoCand->pt(), diphoCand->eta(), diphoCand->phi(), diphoCand->energy());
        TLorentzVector HHforBoost(0,0,0,0);
        HHforBoost.SetPtEtaPhiE( diHiggsCandidate.pt(), diHiggsCandidate.eta(), diHiggsCandidate.phi(), diHiggsCandidate.energy());
        TVector3 HHBoostVector = HHforBoost.BoostVector();
        BoostedHgg.Boost( -HHBoostVector.x(), -HHBoostVector.y(), -HHBoostVector.z() );
        CosThetaStar = BoostedHgg.CosTheta();
        
        		
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

        std::vector<edm::Ptr<flashgg::DiPhotonCandidate>> passedSelection = tools_.DiPhotonKinematicSelection(diphoVec);
        if(passedSelection.size() == 0) return;
	edm::Ptr<flashgg::DiPhotonCandidate> diphoCand = tools_.GetSelected_diphoCandidate();

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
//	    if(tjet.pt() < 25 || fabs(tjet.eta()) > 2.4)
//		continue;
	    for( unsigned int jet2Index = jetIndex+1; jetIndex < theJetsCols[jetCollectionIndex]->size(); jetIndex++ )
	    {
            	edm::Ptr<flashgg::Jet> jet2 = theJetsCols[jetCollectionIndex]->ptrAt( jetIndex );
            	LorentzVector tjet2 = jet2->p4(); //LorentzVector(jet->px(), jet->py(), jet->pz(), jet->energy());
//	    	if(tjet2.pt() < 25 || fabs(tjet2.eta()) > 2.4)
//			continue;

		unsigned int n_lJet = ( tjet.pt() > tjet2.pt() ) ? jetIndex : jet2Index;
		LorentzVector lJet = theJetsCols[jetCollectionIndex]->ptrAt( n_lJet )->p4();
		double bdis_lJet = theJetsCols[jetCollectionIndex]->ptrAt( n_lJet )->bDiscriminator(bTagType);
		int flavour_lJet = theJetsCols[jetCollectionIndex]->ptrAt( n_lJet )->partonFlavour();
		double DR_lpho_lJet = tools_.DeltaR(lJet, diphoCand->leadingPhoton()->p4());
		double DR_spho_lJet = tools_.DeltaR(lJet, diphoCand->subLeadingPhoton()->p4());

		unsigned int n_sJet = ( tjet.pt() > tjet2.pt() ) ? jet2Index : jetIndex;
		LorentzVector sJet = theJetsCols[jetCollectionIndex]->ptrAt( n_sJet )->p4();
		double bdis_sJet = theJetsCols[jetCollectionIndex]->ptrAt( n_sJet )->bDiscriminator(bTagType);
		int flavour_sJet = theJetsCols[jetCollectionIndex]->ptrAt( n_sJet )->partonFlavour();
		double DR_lpho_sJet = tools_.DeltaR(sJet, diphoCand->leadingPhoton()->p4());
		double DR_spho_sJet = tools_.DeltaR(sJet, diphoCand->subLeadingPhoton()->p4());

		double minDRpho_ljet = min(DR_lpho_lJet, DR_spho_lJet);
		double minDRpho_sjet = min(DR_lpho_sJet, DR_spho_sJet);
		double minDRpho_jet = min(minDRpho_ljet, minDRpho_sjet);

		LorentzVector dijet = tjet+tjet2;
		double DR_jets = tools_.DeltaR(sJet, lJet);
		double DR_dipho_dijet = tools_.DeltaR(dijet, diphoCand->p4());

		leadingjets.push_back(lJet);
		subleadingjets.push_back(sJet);
		dijets.push_back(dijet);

		leadingjets_bDiscriminant.push_back(bdis_lJet);
		subleadingjets_bDiscriminant.push_back(bdis_sJet);

		leadingjets_partonID.push_back(flavour_lJet);
		subleadingjets_partonID.push_back(flavour_sJet);

		leadingJets_DRleadingPho.push_back(DR_lpho_lJet);
		leadingJets_DRsubleadingPho.push_back(DR_spho_lJet);
		subleadingJets_DRleadingPho.push_back(DR_lpho_sJet);
		subleadingJets_DRsubleadingPho.push_back(DR_spho_sJet);

		Jets_minDRPho.push_back(minDRpho_jet);

		Jets_DR.push_back(DR_jets);
		DiJets_DRDiPho.push_back(DR_dipho_dijet);
	    }
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
	tree->Branch("diphotonCandidate", &diphotonCandidate);
        tree->Branch("leadingJets", &leadingjets);
        tree->Branch("subleadingJets", &subleadingjets);
        tree->Branch("diJets", &dijets);
        tree->Branch("leadingJets_bDiscriminant", &leadingjets_bDiscriminant);
        tree->Branch("subleadingJets_bDiscriminant", &subleadingjets_bDiscriminant);
	tree->Branch("leadingJets_partonID", &leadingjets_partonID);
	tree->Branch("subleadingJets_partonID", &subleadingjets_partonID);
	tree->Branch("leadingJets_DRleadingPho", &leadingJets_DRleadingPho);
	tree->Branch("leadingJets_DRsubleadingPho", &leadingJets_DRsubleadingPho);
	tree->Branch("subleadingJets_DRleadingPho", &subleadingJets_DRleadingPho);
	tree->Branch("subleadingJets_DRsubleadingPho", &subleadingJets_DRsubleadingPho);
	tree->Branch("Jets_minDRPho", &Jets_minDRPho);
	tree->Branch("Jets_DR", &Jets_DR);
	tree->Branch("DiJets_DRDiPho", &DiJets_DRDiPho);
//	tree->Branch("nvtx", &nvtx);
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
        tree->Branch("leadingJet_KF", &leadingJet_KF);
        tree->Branch("leadingJet_Reg", &leadingJet_Reg);
        tree->Branch("leadingJet_RegKF", &leadingJet_RegKF);
        tree->Branch("leadingJet_bDis", &leadingJet_bDis, "leadingJet_bDis/F");
	tree->Branch("leadingJet_flavour", &leadingJet_flavour, "leadingJet_flavour/I");
        tree->Branch("subleadingJet", &subleadingJet);
        tree->Branch("subleadingJet_KF", &subleadingJet_KF);
        tree->Branch("subleadingJet_Reg", &subleadingJet_Reg);
        tree->Branch("subleadingJet_RegKF", &subleadingJet_RegKF);
        tree->Branch("subleadingJet_bDis", &subleadingJet_bDis, "subleadingJet_bDis/F");
	tree->Branch("subleadingJet_flavour", &subleadingJet_flavour, "subleadingJet_flavour/I");
        tree->Branch("dijetCandidate", &dijetCandidate);
        tree->Branch("dijetCandidate_KF", &dijetCandidate_KF);
        tree->Branch("dijetCandidate_Reg", &dijetCandidate_Reg);
        tree->Branch("dijetCandidate_RegKF", &dijetCandidate_RegKF);
        tree->Branch("diHiggsCandidate", &diHiggsCandidate);
        tree->Branch("diHiggsCandidate_KF", &diHiggsCandidate_KF);
        tree->Branch("diHiggsCandidate_Reg", &diHiggsCandidate_Reg);
        tree->Branch("diHiggsCandidate_RegKF", &diHiggsCandidate_RegKF);
	tree->Branch("isSignal", &isSignal, "isSignal/I");
	tree->Branch("isPhotonCR", &isPhotonCR, "isPhotonCR/I");
        tree->Branch("jet1PtRes", &jet1PtRes, "jet1PtRes/F");
        tree->Branch("jet1EtaRes", &jet1EtaRes, "jet1EtaRes/F");
        tree->Branch("jet1PhiRes", &jet1PhiRes, "jet1PhiRes/F");
        tree->Branch("jet2PtRes", &jet2PtRes, "jet2PtRes/F");
        tree->Branch("jet2EtaRes", &jet2EtaRes, "jet2EtaRes/F");
        tree->Branch("jet2PhiRes", &jet2PhiRes, "jet2PhiRes/F");
	tree->Branch("CosThetaStar", &CosThetaStar, "CosThetaStar/D");
		
    }
    std::map<std::string, std::string> replacements;
    globVar_->bookTreeVariables(tree, replacements);

    if(DEBUG) std::cout << "[bbggTree::beginJob] Output tree set!" << std::endl;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
    bbggTree::endJob() 
{
    outFile->cd();
//    tree->Write();
    outFile->Write();
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
