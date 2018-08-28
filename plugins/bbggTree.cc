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
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

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
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"

//Local
#include "flashgg/bbggTools/interface/bbggTools.h"
#include "flashgg/bbggTools/interface/bbggMC.h"
#include "flashgg/bbggTools/interface/bbggKinFit.h"
#include "flashgg/bbggTools/interface/bbggJetRegression.h"
#include "flashgg/bbggTools/interface/bbggJetSystematics.h"
#include "flashgg/bbggTools/interface/bbggPhotonCorrector.h"
#include "flashgg/bbggTools/interface/bbggNonResMVA.h"
#include "flashgg/bbggTools/interface/NonResWeights.h"
#include "flashgg/bbggTools/interface/genericMVA.h"

//New MVA
#include "flashgg/bbggTools/interface/bbggNonResMVA2017.h"

//needed for sigmaMOverM                                                                                                                                                                                                                  
#include "flashgg/Taggers/interface/FunctionHelpers.h"

//needed for diphoton mva                                                                                                                                                                                                                 
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"


//
// class declaration
//
/*
#ifdef _CINT_
#pragma link C++ class std::vector<std::pair<flashgg::Jet,flashgg::Jet>>;
#pragma link C++ class std::vector<std::pair<flashgg::Jet,flashgg::Jet>>;
#endif
*/

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
    bbggJetSystematics jetSys_;
    bbggPhotonCorrector phoCorr_;
    bbggNonResMVA nonresMVA_;
    bbggNonResMVA2017 nonresMVA2017_;
    bbggNonResMVA resMVA_;
    genericMVA ttHMVA_;
    flashgg::GlobalVariablesDumper* globVar_;
    //Parameter tokens
    edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diPhotonToken_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > genToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > genToken2_;
    std::vector<edm::InputTag> inputTagJets_;
    std::vector<edm::EDGetTokenT<edm::View<flashgg::Jet> > > tokenJets_;
    edm::InputTag genInfo_;
    edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
    edm::EDGetTokenT<edm::TriggerResults> triggerToken_;

    edm::EDGetTokenT<edm::View<reco::Vertex> > vertexToken_;
    edm::EDGetTokenT<edm::View<flashgg::Met> > METToken_;
/// ttH variables
    edm::EDGetTokenT<double> rhoToken_;
    edm::EDGetTokenT<edm::View<flashgg::Electron> > electronToken_;
    edm::EDGetTokenT<edm::View<flashgg::Muon> > muonToken_;



  //needed for diphoton mva
    edm::EDGetTokenT<edm::View<flashgg::DiPhotonMVAResult> > mvaResultToken_;


    //Efficiency histogram
    TH1F* h_Efficiencies = new TH1F("h_Efficiencies", "Efficiencies;Cut Level;Number of events", 10, 0, 10);
    TH1F* h_hggid = new TH1F("h_hggid", "", 100, -1, 1);
    //Tree objects
    LorentzVector leadingPhoton, subleadingPhoton, diphotonCandidate;
    LorentzVector leadingJet, subleadingJet, dijetCandidate;
    LorentzVector leadingJet_VBF, subleadingJet_VBF, DijetVBF;
    LorentzVector leadingMuon, subleadingMuon, leadingElectron, subleadingElectron; 
    LorentzVector p4MET;
    float dEta_VBF, Mjj_VBF;
    LorentzVector leadingJet_KF, subleadingJet_KF, dijetCandidate_KF;
    LorentzVector leadingJet_Reg, subleadingJet_Reg, dijetCandidate_Reg;
    LorentzVector leadingJet_RegKF, subleadingJet_RegKF, dijetCandidate_RegKF;
    LorentzVector diHiggsCandidate, diHiggsCandidate_KF, diHiggsCandidate_Reg,diHiggsCandidate_RegKF;
    vector<int> leadingPhotonID, leadingPhotonISO, subleadingPhotonID, subleadingPhotonISO;
    vector<double> genWeights;
    float leadingJet_bDis, subleadingJet_bDis, jet1PtRes, jet1EtaRes, jet1PhiRes, jet2PtRes, jet2EtaRes, jet2PhiRes;
    float leadingJet_CSVv2, leadingJet_cMVA, subleadingJet_CSVv2, subleadingJet_cMVA;
    float leadingPhotonIDMVA, subleadingPhotonIDMVA, DiJetDiPho_DR, PhoJetMinDr;
    float leadingPhotonSigOverE, subleadingPhotonSigOverE, sigmaMOverM, sigmaMOverMDecorr;
    float diphoMVA;
    float CosThetaStar, CosThetaStar_CS, CosTheta_bb, CosTheta_gg, CosTheta_bbgg, CosTheta_ggbb, Phi0, Phi1;
    std::map<std::string, int> myTriggerResults;
    float leadingPhotonR9full5x5, subleadingPhotonR9full5x5, customLeadingPhotonMVA, customSubLeadingPhotonMVA;
    int leadingPhotonHasGain1, leadingPhotonHasGain6;
    int subLeadingPhotonHasGain1, subLeadingPhotonHasGain6;
    float HHTagger, HHTagger2017, HHTagger_LM, HHTagger_HM, HHTagger2017_transform;
    float ResHHTagger, ResHHTagger_LM, ResHHTagger_HM;
    float ttHTagger;
    float MX, sumEt;
    TGraph *HHTagger2017_cumulative;
    
    double genTotalWeight;
    unsigned int nPromptInDiPhoton;
    int leadingPhotonEVeto, subleadingPhotonEVeto;
    int leadingJet_flavour, subleadingJet_flavour;
    int leadingJet_hadFlavour, subleadingJet_hadFlavour;
    int isSignal, isPhotonCR;
    //    int nvtx;

    Double_t gen_mHH;
    Double_t gen_cosTheta;
  Double_t gen_NRW;
  // -- End of Tree objects --
  // --    ---        --
    float leadingJet_genPtb, leadingJet_genPartonidb, leadingJet_genFlavourb,  leadingJet_genPartonFlavourb, leadingJet_genHadronFlavourb, leadingJet_genNbHadronsb, leadingJet_genNcHadronsb;
    float subleadingJet_genPtb, subleadingJet_genPartonidb, subleadingJet_genFlavourb,  subleadingJet_genPartonFlavourb, subleadingJet_genHadronFlavourb, subleadingJet_genNbHadronsb, subleadingJet_genNcHadronsb;

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
    unsigned int DoMVAPhotonID;
    edm::FileInPath bRegFileLeading, bRegFileSubLeading;
    unsigned int is2016, doCustomPhotonMVA;
    unsigned int doTnP;
    int doPhotonScale, doPhotonExtraScale, doPhotonSmearing;
    std::string PhotonScaleFile;
    int addNonResMVA, addNonResMVA2017, addttHMVA;
    edm::FileInPath NonResMVAWeights_LowMass, NonResMVAWeights_HighMass;
    edm::FileInPath ResMVAWeights_LowMass, ResMVAWeights_HighMass;
    edm::FileInPath NonResMVA2017Weights;

    edm::FileInPath ttHMVAWeights;
    edm::FileInPath NonResMVA2017Transformation;

    std::vector<std::string> NonResMVAVars;
    std::vector<std::string> NonResMVA2017Vars;
    std::vector<std::string> ttHMVAVars;

    double muPtThreshold, muEtaThreshold, muPFIsoSumRelThreshold;// deltaRMuonPhoThreshold;
    double dRPhoLeptonThreshold, dRJetLeptonThreshold;

    double elecPtThreshold;
    bool useElecMVARecipe, useElecLooseId;
    std::vector<double> elecEtaThresholds;

  int njets;
    vector<float> Xtt;
    float Xtt0, Xtt1, MjjW0, MjjW1, Mjjbt0, Mjjbt1;

    int jetSmear;
    int jetScale;
    std::string randomLabel;
    edm::FileInPath resFile, sfFile, scaleFile;

    //sigmaMoverM
    unsigned int  doDecorr;
    unsigned int  def_doDecorr;
    edm::FileInPath sigmaMdecorrFile;
    edm::FileInPath def_sigmaMdecorrFile;
    DecorrTransform* transfEBEB_;
    DecorrTransform* transfNotEBEB_;
    TH2D* h_decorrEBEB_;
    TH2D* h_decorrNotEBEB_;



  Bool_t getNonResGenInfo;
  // Class for NonRes re-weighting:
  NonResWeights *NRW;
  // Array to store the weights of 12 benchmarks and [0] is always 1:
  float  NRWeights[13];
  unsigned int BenchNum;

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
genToken2_( consumes<edm::View<reco::GenParticle> >( iConfig.getUntrackedParameter<edm::InputTag>( "GenTag2", edm::InputTag( "flashggPrunedGenParticles" ) ) ) ),
inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) )
{
  iConfig.dump();

    //now do what ever initialization is needed
    tools_ = bbggTools();
    jetReg_ = bbggJetRegression();
    jetSys_ = bbggJetSystematics();
    phoCorr_ = bbggPhotonCorrector();
    nonresMVA_ = bbggNonResMVA();
    resMVA_ = bbggNonResMVA();
    ttHMVA_ = genericMVA();

    NRW    = new NonResWeights();
    globVar_ = new flashgg::GlobalVariablesDumper(iConfig, consumesCollector() );
    //Lumi weight
    double lumiWeight_ = ( iConfig.getParameter<double>( "lumiWeight" ) );
    globVar_->dumpLumiFactor(lumiWeight_);
    EvtCount = 0;
    //Default values for thresholds
    std::string def_bTagType, def_PhotonMVAEstimator;
    unsigned int def_is2016 = 1;
    unsigned int def_doTnP = 0;
    std::vector<double> def_ph_pt, def_ph_eta, def_ph_r9;
    std::vector<double> def_diph_pt, def_diph_eta, def_diph_mass, def_MVAPhotonID;
    std::vector<double> def_jt_pt, def_jt_eta, def_jt_drPho, def_jt_bDis;
    std::vector<double> def_dijt_pt, def_dijt_eta, def_dijt_mass, def_cand_pt, def_cand_eta, def_cand_mass, def_dr_cands;
    std::vector<int> def_ph_elVeto, def_ph_doelVeto, def_ph_doID;
    std::vector<int> def_ph_doISO;
    std::vector<int> def_jt_doPU, def_jt_doID;
    std::vector<std::string> def_ph_whichID, def_ph_whichISO;
    unsigned int def_diph_onlyfirst;
    unsigned int def_n_bJets, def_doCustomPhotonMVA;
    int def_doPhotonScale, def_doPhotonExtraScale, def_doPhotonSmearing;
    std::string def_PhotonScaleFile;
    int def_addNonResMVA, def_addNonResMVA2017, def_addttHMVA;
    edm::FileInPath def_NonResMVA2017Weights;
    
    std::vector<std::string> def_NonResMVAVars, def_NonResMVA2017Vars, def_ttHMVAVars;

    edm::FileInPath def_bRegFileLeading, def_bRegFileSubLeading;

    //init values
    def_ph_pt.push_back(10.); def_ph_pt.push_back(10.); def_ph_eta.push_back(20.); def_ph_eta.push_back(20.);
    def_ph_r9.push_back(-1.); def_ph_r9.push_back(-1.); def_ph_elVeto.push_back(-1); def_ph_elVeto.push_back(-1);
    def_ph_doelVeto.push_back(0); def_ph_elVeto.push_back(0); def_ph_doID.push_back(0); def_ph_doID.push_back(0);
    def_ph_whichID.push_back("loose"); def_ph_whichID.push_back("loose"); def_ph_doISO.push_back(0); def_ph_doISO.push_back(0);
    def_ph_whichISO.push_back("loose"); def_ph_whichISO.push_back("loose"); def_MVAPhotonID.push_back(1.); def_MVAPhotonID.push_back(1.);
    def_diph_pt.push_back(10.); def_diph_pt.push_back(10.); def_diph_eta.push_back(0.); def_diph_eta.push_back(0.);
    def_diph_mass.push_back(0.); def_diph_mass.push_back(1000.); def_diph_onlyfirst = 0; def_jt_pt.push_back(10.);
    def_jt_pt.push_back(10.); def_jt_eta.push_back(20.); def_jt_eta.push_back(20.); def_jt_bDis.push_back(0.);
    def_jt_bDis.push_back(0.); def_jt_doPU.push_back(0); def_jt_doPU.push_back(0); def_jt_doID.push_back(0);
    def_jt_doID.push_back(0); def_n_bJets = 0; def_dijt_pt.push_back(10.); def_dijt_pt.push_back(10.);
    def_dijt_eta.push_back(20.); def_dijt_eta.push_back(20.); def_dijt_mass.push_back(0.); def_dijt_mass.push_back(1000.);
    def_jt_drPho.push_back(0.5); def_cand_pt.push_back(0.); def_cand_eta.push_back(20.); def_cand_mass.push_back(0.);
    def_cand_mass.push_back(2000.); def_dr_cands.push_back(0.11); def_doCustomPhotonMVA = 1;
    def_doPhotonScale = -10; def_doPhotonExtraScale = -10; def_doPhotonSmearing = -10;
    def_PhotonScaleFile = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/80X_ichepV2_2016_pho";

    unsigned int def_nPromptPhotons = 0;
    unsigned int def_doDoubleCountingMitigation = 0;
    unsigned int def_doPhotonCR = 0;
    unsigned int def_doJetRegression = 0;
    unsigned int def_DoMVAPhotonID = 0;
    int def_jetSmear = 0;
    int def_jetScale = 0;
    std::string def_randomLabel = "";
    def_addNonResMVA = 0;
    def_addNonResMVA2017 = 0;
    def_addttHMVA = 0;

    def_NonResMVAVars.push_back("");
    def_NonResMVA2017Vars.push_back("");
    def_ttHMVAVars.push_back("");

    edm::FileInPath def_resFile = edm::FileInPath("flashgg/Systematics/data/JER/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt");
    edm::FileInPath def_sfFile = edm::FileInPath("flashgg/Systematics/data/JER/Spring16_25nsV10_MC_SF_AK4PFchs.txt");
    edm::FileInPath def_scaleFile = edm::FileInPath("flashgg/bbggTools/data/JetSystematics/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt");

    std::vector<std::string> def_myTriggers;

    std::string def_fileName;

    def_bTagType = "pfCombinedInclusiveSecondaryVertexV2BJetTags";
    def_fileName =  "out.root";

    def_bRegFileLeading = edm::FileInPath("flashgg/bbggTools/data/BRegression/2016/BDTG_16plus3_jetGenJet_nu_80X_leading_9_13.weights.xml");
    def_bRegFileSubLeading = edm::FileInPath("flashgg/bbggTools/data/BRegression/2016/BDTG_16plus3_jetGenJet_nu_80X_trailing_9_13.weights.xml");

      //sigmaMoverM
      def_doDecorr=0;
      def_sigmaMdecorrFile = edm::FileInPath("flashgg/Taggers/data/diphoMVA_sigmaMoMdecorr_split_Mgg40_180.root");


    //Get photon ID thresholds from config file
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

    is2016 = iConfig.getUntrackedParameter<unsigned int>("is2016", def_is2016);
    doTnP = iConfig.getUntrackedParameter<unsigned int>("doTnP",def_doTnP);

    BenchNum = iConfig.getUntrackedParameter<unsigned int>("benchmark", 0);

    //photon selection parameters
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
    nPromptPhotons = iConfig.getUntrackedParameter<unsigned int>("nPromptPhotons", def_nPromptPhotons);
    doDoubleCountingMitigation = iConfig.getUntrackedParameter<unsigned int>("doDoubleCountingMitigation", def_doDoubleCountingMitigation);
    doPhotonCR = iConfig.getUntrackedParameter<unsigned int>("doPhotonCR", def_doPhotonCR);
    DoMVAPhotonID = iConfig.getUntrackedParameter<unsigned int>("DoMVAPhotonID", def_DoMVAPhotonID);
    MVAPhotonID = iConfig.getUntrackedParameter<std::vector<double>>("MVAPhotonID", def_MVAPhotonID);
    PhotonMVAEstimator = iConfig.getUntrackedParameter<std::string>("PhotonMVAEstimator", def_PhotonMVAEstimator);
    doPhotonScale = iConfig.getUntrackedParameter<int>("doPhotonScale", def_doPhotonScale);
    doPhotonExtraScale = iConfig.getUntrackedParameter<int>("doPhotonExtraScale", def_doPhotonExtraScale);
    doPhotonSmearing = iConfig.getUntrackedParameter<int>("doPhotonSmearing", def_doPhotonSmearing);
    PhotonScaleFile = iConfig.getUntrackedParameter<std::string >("PhotonCorrectionFile", def_PhotonScaleFile);
    doCustomPhotonMVA = iConfig.getUntrackedParameter<unsigned int>("doCustomPhotonMVA", def_doCustomPhotonMVA);

    //jet selection parameters
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
    doJetRegression = iConfig.getUntrackedParameter<unsigned int>("doJetRegression", def_doJetRegression);
    bTagType = iConfig.getUntrackedParameter<std::string>( "bTagType", def_bTagType );
    bRegFileLeading = iConfig.getUntrackedParameter<edm::FileInPath>("bRegFileLeading", def_bRegFileLeading);
    bRegFileSubLeading = iConfig.getUntrackedParameter<edm::FileInPath>("bRegFileSubLeading", def_bRegFileSubLeading);
    jetSmear = iConfig.getUntrackedParameter<int>("jetSmear", def_jetSmear);
    jetScale = iConfig.getUntrackedParameter<int>("jetScale", def_jetScale);

    //extra selection
    cand_pt 	= iConfig.getUntrackedParameter<std::vector<double > >("CandidatePt", def_cand_mass);
    cand_eta 	= iConfig.getUntrackedParameter<std::vector<double > >("CandidateEta", def_cand_mass);
    cand_mass = iConfig.getUntrackedParameter<std::vector<double > >("CandidateMassWindow", def_cand_mass);
    dr_cands  = iConfig.getUntrackedParameter<std::vector<double > >("CandidatesDeltaR", def_dr_cands);

    //leptons selection
    muPtThreshold = iConfig.getParameter<double>("muPtThreshold");
    muEtaThreshold = iConfig.getParameter<double>("muEtaThreshold");
    muPFIsoSumRelThreshold = iConfig.getParameter<double>("muPFIsoSumRelThreshold");

    dRPhoLeptonThreshold = iConfig.getParameter<double>("dRPhoLeptonThreshold");
    dRJetLeptonThreshold = iConfig.getParameter<double>("dRJetLeptonThreshold");

    elecPtThreshold  = iConfig.getParameter<double>("elecPtThreshold");
    useElecMVARecipe = iConfig.getParameter<bool>("useElecMVARecipe"); 
    useElecLooseId = iConfig.getParameter<bool>("useElecLooseId");
    elecEtaThresholds = iConfig.getParameter<std::vector<double > >("elecEtaThresholds");
    // Xtt0                MW0                 Mt0                Xtt1                 MW1               Mt1
    Xtt.push_back(1000); Xtt.push_back(0); Xtt.push_back(0); Xtt.push_back(1000); Xtt.push_back(0); Xtt.push_back(0); 

    fileName = iConfig.getUntrackedParameter<std::string>( "OutFileName", def_fileName );

     //sigmaMOverM
    doDecorr = iConfig.getUntrackedParameter<unsigned int>("doSigmaMdecorr", def_doDecorr);
    sigmaMdecorrFile = iConfig.getUntrackedParameter<edm::FileInPath>("sigmaMdecorrFile", def_sigmaMdecorrFile); 


    addNonResMVA = iConfig.getUntrackedParameter<unsigned int>("addNonResMVA", def_addNonResMVA);
    addNonResMVA2017 = iConfig.getUntrackedParameter<unsigned int>("addNonResMVA2017", def_addNonResMVA2017);
    addttHMVA = iConfig.getUntrackedParameter<unsigned int>("addttHMVA", def_addttHMVA);

    NonResMVAWeights_LowMass = iConfig.getUntrackedParameter<edm::FileInPath>("NonResMVAWeights_LowMass");
    NonResMVAWeights_HighMass = iConfig.getUntrackedParameter<edm::FileInPath>("NonResMVAWeights_HighMass");
    ResMVAWeights_LowMass = iConfig.getUntrackedParameter<edm::FileInPath>("ResMVAWeights_LowMass");
    ResMVAWeights_HighMass = iConfig.getUntrackedParameter<edm::FileInPath>("ResMVAWeights_HighMass");
    NonResMVAVars = iConfig.getUntrackedParameter<std::vector<std::string > >("NonResMVAVars", NonResMVAVars);

    NonResMVA2017Weights = iConfig.getUntrackedParameter<edm::FileInPath>("NonResMVA2017Weights");
    NonResMVA2017Vars = iConfig.getUntrackedParameter<std::vector<std::string > >("NonResMVA2017Vars", NonResMVA2017Vars);
    NonResMVA2017Transformation = iConfig.getUntrackedParameter<edm::FileInPath>("NonResMVA2017Transformation");

    ttHMVAWeights = iConfig.getUntrackedParameter<edm::FileInPath>("ttHMVAWeights");
    ttHMVAVars = iConfig.getUntrackedParameter<std::vector<std::string > >("ttHMVAVars", def_ttHMVAVars);

    //tokens and labels
    genInfo_ = iConfig.getUntrackedParameter<edm::InputTag>( "genInfo", edm::InputTag("generator") );
    genInfoToken_ = consumes<GenEventInfoProduct>( genInfo_ );
    myTriggers = iConfig.getUntrackedParameter<std::vector<std::string> >("myTriggers", def_myTriggers);
    triggerToken_ = consumes<edm::TriggerResults>( iConfig.getParameter<edm::InputTag>( "triggerTag" ) );

    if (is2016){
      triggerToken_ = consumes<edm::TriggerResults>( edm::InputTag("TriggerResults", "", "HLT" )) ;
    }

    rhoToken_ = consumes<double>( iConfig.getParameter<edm::InputTag>( "rhoTag" ) );
    vertexToken_ = consumes<edm::View<reco::Vertex> >( iConfig.getParameter<edm::InputTag> ( "VertexTag" ) );
    METToken_ = consumes<edm::View<flashgg::Met> >( iConfig.getParameter<edm::InputTag> ( "metTag" ) );
    electronToken_ = consumes<edm::View<flashgg::Electron> >( iConfig.getParameter<edm::InputTag> ( "ElectronTag" ) );
    muonToken_ = consumes<edm::View<flashgg::Muon> >( iConfig.getParameter<edm::InputTag>( "MuonTag" ) );


    randomLabel = iConfig.getUntrackedParameter<std::string>("randomLabel", def_randomLabel);

    resFile = iConfig.getUntrackedParameter<edm::FileInPath>("resFile", def_resFile);
    sfFile = iConfig.getUntrackedParameter<edm::FileInPath>("sfFile", def_sfFile);
    scaleFile = iConfig.getUntrackedParameter<edm::FileInPath>("scaleFile", def_scaleFile);

    getNonResGenInfo = iConfig.getUntrackedParameter<bool>("getNonResGenInfo", false);

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

    mvaResultToken_ = consumes<edm::View<flashgg::DiPhotonMVAResult> >( iConfig.getParameter<edm::InputTag> ( "MVAResultTag" ) ) ;

    kinFit_ = bbggKinFit();
    kinFit_.SetJetResolutionParameters(etaBins, ptRes, etaRes, phiRes);
    jetSys_.SetupScale(scaleFile.fullPath().data());

    for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
        auto token = consumes<edm::View<flashgg::Jet> >(inputTagJets_[i]);
        tokenJets_.push_back(token);
    }

    if(doJetRegression) {
        std::cout << "Regressing jets with the following weights files: " << std::endl;
        std::cout << "\t Leading jet: " << bRegFileLeading.fullPath().data() << std::endl;
        std::cout << "\t SubLeading jet: " << bRegFileSubLeading.fullPath().data() << std::endl;
        jetReg_.SetupRegression("BDTG method",  bRegFileLeading.fullPath().data(), bRegFileSubLeading.fullPath().data());
    }

    //Prepare scales and smearings
    if (doPhotonSmearing > -10 || doPhotonScale > -10) {
        phoCorr_.SetupCorrector(PhotonScaleFile);
        phoCorr_.setRandomLabel(std::string("rnd_g_E"));
    }
    
    if(addttHMVA){
      std::cout << "Adding ttH MVA with files: " << ttHMVAWeights.fullPath() << std::endl;
      ttHMVA_.setupMVA( ttHMVAWeights.fullPath().data(), ttHMVAVars);
    }

    if(addNonResMVA) {
        std::cout << "Adding nonres MVA with files: " << NonResMVAWeights_LowMass.fullPath() << " and " << NonResMVAWeights_HighMass.fullPath() << std::endl;
        nonresMVA_.SetupNonResMVA( NonResMVAWeights_LowMass.fullPath().data(), NonResMVAWeights_HighMass.fullPath().data(), NonResMVAVars);
        resMVA_.SetupNonResMVA( ResMVAWeights_LowMass.fullPath().data(), ResMVAWeights_HighMass.fullPath().data(), NonResMVAVars);
    }

   if(addNonResMVA2017) {
        std::cout << "Adding nonres MVA 207 with files: " << NonResMVA2017Weights.fullPath() <<  std::endl;
 	 nonresMVA2017_.SetupNonResMVA( NonResMVA2017Weights.fullPath().data(), NonResMVA2017Vars);
        TFile* fHHTagger2017Transformation = new TFile((NonResMVA2017Transformation.fullPath()).c_str(), "READ");
        HHTagger2017_cumulative = (TGraph*)fHHTagger2017Transformation->Get("cumulativeGraph"); 
      }
 
      if (doDecorr){
        TFile* f_decorr = new TFile((sigmaMdecorrFile.fullPath()).c_str(), "READ");
        h_decorrEBEB_ = (TH2D*)f_decorr->Get("hist_sigmaM_M_EBEB"); 
        h_decorrNotEBEB_ = (TH2D*)f_decorr->Get("hist_sigmaM_M_notEBEB");
        if(h_decorrEBEB_ && h_decorrNotEBEB_){
 	 transfEBEB_ = new DecorrTransform(h_decorrEBEB_ , 125., 1, 0);
 	 transfNotEBEB_ = new DecorrTransform(h_decorrNotEBEB_ , 125., 1, 0);
        }
        else {
 	 throw cms::Exception( "Configuration" ) << "The file "<<sigmaMdecorrFile.fullPath()<<" provided for sigmaM/M decorrelation does not contain the expected histograms."<<std::endl;
        }
      }

      

    if (getNonResGenInfo){
      std::string fileNameWei1 = edm::FileInPath("flashgg/bbggTools/data/NonResReWeight/weights_v1_1507_points.root").fullPath();
      std::string fileNameWei2 = edm::FileInPath("flashgg/bbggTools/data/NonResReWeight/weights_v3_bench12_points.root").fullPath();
      NRW->LoadHists(fileNameWei1, fileNameWei2);
      for (UInt_t n=0; n<13; n++) 
	NRWeights[n]=1;

      // Check that provided benchmark number is in range 1-12
      if (BenchNum>12){
	std::cout<<"\t ** warning** Provided BenchNum is out of range (1-12): "<<BenchNum<<std::endl;
	std::cout<<"I'm gonna set it to 0 to avoid problems"<<std::endl;
	BenchNum=0;
      }
    }

    std::cout << "Parameters initialized... \n ############ Doing selection tree!" <<  std::endl;

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

    globVar_->fill(iEvent);

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
    CosThetaStar_CS = -999;
    CosTheta_bb = -999;
    CosTheta_gg = -999;
    CosTheta_bbgg = -999;
    CosTheta_ggbb = -999;
    Phi0 = -999;
    Phi1 = -999;
    genTotalWeight = 1.0;
    leadingPhotonR9full5x5 = -999;
    subleadingPhotonR9full5x5 = -999;
    customLeadingPhotonMVA = -999;
    customSubLeadingPhotonMVA = -999;
    leadingPhotonHasGain1 = -10;
    leadingPhotonHasGain6 = -10;
    subLeadingPhotonHasGain1 = -10;
    subLeadingPhotonHasGain6 = -10;
    HHTagger = -10;
    HHTagger2017 = -10;
    HHTagger2017_transform = -10;
    HHTagger_LM = -10;
    HHTagger_HM = -10;
    ResHHTagger = -10;
    ResHHTagger_LM = -10;
    ResHHTagger_HM = -10;
    ttHTagger = -10;


    dEta_VBF = -999; 
    Mjj_VBF = 0;

    sumEt = 0;

    njets = 0;

    leadingMuon.SetPxPyPzE(0,0,0,0);          leadingElectron.SetPxPyPzE(0,0,0,0);
    subleadingMuon.SetPxPyPzE(0,0,0,0);    subleadingElectron.SetPxPyPzE(0,0,0,0);


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
    leadingJet_CSVv2 = 0;
    leadingJet_cMVA = 0;
    subleadingJet_CSVv2 = 0;
    subleadingJet_cMVA = 0;

    p4MET.SetPxPyPzE(0,0,0,0);// MET p4

    leadingJet_VBF.SetPxPyPzE(0,0,0,0);// = LeadingJet->p4();
    subleadingJet_VBF.SetPxPyPzE(0,0,0,0);// = LeadingJet->p4();

    

    gen_mHH = 0;
    gen_cosTheta = -99;

    gen_NRW = 1;

      //gen jets info
      leadingJet_genPtb=-999; leadingJet_genPartonidb=-999; leadingJet_genFlavourb=-999;  leadingJet_genPartonFlavourb=-999; leadingJet_genHadronFlavourb=-999; leadingJet_genNbHadronsb=-999; leadingJet_genNcHadronsb=-999;
      subleadingJet_genPtb=-999; subleadingJet_genPartonidb=-999; subleadingJet_genFlavourb=-999; subleadingJet_genPartonFlavourb=-999; subleadingJet_genHadronFlavourb=-999; subleadingJet_genNbHadronsb=-999; subleadingJet_genNcHadronsb=-999;


    //Get Jets collections!
    JetCollectionVector theJetsCols( inputTagJets_.size() );
    for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
        iEvent.getByToken( tokenJets_[j], theJetsCols[j] );
    }

    if (DEBUG) std::cout << "Number of jet collections!!!! " << theJetsCols.size() << std::endl;

    Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
    iEvent.getByToken( diPhotonToken_, diPhotons );


    const double rhoFixedGrd = globVar_->valueOf(globVar_->indexOf("rho"));
    tools_.setRho(rhoFixedGrd);
    const double nPVs = globVar_->valueOf(globVar_->indexOf("nvtx"));
    const int runNumber = globVar_->valueOf(globVar_->indexOf("run"));

    Handle<View<pat::PackedGenParticle> > genParticles;
    Handle<View<reco::GenParticle> > genParticles2;

    //Trigger
    if(myTriggers.size() > 0 && !is2016){
        Handle<edm::TriggerResults> trigResults;
        iEvent.getByToken(triggerToken_, trigResults);
        const edm::TriggerNames &names = iEvent.triggerNames(*trigResults);
        myTriggerResults = tools_.TriggerSelection(myTriggers, names, trigResults);
    }

    //MC Weights
    Handle<GenEventInfoProduct> genInfo;
    if( ! iEvent.isRealData() ) {
        iEvent.getByToken(genInfoToken_, genInfo);
        genTotalWeight = genInfo->weight();
        iEvent.getByToken( genToken_, genParticles);
        iEvent.getByToken( genToken2_, genParticles2);
    } else {
        genTotalWeight = 1;
    }

    //Remove MC duplicates (if needed)
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
            if (nPromptPhotons > 1 && nPrompt > 1) diphoVec.push_back(dipho);
            if (nPromptPhotons < 2 && nPrompt < 2) diphoVec.push_back(dipho);
        } else {
            diphoVec.push_back(dipho);
        }
    }

    //Eff step 0
    if(DEBUG) std::cout << "Number of diphoton candidates: " << diphoVec.size() << std::endl;
    h_Efficiencies->Fill(0.0, genTotalWeight);


    if (getNonResGenInfo){
      
      // ---- Gen HH info
      // Here we get gen level mHH and costheta for non-resonant samples
      TLorentzVector H1, H2;
      UInt_t nH = 0;

      for( unsigned int igen = 0; igen < genParticles2->size(); igen++) {
	edm::Ptr<reco::GenParticle> genPar = genParticles2->ptrAt(igen);

	if (genPar->pdgId()==25 && genPar->isHardProcess()){

	  if (nH==0)
	    H1.SetXYZM(genPar->px(), genPar->py(), genPar->pz(), genPar->mass());
	  if (nH==1)
	    H2.SetXYZM(genPar->px(), genPar->py(), genPar->pz(), genPar->mass());

	  nH++;
	  if (nH==2) break;
	}
      }

      if (nH==2){

	  gen_mHH  = (H1+H2).M();
	  gen_cosTheta = tools_.getCosThetaStar_CS(H1,H2,6500);

	  // Now, lets fill in the weigts for the 12 benchmarks.
	  for (UInt_t n=1; n<13; n++)
	    NRWeights[n] = NRW->GetWeight(1506+n, gen_mHH, gen_cosTheta);
	  
	  gen_NRW =  NRWeights[BenchNum];
      }

      // --- End of gen HH info
      
    }

    //Eff step 0
    if(DEBUG) std::cout << "Number of diphoton candidates: " << diphoVec.size() << std::endl;
    h_Efficiencies->Fill(0.0, genTotalWeight*gen_NRW);

    if (diphoVec.size() < 1) return;
    //Eff step 1
    h_Efficiencies->Fill(1, genTotalWeight*gen_NRW);
    if (theJetsCols.size() < 1) return;

    bool cutsChecked = tools_.CheckCuts();
    if(!cutsChecked) {
        std::cout << "You haven't filled all the cuts correctly!" << std::endl;
        return;
    }

    if(DEBUG) std::cout << "[bbggTree::analyze] About to do event selection! " << std::endl;

    /////////////////////////////////////////////////////
    /// DOING SELECTION HERE NOW ////////////////////////
    /////////////////////////////////////////////////////

    /////////////
    // Photons //
    /////////////
 
    //Prepare reference collection
    std::vector<flashgg::DiPhotonCandidate> diphotonCollection;
    for( unsigned int dpIndex = 0; dpIndex < diphoVec.size(); dpIndex++ )
    {
        edm::Ptr<flashgg::DiPhotonCandidate> thisDPPtr = diphoVec[ dpIndex ];
        flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(thisDPPtr.get());
        diphotonCollection.push_back(*thisDPPointer);
        diphotonCollection[dpIndex].makePhotonsPersistent();
    }

    ////Do scale and smearing
    //Smear MC
    if (!iEvent.isRealData() && doPhotonSmearing > -10) {
        phoCorr_.setVariation(doPhotonSmearing);
        phoCorr_.SmearPhotonsInDiPhotons(diphotonCollection, runNumber);
    }
    //Scale Data
    if (iEvent.isRealData() && doPhotonScale > -10) {
        phoCorr_.setVariation(doPhotonScale);
        phoCorr_.ScalePhotonsInDiPhotons(diphotonCollection, runNumber);
    }
    //Extra scale on data
    if (iEvent.isRealData() && doPhotonExtraScale > -10) {
        phoCorr_.setVariation(doPhotonScale);
        phoCorr_.ExtraScalePhotonsInDiPhotons(diphotonCollection);
    }

    //Trigger preselection on diphoton candidates:
    vector<flashgg::DiPhotonCandidate> PreSelDipho;
    if(is2016) PreSelDipho = tools_.DiPhotonPreselection( diphotonCollection );
    if(!is2016) PreSelDipho = tools_.DiPhoton76XPreselection( diphotonCollection, myTriggerResults);
    if(doTnP && is2016){
       Handle<edm::TriggerResults> trigResults;
       iEvent.getByToken(triggerToken_, trigResults);
	const edm::TriggerNames &names = iEvent.triggerNames(*trigResults);
       myTriggerResults = tools_.TriggerSelection(myTriggers, names, trigResults);
       PreSelDipho = tools_.DiPhotonPreselectionTnP2016( diphotonCollection, myTriggerResults);
     }
    if(DEBUG) std::cout << "[bbggTree::analyze] Number of pre-selected diphotons: " << PreSelDipho.size() << std::endl;
    //If no diphoton passed presel, skip event
    if ( PreSelDipho.size() < 1 ) return;
    h_Efficiencies->Fill(2, genTotalWeight*gen_NRW);

    
    //Kinematic selection
    std::vector<flashgg::DiPhotonCandidate> KinDiPhoton;
    if(!doTnP)  KinDiPhoton  = tools_.DiPhotonKinematicSelection( PreSelDipho, 1);
    else {
      KinDiPhoton = PreSelDipho;
    }
    //    if(DEBUG)
    if( KinDiPhoton.size() < 1) return;
    h_Efficiencies->Fill(3, genTotalWeight*gen_NRW);

    //Evaluate photon IDs
    vector<pair<flashgg::DiPhotonCandidate, int> > KinDiPhotonWithID;
      if(!doTnP) KinDiPhotonWithID = tools_.EvaluatePhotonIDs( KinDiPhoton, doCustomPhotonMVA );
      else KinDiPhotonWithID = tools_.EvaluatePhotonIDs( KinDiPhoton, doCustomPhotonMVA, 1 );
    vector<flashgg::DiPhotonCandidate> SignalDiPhotons = tools_.GetDiPhotonsInCategory( KinDiPhotonWithID, 2 );
    vector<flashgg::DiPhotonCandidate> CRDiPhotons = tools_.GetDiPhotonsInCategory( KinDiPhotonWithID, 1 );

    //needed for diphoton mva
    if(tools_.indexSel_>-1){
      Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;
      iEvent.getByToken( mvaResultToken_, mvaResults );
      edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt(tools_.indexSel_);//get mva of selected diphoton
      diphoMVA=mvares->result;
    }else{
      diphoMVA=-999;
    }

    //Select diphoton candidate
    flashgg::DiPhotonCandidate diphoCandidate;
    if(SignalDiPhotons.size() > 0) {
        diphoCandidate = tools_.PtSumDiPhotonSelection(SignalDiPhotons);//SignalDiPhotons[0];
        isSignal = 1;
        isPhotonCR = 0;
    }
    if(SignalDiPhotons.size() < 1 && doPhotonCR && CRDiPhotons.size() > 0) {
        diphoCandidate = tools_.PtSumDiPhotonSelection(CRDiPhotons);//CRDiPhotons[0];
        isSignal = 0;
        isPhotonCR = 1;
    }
    if(SignalDiPhotons.size() < 1 && CRDiPhotons.size() < 1) return; //if event is not signal and is not photon control region, skip
    if(SignalDiPhotons.size() < 1 && !doPhotonCR) return; //if event is not signal and you don't want to save photon control region, skip

    if(isSignal) h_Efficiencies->Fill(4, genTotalWeight*gen_NRW);
    if(DEBUG) std::cout << "[bbggTree::analyze] Is signal region: " << isSignal << "; Is control region: " << isPhotonCR << std::endl;

    ///////////// JETS
    bool hasLeadJet = 0;
    bool hasSubJet = 0;
    flashgg::Jet LeadingJet, SubLeadingJet;
    unsigned int jetCollectionIndex = diphoCandidate.jetCollectionIndex();
    if( theJetsCols[jetCollectionIndex]->size() < 2) return;
    if(isSignal) h_Efficiencies->Fill(5, genTotalWeight*gen_NRW);

    std::vector<flashgg::Jet> testCollection;
    for( unsigned int jetIndex = 0; jetIndex < theJetsCols[jetCollectionIndex]->size(); jetIndex++ )
    {
        edm::Ptr<flashgg::Jet> thisJetPtr = theJetsCols[jetCollectionIndex]->ptrAt( jetIndex );
        flashgg::Jet * thisJetPointer = const_cast<flashgg::Jet *>(thisJetPtr.get());
        testCollection.push_back(*thisJetPointer);
    }
    std::vector<std::pair<flashgg::Jet,flashgg::Jet>> KinJets;
    std::vector<flashgg::Jet> SelJets;
    std::vector<flashgg::Jet> SelVBFJets;
    //Here I can apply smearing to my jets
    if(jetSmear!=0){
      jetSys_.SetupSmear(resFile.fullPath().data(), sfFile.fullPath().data());
      jetSys_.SmearJets(testCollection, rhoFixedGrd, jetSmear, randomLabel);
    }
    if(jetScale!=0){
        jetSys_.ScaleJets(testCollection, jetScale);
    }

    //Make dijet candidates
    std::vector<std::pair<flashgg::Jet,flashgg::Jet>> myDiJetCandidates;
    for (unsigned int j1 = 0; j1 < testCollection.size(); j1++) {
      for (unsigned int j2 = j1+1; j2 < testCollection.size(); j2++) {
        if ( testCollection[j1].pt() > testCollection[j2].pt() ) myDiJetCandidates.push_back(std::make_pair(testCollection[j1], testCollection[j2]));
        else myDiJetCandidates.push_back(std::make_pair(testCollection[j2], testCollection[j1]));
      }
    }

    std::vector<flashgg::Jet> collectionForVBF;
    std::vector<flashgg::Jet> collectionForCounting;
    for( unsigned int jetIndex = 0; jetIndex < testCollection.size(); jetIndex++ ){
      collectionForVBF.push_back(testCollection[jetIndex]);
      collectionForCounting.push_back(testCollection[jetIndex]);
    }


    Handle<View<reco::Vertex> > vertices;
    iEvent.getByToken( vertexToken_, vertices );

    //Regression bzns
    if(doJetRegression!=0) {
//        Handle<View<pat::MET> > METs;
        Handle<View<flashgg::Met> > METs;
        iEvent.getByToken( METToken_, METs );
	if( METs->size() != 1 )
        { std::cout << "WARNING number of MET is not equal to 1" << std::endl; }
        Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );
//        Ptr<pat::MET> theMET = METs->ptrAt( 0 );

	p4MET = theMET->p4();
        if(DEBUG) std::cout << "DOING REGRESSION! JetCol before: " << testCollection.size() << std::endl;
        jetReg_.RegressedJets(myDiJetCandidates, nPVs, p4MET );

        if(DEBUG) std::cout << "DOING REGRESSION! JetCol after: " << testCollection.size() << std::endl;
    }

//    if(jetSmear!=0 || jetScale!=0 || doJetRegression!=0){
    if(DEBUG) std::cout << "DOING SELECTION AFTER JET MANIPULATION - PreSel" << std::endl;

    KinJets = tools_.JetPreSelection(myDiJetCandidates, diphoCandidate);
    if( KinJets.size() < 1 ) return;
    if(isSignal) h_Efficiencies->Fill(6, genTotalWeight*gen_NRW);

    if(DEBUG) std::cout << "DOING SELECTION AFTER JET MANIPULATION - KinSel" << KinJets.size() << std::endl;
    SelJets = tools_.DiJetSelection(KinJets, 1);
    if( SelJets.size() < 1 ) return;
    if(isSignal) h_Efficiencies->Fill(7, genTotalWeight*gen_NRW);

    if(DEBUG) std::cout << "DOING SELECTION AFTER JET MANIPULATION - DiJetSel" << SelJets.size() << std::endl;

    if( SelJets.size() > 1 ) {
        hasLeadJet = 1;
        hasSubJet = 1;
        LeadingJet = SelJets[0];
        SubLeadingJet = SelJets[1];
    }


    /////////////////////////////////////////////////////
    /// DOING SELECTION HERE NOW ////////////////////////
    /////////////////////////////////////////////////////


    if(!hasLeadJet || !hasSubJet) return;

    flashgg::DiPhotonCandidate diphoCand = diphoCandidate;//tools_.GetSelected_diphoCandidate();

    if(DEBUG && doJetRegression){
        std::cout << "Userflots of leading jet: " << std::endl;
        std::cout << "\t nSecVertices: " << LeadingJet.userFloat("nSecVertices") << std::endl;
        std::cout << "\t vtxNTracks: " << LeadingJet.userFloat("vtxNTracks") << std::endl;
        std::cout << "\t vtxMass: " << LeadingJet.userFloat("vtxMass") << std::endl;
        std::cout << "\t vtxPx: " << LeadingJet.userFloat("vtxPx") << std::endl;
        std::cout << "\t vtxPy: " << LeadingJet.userFloat("vtxPy") << std::endl;
        std::cout << "\t vtx3DVal: " << LeadingJet.userFloat("vtx3DVal") << std::endl;
        std::cout << "\t vtx3DSig: " << LeadingJet.userFloat("vtx3DSig") << std::endl;
        std::cout << "\t leadTrackPt: " << LeadingJet.userFloat("leadTrackPt") << std::endl;
        std::cout << "\t softLepPt: " << LeadingJet.userFloat("softLepPt") << std::endl;
        std::cout << "\t softLepRatio: " << LeadingJet.userFloat("softLepRatio") << std::endl;
        std::cout << "\t softLepDr: " << LeadingJet.userFloat("softLepDr") << std::endl;
    }

    if (DEBUG) std::cout << "Jet collection picked: " << diphoCand.jetCollectionIndex() << std::endl;

    nPromptInDiPhoton = 999;
    if ( genInfo.isValid() ){
        bbggMC _mcTools = bbggMC();
        nPromptInDiPhoton = _mcTools.CheckNumberOfPromptPhotons(diphoCand, genParticles);
    }

    if ( iEvent.isRealData() ) {
        leadingPhotonHasGain1 = diphoCand.leadingPhoton()->checkStatusFlag( flashgg::Photon::kHasSwitchToGain1 );
        leadingPhotonHasGain6 = diphoCand.leadingPhoton()->checkStatusFlag( flashgg::Photon::kHasSwitchToGain6 );
        subLeadingPhotonHasGain1 = diphoCand.subLeadingPhoton()->checkStatusFlag( flashgg::Photon::kHasSwitchToGain1 );
        subLeadingPhotonHasGain6 = diphoCand.subLeadingPhoton()->checkStatusFlag( flashgg::Photon::kHasSwitchToGain6 );

    }

    leadingPhotonR9full5x5 = diphoCand.leadingPhoton()->full5x5_r9();
    subleadingPhotonR9full5x5 = diphoCand.subLeadingPhoton()->full5x5_r9();

    //    diphotonCandidate = diphoCand.p4(); //this is wrong since diphoton collection is not updated by the smearer
    diphotonCandidate = (diphoCand.leadingPhoton()->p4()+diphoCand.subLeadingPhoton()->p4());
    leadingPhoton = diphoCand.leadingPhoton()->p4();
    subleadingPhoton = diphoCand.subLeadingPhoton()->p4();

    leadingJet = LeadingJet.p4();
    leadingJet_bDis = LeadingJet.bDiscriminator(bTagType);
    leadingJet_flavour = LeadingJet.partonFlavour();
    leadingJet_hadFlavour = LeadingJet.hadronFlavour();
    subleadingJet = SubLeadingJet.p4();
    subleadingJet_bDis = SubLeadingJet.bDiscriminator(bTagType);
    subleadingJet_flavour = SubLeadingJet.partonFlavour();
    subleadingJet_hadFlavour = SubLeadingJet.hadronFlavour();

    leadingJet_CSVv2 = LeadingJet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    leadingJet_cMVA = LeadingJet.bDiscriminator("pfCombinedMVAV2BJetTags");
    subleadingJet_CSVv2 = SubLeadingJet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    subleadingJet_cMVA = SubLeadingJet.bDiscriminator("pfCombinedMVAV2BJetTags");;

     //..... gen jets info 
 
     int cflavLeading = 0; //~correct flavour definition
     int cflavSubLeading = 0; //~correct flavour definition
     if ( !iEvent.isRealData() ) {
       int hflav = LeadingJet.hadronFlavour();//4 if c, 5 if b, 0 if light jets
       int pflav = LeadingJet.partonFlavour();
 
 	if( hflav != 0 ) {
 	  cflavLeading = hflav;
 	} else { //not a heavy jet
 	  cflavLeading = std::abs(pflav) == 4 || std::abs(pflav) == 5 ? 0 : pflav;
 	}
 	
 	hflav = SubLeadingJet.hadronFlavour();
 	pflav = SubLeadingJet.partonFlavour();
 
 	if( hflav != 0 ) {
 	  cflavSubLeading = hflav;
 	} else { //not a heavy jet
 	  cflavSubLeading = std::abs(pflav) == 4 || std::abs(pflav) == 5 ? 0 : pflav;
 	}
 
 
       leadingJet_genPtb = ( LeadingJet.genJet()!=0 ? LeadingJet.genJet()->pt() : -1. );
       leadingJet_genPartonidb = LeadingJet.genParton() ? LeadingJet.genParton()->pdgId() : 0 ;
       leadingJet_genFlavourb = cflavLeading;
       leadingJet_genPartonFlavourb = LeadingJet.partonFlavour();
       leadingJet_genHadronFlavourb = LeadingJet.hadronFlavour();
       leadingJet_genNbHadronsb = LeadingJet.jetFlavourInfo().getbHadrons().size();
       leadingJet_genNcHadronsb = LeadingJet.jetFlavourInfo().getcHadrons().size();
       subleadingJet_genPtb = ( SubLeadingJet.genJet()!=0 ? SubLeadingJet.genJet()->pt() : -1. );
       subleadingJet_genPartonidb = SubLeadingJet.genParton() ? SubLeadingJet.genParton()->pdgId() : 0 ;
       subleadingJet_genFlavourb = cflavSubLeading;
       subleadingJet_genPartonFlavourb = SubLeadingJet.partonFlavour();
       subleadingJet_genHadronFlavourb = SubLeadingJet.hadronFlavour();
       subleadingJet_genNbHadronsb = SubLeadingJet.jetFlavourInfo().getbHadrons().size();
       subleadingJet_genNcHadronsb = SubLeadingJet.jetFlavourInfo().getcHadrons().size();
       
     }//geninfo



    dijetCandidate = leadingJet + subleadingJet;
    diHiggsCandidate = diphotonCandidate + dijetCandidate;

    MX = diHiggsCandidate.M() - diphotonCandidate.M() - dijetCandidate.M() + 250;

 
   // ================================ VBF definition ======================

    collectionForVBF = tools_.JetVBFPreSelection(collectionForVBF, diphoCandidate, SelJets);

    if (collectionForVBF.size() > 1) 
      SelVBFJets = tools_.DiJetVBFSelection(collectionForVBF, SelJets);


    if (SelVBFJets.size() > 1){
      leadingJet_VBF = SelVBFJets[0].p4();
      subleadingJet_VBF = SelVBFJets[1].p4();
      DijetVBF = leadingJet_VBF + subleadingJet_VBF;
      
      dEta_VBF = fabs(leadingJet_VBF.eta() - subleadingJet_VBF.eta());
      Mjj_VBF = DijetVBF.M();
    }


    // ==================== Kinematic fit =======

    if(DEBUG) std::cout << "[bbggTree::analyze] Doing kinematic fit!" << std::endl;
    kinFit_.KinematicFit(leadingJet, subleadingJet);
    leadingJet_KF = kinFit_.GetJet1();
    subleadingJet_KF = kinFit_.GetJet2();
    dijetCandidate_KF = leadingJet_KF + subleadingJet_KF;
    diHiggsCandidate_KF = dijetCandidate_KF + diphotonCandidate;
    jet1PtRes = kinFit_.GetPtResolution(leadingJet);
    jet1EtaRes = kinFit_.GetEtaResolution(leadingJet);
    jet1PhiRes= kinFit_.GetPhiResolution(leadingJet);

    // ===================== ANGLES ==============

    vector<float> CosThetaAngles  = tools_.CosThetaAngles(&diphoCand, LeadingJet, SubLeadingJet);
    TLorentzVector djc, dpc;
    djc.SetPtEtaPhiE( dijetCandidate.pt(), dijetCandidate.eta(), dijetCandidate.phi(), dijetCandidate.energy());
    dpc.SetPtEtaPhiE( diphotonCandidate.pt(), diphotonCandidate.eta(), diphotonCandidate.phi(), diphotonCandidate.energy());
    CosThetaStar_CS = tools_.getCosThetaStar_CS(djc,dpc,6500); //Colin Sopper Frame CTS
    CosThetaStar = CosThetaAngles[0];
    CosTheta_bb = CosThetaAngles[1];
    CosTheta_gg = CosThetaAngles[2];
    CosTheta_bbgg = CosThetaAngles[3];
    CosTheta_ggbb = CosThetaAngles[4];
    vector<double> PhiAngles = tools_.getPhi(&diphoCand, LeadingJet, SubLeadingJet);
    Phi0 = PhiAngles[0];
    Phi1 = PhiAngles[1];

    DiJetDiPho_DR = tools_.DeltaR(diphotonCandidate, dijetCandidate);

    PhoJetMinDr = min( min( tools_.DeltaR( leadingPhoton, leadingJet ), tools_.DeltaR( leadingPhoton, subleadingJet ) ),
                       min( tools_.DeltaR( subleadingPhoton, leadingJet ), tools_.DeltaR( subleadingPhoton, subleadingJet ) ) );


    // ================== PHOTON MVA ============
    
    customLeadingPhotonMVA = diphoCand.leadingPhoton()->phoIdMvaDWrtVtx( diphoCand.vtx() );
    h_hggid->Fill(diphoCand.leadingPhoton()->phoIdMvaDWrtVtx( diphoCand.vtx() ));
    customSubLeadingPhotonMVA = diphoCand.subLeadingPhoton()->phoIdMvaDWrtVtx( diphoCand.vtx() );
    leadingPhotonIDMVA = diphoCand.leadingPhoton()->userFloat(PhotonMVAEstimator);
    subleadingPhotonIDMVA = diphoCand.subLeadingPhoton()->userFloat(PhotonMVAEstimator);
     leadingPhotonSigOverE = diphoCand.leadingPhoton()->sigEOverE();
     subleadingPhotonSigOverE = diphoCand.subLeadingPhoton()->sigEOverE();
     sigmaMOverM = 0.5 * TMath::Sqrt(leadingPhotonSigOverE *leadingPhotonSigOverE + subleadingPhotonSigOverE * subleadingPhotonSigOverE );


    if(DEBUG) std::cout << "customLeadingPhotonMVA: " << diphoCand.leadingPhoton()->phoIdMvaDWrtVtx( diphoCand.vtx() ) << std::endl;
    if(DEBUG) std::cout << "leadingPhotonIDMVA: " << diphoCand.leadingPhoton()->userFloat(PhotonMVAEstimator) << std::endl;

    //Fill HHTagger
    if(addNonResMVA) {
        std::map<std::string, float> HHVars;
        HHVars["leadingJet_bDis"] = leadingJet_bDis;
        HHVars["subleadingJet_bDis"] = subleadingJet_bDis;
        HHVars["diphotonCandidate.Pt()/(diHiggsCandidate.M())"] = diphotonCandidate.Pt()/(diHiggsCandidate.M());
        HHVars["fabs(CosThetaStar_CS)"] = fabs(CosThetaStar_CS);
        HHVars["fabs(CosTheta_bb)"] = fabs(CosTheta_bb);
        HHVars["fabs(CosTheta_gg)"] = fabs(CosTheta_gg);
        HHVars["dijetCandidate.Pt()/(diHiggsCandidate.M())"] = dijetCandidate.Pt()/(diHiggsCandidate.M());
        HHVars["diphotonCandidate.Pt()/dijetCandidate.Pt()"] = diphotonCandidate.Pt()/dijetCandidate.Pt();
        std::vector<float> myHHTagger = nonresMVA_.mvaDiscriminants(HHVars);
        std::vector<float> myResHHTagger = resMVA_.mvaDiscriminants(HHVars);
        HHTagger_LM = myHHTagger[0];
        HHTagger_HM = myHHTagger[1];
        HHTagger = ( (diHiggsCandidate.M()-dijetCandidate.M()-diphotonCandidate.M()+250)<350 ) ? HHTagger_LM : HHTagger_HM;
        ResHHTagger_LM = myResHHTagger[0];
        ResHHTagger_HM = myResHHTagger[1];
        ResHHTagger = ( (diHiggsCandidate.M()-dijetCandidate.M()-diphotonCandidate.M()+250)<500 ) ? ResHHTagger_LM : ResHHTagger_HM;
    }


     //Fill HHTagger 2017
     if(addNonResMVA2017) {
         std::map<std::string, float> HHVars2017;
         HHVars2017["leadingJet_bDis"] = leadingJet_bDis;
         HHVars2017["subleadingJet_bDis"] = subleadingJet_bDis;
         HHVars2017["fabs(CosThetaStar_CS)"] = fabs(CosThetaStar_CS);
         HHVars2017["fabs(CosTheta_bb)"] = fabs(CosTheta_bb);
         HHVars2017["fabs(CosTheta_gg)"] = fabs(CosTheta_gg);
         HHVars2017["diphotonCandidate.Pt()/(diHiggsCandidate.M())"] = diphotonCandidate.Pt()/(diHiggsCandidate.M());
         HHVars2017["dijetCandidate.Pt()/(diHiggsCandidate.M())"] = dijetCandidate.Pt()/(diHiggsCandidate.M());
 	HHVars2017["customLeadingPhotonIDMVA"] = customLeadingPhotonMVA;
 	HHVars2017["customSubLeadingPhotonIDMVA"] = customSubLeadingPhotonMVA;
 	HHVars2017["leadingPhotonSigOverE"] = leadingPhotonSigOverE;
 	HHVars2017["subleadingPhotonSigOverE"] = subleadingPhotonSigOverE;
 	HHVars2017["sigmaMOverMDecorr"] = sigmaMOverMDecorr;
 	HHVars2017["PhoJetMinDr"] = PhoJetMinDr;
 
         std::vector<float> myHHTagger2017 = nonresMVA2017_.mvaDiscriminants(HHVars2017);
         HHTagger2017 = myHHTagger2017[0];

         HHTagger2017_transform = HHTagger2017_cumulative->Eval(HHTagger2017) ;
 
     }


     //sigmaMOverM
     if(doDecorr){
       //                std::cout<<"sigmaMdecorrFile is set, so we evaluate the transf"<<std::endl;
       double mass_sigma[2]={0.,0.};
       double dummy[1]={0.};
       mass_sigma[0]=diphotonCandidate.M();
       mass_sigma[1]=sigmaMOverM;
                 
       //splitting EBEB and !EBEB, using cuts as in preselection
       if(abs(diphoCand.leadingPhoton()->superCluster()->eta())<1.4442 && abs(diphoCand.subLeadingPhoton()->superCluster()->eta())<1.4442){
 	sigmaMOverMDecorr = (*transfEBEB_)(mass_sigma,dummy);
       }
       else{
 	sigmaMOverMDecorr = (*transfNotEBEB_)(mass_sigma,dummy);
       }
     }
 


    if(DEBUG) std::cout << "[bbggTree::analyze] After filling candidates" << std::endl;


    leadingPhotonEVeto = diphoCand.leadingPhoton()->passElectronVeto();
    subleadingPhotonEVeto = diphoCand.subLeadingPhoton()->passElectronVeto();



    // ============== ttH hadronic ==============
   
    collectionForCounting = tools_.JetPreSelection(collectionForCounting, diphoCandidate);
    njets = collectionForCounting.size();


    if (collectionForCounting.size() > 2 && SelJets.size() > 1){
      Xtt = tools_.XttCalculation(collectionForCounting, SelJets);
    }
    else {
      Xtt[0] = 1000, Xtt[1] = 0, Xtt[2] = 0, Xtt[3] = 1000, Xtt[4] = 0, Xtt[5] = 0;
    }

    Xtt0 = Xtt[0], Xtt1 = Xtt[3]; MjjW0 = Xtt[1],  MjjW1 = Xtt[4]; Mjjbt0 = Xtt[2],  Mjjbt1 = Xtt[5];

    for( unsigned int jetIndex = 0; jetIndex < collectionForVBF.size(); jetIndex++ )
      sumEt += collectionForVBF[jetIndex].pt();
    


    // ========== ttH leptonic ===========

    edm::Handle<double>  rho;
    iEvent.getByToken(rhoToken_,rho);


    double looseLeptonPtThreshold = 10;

    // =========== electrons =============

    Handle<View<flashgg::Electron> > theElectrons;
    iEvent.getByToken( electronToken_, theElectrons );

    std::vector<edm::Ptr<flashgg::Electron> > selectedElectrons = selectStdAllElectrons( theElectrons->ptrs(), vertices->ptrs(), looseLeptonPtThreshold, elecEtaThresholds, useElecMVARecipe, useElecLooseId, *rho, iEvent.isRealData() );

    std::vector<edm::Ptr<flashgg::Electron> > tagElectrons = tools_.filterElectrons( selectedElectrons, diphoCandidate, leadingJet, subleadingJet, dRPhoLeptonThreshold, dRJetLeptonThreshold);

    if (tagElectrons.size() > 0) leadingElectron = tagElectrons.at( 0 )->p4();
    if (tagElectrons.size() > 1) subleadingElectron = tagElectrons.at( 1 )->p4();

    // =========== muons =============

    Handle<View<flashgg::Muon> > theMuons;
    iEvent.getByToken( muonToken_, theMuons );
    //diphoCandidate
    std::vector<edm::Ptr<flashgg::Muon> > selectedMuons = selectAllMuons( theMuons->ptrs(), vertices->ptrs(), muEtaThreshold, looseLeptonPtThreshold, muPFIsoSumRelThreshold);
    std::vector<edm::Ptr<flashgg::Muon> > tagMuons = tools_.filterMuons( selectedMuons, diphoCandidate, leadingJet, subleadingJet, dRPhoLeptonThreshold, dRJetLeptonThreshold);

    if (tagMuons.size() > 0) leadingMuon = tagMuons.at( 0 )->p4();
    if (tagMuons.size() > 1) subleadingMuon = tagMuons.at( 1 )->p4();



   // ============== ttH BDT ==============

     if(addttHMVA) {
       std::map<std::string, float> ttHVars;
       ttHVars["sumEt"] = sumEt;
       ttHVars["MET"] = p4MET.Pt();

       ttHVars["dPhi1"] = deltaPhi(p4MET.Phi(), leadingJet.Phi());
       ttHVars["dPhi2"] = deltaPhi(p4MET.Phi(), subleadingJet.Phi());

       ttHVars["PhoJetMinDr"] =  PhoJetMinDr;

       ttHVars["njets>8"] =  njets > 8;
       ttHVars["Xtt0"] = Xtt0;
       ttHVars["Xtt1"] = Xtt1;

       ttHVars["pte1"] = leadingElectron.Pt();
       ttHVars["pte2"] = subleadingElectron.Pt();
       ttHVars["ptmu1"] = leadingMuon.Pt();
       ttHVars["ptmu2"] = subleadingMuon.Pt();

       ttHVars["fabs_CosThetaStar_CS"] = fabs(CosThetaStar_CS);
       ttHVars["fabs_CosTheta_bb"] = fabs(CosTheta_bb);
       
    
       ttHTagger = ttHMVA_.mvaDiscriminants(ttHVars);
     }
       
       if(DEBUG) std::cout << "GOT TO THE END!!" << std::endl;
       tree->Fill();
    
    if(DEBUG) std::cout << "Histograms filled!" << std::endl;

}


// ------------ method called once each job just before starting event loop  ------------
void
bbggTree::beginJob()
{
    if(DEBUG) std::cout << "[bbggTree::beginJob] Setting up output tree..." << std::endl;

    outFile = new TFile(fileName.c_str(), "RECREATE");
    tree = new TTree("bbggSelectionTree", "Flat tree for HH->bbgg analyses (after pre selection)");
    tree->Branch("genWeights", &genWeights);
    tree->Branch("genTotalWeight", &genTotalWeight, "genTotalWeight/D");
    tree->Branch("gen_mHH", &gen_mHH, "gen_mHH/D");
    tree->Branch("gen_cosTheta", &gen_cosTheta, "cosTheta/D");
    tree->Branch("gen_NRW", &gen_NRW, "NRW/D");
    tree->Branch("leadingPhoton", &leadingPhoton);
    tree->Branch("leadingPhotonID", &leadingPhotonID);
    tree->Branch("leadingPhotonISO", &leadingPhotonISO);
    tree->Branch("leadingPhotonEVeto", &leadingPhotonEVeto, "leadingPhotonEVeto/I");
    tree->Branch("leadingPhotonIDMVA", &leadingPhotonIDMVA, "leadingPhotonIDMVA/F");
    tree->Branch("customLeadingPhotonIDMVA", &customLeadingPhotonMVA, "customLeadingPhotonIDMVA/F");
    tree->Branch("leadingPhotonR9full5x5", &leadingPhotonR9full5x5, "leadingPhotonR9full5x5/F");
    tree->Branch("leadingPhotonHasGain1", &leadingPhotonHasGain1, "leadingPhotonHasGain1/I");
    tree->Branch("leadingPhotonHasGain6", &leadingPhotonHasGain6, "leadingPhotonHasGain6/I");
    tree->Branch("subleadingPhoton", &subleadingPhoton);
    tree->Branch("subleadingPhotonID", &subleadingPhotonID);
    tree->Branch("subleadingPhotonISO", &subleadingPhotonISO);
    tree->Branch("subleadingPhotonEVeto", &subleadingPhotonEVeto, "subleadingPhotonEVeto/I");
    tree->Branch("subleadingPhotonIDMVA", &subleadingPhotonIDMVA, "subleadingPhotonIDMVA/F");
    tree->Branch("customSubLeadingPhotonIDMVA", &customSubLeadingPhotonMVA, "customSubLeadingPhotonIDMVA/F");
    tree->Branch("subleadingPhotonR9full5x5", &subleadingPhotonR9full5x5, "subleadingPhotonR9full5x5/F");
    tree->Branch("subLeadingPhotonHasGain1", &subLeadingPhotonHasGain1, "subLeadingPhotonHasGain1/I");
    tree->Branch("subLeadingPhotonHasGain6", &subLeadingPhotonHasGain6, "subLeadingPhotonHasGain6/I");
    tree->Branch("leadingPhotonSigOverE", &leadingPhotonSigOverE, "leadingPhotonSigOverE/F");
    tree->Branch("subleadingPhotonSigOverE", &subleadingPhotonSigOverE, "subleadingPhotonSigOverE/F");
    tree->Branch("sigmaMOverM", &sigmaMOverM, "sigmaMOverM/F");
    tree->Branch("sigmaMOverMDecorr", &sigmaMOverMDecorr, "sigmaMOverMDecorr/F");
     tree->Branch("diphoMVA", &diphoMVA, "diphoMVA/F");
    tree->Branch("diphotonCandidate", &diphotonCandidate);
    tree->Branch("nPromptInDiPhoton", &nPromptInDiPhoton, "nPromptInDiPhoton/I");
    tree->Branch("leadingJet", &leadingJet);
    tree->Branch("leadingJet_KF", &leadingJet_KF);
    tree->Branch("leadingJet_Reg", &leadingJet_Reg);
    tree->Branch("leadingJet_RegKF", &leadingJet_RegKF);
    tree->Branch("leadingJet_bDis", &leadingJet_bDis, "leadingJet_bDis/F");
    tree->Branch("leadingJet_CSVv2", &leadingJet_CSVv2, "leadingJet_CSVv2/F");
    tree->Branch("leadingJet_cMVA", &leadingJet_cMVA, "leadingJet_cMVA/F");
    tree->Branch("leadingJet_flavour", &leadingJet_flavour, "leadingJet_flavour/I");
    tree->Branch("leadingJet_hadFlavour", &leadingJet_hadFlavour, "leadingJet_hadFlavour/I");
    tree->Branch("subleadingJet", &subleadingJet);
    tree->Branch("subleadingJet_KF", &subleadingJet_KF);
    tree->Branch("subleadingJet_Reg", &subleadingJet_Reg);
    tree->Branch("subleadingJet_RegKF", &subleadingJet_RegKF);
    tree->Branch("subleadingJet_bDis", &subleadingJet_bDis, "subleadingJet_bDis/F");
    tree->Branch("subleadingJet_CSVv2", &subleadingJet_CSVv2, "subleadingJet_CSVv2/F");
    tree->Branch("subleadingJet_cMVA", &subleadingJet_cMVA, "subleadingJet_cMVA/F");
    tree->Branch("subleadingJet_flavour", &subleadingJet_flavour, "subleadingJet_flavour/I");
    tree->Branch("subleadingJet_hadFlavour", &subleadingJet_hadFlavour, "subleadingJet_hadFlavour/I");
         //gen jets info
     tree->Branch("leadingJet_genPtb",&leadingJet_genPtb, "leadingJet_genPtb/F");
     tree->Branch("leadingJet_genPartonidb",&leadingJet_genPartonidb, "leadingJet_genPartonidb/F");
     tree->Branch("leadingJet_genFlavourb",&leadingJet_genFlavourb, "leadingJet_genFlavourb/F");
     tree->Branch("leadingJet_genPartonFlavourb",&leadingJet_genPartonFlavourb, "leadingJet_genPartonFlavourb/F");
     tree->Branch("leadingJet_genHadronFlavourb",&leadingJet_genHadronFlavourb, "leadingJet_genHadronFlavourb/F");
     tree->Branch("leadingJet_genNbHadronsb",&leadingJet_genNbHadronsb, "leadingJet_genNbHadronsb/F");
     tree->Branch("leadingJet_genNcHadronsb",&leadingJet_genNcHadronsb, "leadingJet_genNcHadronsb/F");
 
     tree->Branch("subleadingJet_genPtb",&subleadingJet_genPtb, "subleadingJet_genPtb/F");
     tree->Branch("subleadingJet_genPartonidb",&subleadingJet_genPartonidb, "subleadingJet_genPartonidb/F");
     tree->Branch("subleadingJet_genFlavourb",&subleadingJet_genFlavourb, "subleadingJet_genFlavourb/F");
     tree->Branch("subleadingJet_genPartonFlavourb",&subleadingJet_genPartonFlavourb, "subleadingJet_genPartonFlavourb/F");
     tree->Branch("subleadingJet_genHadronFlavourb",&subleadingJet_genHadronFlavourb, "subleadingJet_genHadronFlavourb/F");
     tree->Branch("subleadingJet_genNbHadronsb",&subleadingJet_genNbHadronsb, "subleadingJet_genNbHadronsb/F");
     tree->Branch("subleadingJet_genNcHadronsb",&subleadingJet_genNcHadronsb, "subleadingJet_genNcHadronsb/F");

    tree->Branch("HHTagger2017", &HHTagger2017, "HHTagger2017/F"); 
    tree->Branch("HHTagger2017_transform", &HHTagger2017_transform, "HHTagger2017_transform/F"); 

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
    tree->Branch("CosThetaStar", &CosThetaStar, "CosThetaStar/F");
    tree->Branch("CosThetaStar_CS", &CosThetaStar_CS, "CosThetaStar_CS/F");
    tree->Branch("CosTheta_bb", &CosTheta_bb, "CosTheta_bb/F");
    tree->Branch("CosTheta_gg", &CosTheta_gg, "CosTheta_gg/F");
    tree->Branch("CosTheta_bbgg", &CosTheta_bbgg, "CosTheta_bbgg/F");
    tree->Branch("CosTheta_ggbb", &CosTheta_ggbb, "CosTheta_ggbb/F");
    tree->Branch("Phi0", &Phi0, "Phi0/F");
    tree->Branch("Phi1", &Phi1, "Phi1/F");
    tree->Branch("TriggerResults", &myTriggerResults);
    tree->Branch("DiJetDiPho_DR", &DiJetDiPho_DR, "DiJetDiPho_DR/F");
    tree->Branch("PhoJetMinDr", &PhoJetMinDr, "PhoJetMinDr/F");
    tree->Branch("leadingJet_VBF", &leadingJet_VBF);
    tree->Branch("subleadingJet_VBF", &subleadingJet_VBF);
    tree->Branch("dEta_VBF", &dEta_VBF, "dEta_VBF/F");
    tree->Branch("Mjj_VBF", &Mjj_VBF, "Mjj_VBF/F");
    tree->Branch("HHTagger", &HHTagger, "HHTagger/F"); 
    tree->Branch("HHTagger_LM", &HHTagger_LM, "HHTagger_LM/F"); 
    tree->Branch("HHTagger_HM", &HHTagger_HM, "HHTagger_HM/F"); 
    tree->Branch("ResHHTagger", &ResHHTagger, "ResHHTagger/F"); 
    tree->Branch("ResHHTagger_LM", &ResHHTagger_LM, "ResHHTagger_LM/F"); 
    tree->Branch("ResHHTagger_HM", &ResHHTagger_HM, "ResHHTagger_HM/F"); 
    tree->Branch("MX", &MX, "MX/F"); 
    tree->Branch("MET", &p4MET);
    tree->Branch("sumEt", &sumEt, "sumEt/F");
    tree->Branch("njets", &njets, "njets/I");
    tree->Branch("Xtt0", &Xtt0, "Xtt0/F");
    tree->Branch("Xtt1", &Xtt1, "Xtt1/F");
    tree->Branch("MjjW0", &MjjW0, "MjjW0/F");
    tree->Branch("MjjW1", &MjjW1, "MjjW1/F");
    tree->Branch("Mjjbt0", &Mjjbt0, "Mjjbt0/F");
    tree->Branch("Mjjbt1", &Mjjbt1, "Mjjbt1/F");

    //    std::cout << "branch set Xtt" << std::endl;

    tree->Branch("leadingMuon", &leadingMuon);
    tree->Branch("subleadingMuon", &subleadingMuon);
    tree->Branch("leadingElectron", &leadingElectron);
    tree->Branch("subleadingElectron", &subleadingElectron);

    tree->Branch("ttHTagger", &ttHTagger, "ttHTagger/F"); 

    std::map<std::string, std::string> replacements;
    globVar_->bookTreeVariables(tree, replacements);

    std::cout << "done" << std::endl;

    if(DEBUG) std::cout << "[bbggTree::beginJob] Output tree set!" << std::endl;

}

// ------------ method called once each job just after ending the event loop  ------------
void
bbggTree::endJob()
{
outFile->cd();
//tree->Write();
h_Efficiencies->Write();
h_hggid->Write();
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
