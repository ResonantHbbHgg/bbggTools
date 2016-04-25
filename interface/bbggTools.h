#ifndef FLASHgg_bbggTools_h
#define FLASHgg_bbggTools_h

//FLASHgg files
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
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

using namespace std;

class bbggTools{
public:
	bbggTools() : hasDiPho( 0 ), hasLeadJet( 0 ), hasSubJet( 0 ), _isSignal( 0 ), _isPhotonCR( 0 ) { rho_ = -10;}
	~bbggTools() {}
	typedef math::XYZTLorentzVector LorentzVector;
    typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

	std::vector<edm::Ptr<flashgg::DiPhotonCandidate>> DiPhotonKinematicSelection(vector<edm::Ptr<flashgg::DiPhotonCandidate>> diphoCol, bool DoMassCut = 0);
	std::vector<edm::Ptr<flashgg::DiPhotonCandidate>> DiPhotonIDSelection( std::vector<edm::Ptr<flashgg::DiPhotonCandidate>> diphoCol);
    std::vector<edm::Ptr<flashgg::DiPhotonCandidate>> GetDiPhotonsInCategory( std::vector<std::pair<edm::Ptr<flashgg::DiPhotonCandidate>, int > > SelectedDiPhotons, int category );
    std::vector<std::pair<edm::Ptr<flashgg::DiPhotonCandidate>, int > > EvaluatePhotonIDs( std::vector<edm::Ptr<flashgg::DiPhotonCandidate>> diphoCol);
	std::vector<edm::Ptr<flashgg::Jet>> DiJetSelection(std::vector<edm::Ptr<flashgg::Jet>> Jets, bool DoMassCut = 0);
	std::vector<edm::Ptr<flashgg::Jet>> JetPreSelection(JetCollectionVector jetsCol, edm::Ptr<flashgg::DiPhotonCandidate> diphoCandidate);
    std::vector<edm::Ptr<flashgg::DiPhotonCandidate>> DiPhoton76XPreselection(vector<edm::Ptr<flashgg::DiPhotonCandidate>> diphoCol, std::map<std::string, int> myTriggersResults);

	edm::Ptr<flashgg::DiPhotonCandidate> MVAIDDiPhotonSelection( vector<edm::Ptr<flashgg::DiPhotonCandidate>> DiPhotons);
	edm::Ptr<flashgg::DiPhotonCandidate> PtSumDiPhotonSelection( vector<edm::Ptr<flashgg::DiPhotonCandidate>> DiPhotons);
    
    std::map<std::string,int> TriggerSelection(std::vector<std::string> myTriggers, const edm::TriggerNames &names, edm::Handle<edm::TriggerResults> triggerBits);


    std::map<int, vector<double> > getWhichID (std::string wpoint);
    std::map<int, vector<double> > getWhichISO (std::string wpoint);
    
	double getCHisoToCutValue(const flashgg::DiPhotonCandidate * dipho, int whichPho);
	double getCHisoToCutValue(edm::Ptr<flashgg::DiPhotonCandidate> dipho, int whichPho);
    
    double getNHisoToCutValue(const flashgg::Photon* pho);
	double getPHisoToCutValue(const flashgg::Photon* pho);
    
    double getNHisoToCutValue(const flashgg::Photon* pho, vector<double> nhCorr);
	double getPHisoToCutValue(const flashgg::Photon* pho, vector<double> phCorr);
    
	double getEA( float eta, int whichEA);
	double DeltaR( bbggTools::LorentzVector vec1, bbggTools::LorentzVector vec2);
    
	bool isPhoID(edm::Ptr<flashgg::Photon> pho, vector<double> cuts);
	bool isPhoID(const flashgg::Photon* pho, vector<double> cuts);
    
	bool isPhoISO(edm::Ptr<flashgg::DiPhotonCandidate> pho, int whichPho, vector<double> cuts, vector<double> nhCorr, vector<double> phCorr);
	bool isPhoISO(const flashgg::DiPhotonCandidate * pho, int whichPho, vector<double> cuts, vector<double> nhCorr, vector<double> phCorr);
    
    bool isJetID(edm::Ptr<flashgg::Jet> jet);
    
	void setRho(double rho) {rho_ = rho;}
    
	bool IsSignal() { return _isSignal; }
	bool IsPhotonCR() { return _isPhotonCR; }
    
	bool HasLeadJet(){ return hasLeadJet; }
	bool HasSubJet(){ return hasSubJet; }
	
	//Set cuts for selection
	void SetCut_PhotonPtOverDiPhotonMass( vector<double> cuts) { _PhotonPtOverDiPhotonMass = cuts; }
	void SetCut_PhotonEta( vector<double> cuts) { _PhotonEta = cuts; }
	void SetCut_PhotonR9( vector<double> cuts) { _PhotonR9 = cuts; }
	void SetCut_PhotonElectronVeto( vector<int> cuts) { _PhotonElectronVeto = cuts; }
	void SetCut_PhotonDoElectronVeto( vector<int> cuts) { _PhotonDoElectronVeto = cuts; }
	void SetCut_DiPhotonPt( vector<double> cuts) { _DiPhotonPt = cuts; }
	void SetCut_DiPhotonEta( vector<double> cuts) { _DiPhotonEta = cuts; }
	void SetCut_DiPhotonMassWindow( vector<double> cuts) { _DiPhotonMassWindow = cuts; }
	void SetCut_DiPhotonOnlyFirst( unsigned int cuts) { _DiPhotonOnlyFirst = cuts; }
    
	void SetCut_JetPt( vector<double> cuts) { _JetPt = cuts; }
	void SetCut_JetEta( vector<double> cuts) { _JetEta = cuts; }
	void SetCut_JetBDiscriminant( vector<double> cuts) { _JetBDiscriminant = cuts; }
	void SetCut_JetDrPho( vector<double> cuts ) { _JetDrPho = cuts; }
	void SetCut_JetDoPUID( vector<int> cuts) { _JetDoPUID = cuts; }
	void SetCut_JetDoID( vector<int> cuts) { _JetDoID = cuts; }
	void SetCut_n_bJets( unsigned int cuts) { _n_bJets = cuts; }
	void SetCut_bTagType( std::string cuts) { _bTagType = cuts; }

	void SetCut_DiJetPt( vector<double> cuts) { _DiJetPt = cuts; }
	void SetCut_DiJetEta( vector<double> cuts) { _DiJetEta = cuts; }
	void SetCut_DiJetMassWindow( vector<double> cuts) { _DiJetMassWindow = cuts; }
    
	void SetCut_CandidateMassWindow( vector<double> cuts) { _CandidateMassWindow = cuts; }
	void SetCut_CandidatePt( vector<double> cuts) { _CandidatePt = cuts; }
	void SetCut_CandidateEta( vector<double> cuts) { _CandidateEta = cuts; }
	void SetCut_CandidatesDeltaR( vector<double> cuts) { _CandidatesDeltaR = cuts; }
    
	void SetCut_PhotonDoID( vector<int> cuts) { _PhotonDoID = cuts; }
	void SetCut_PhotonDoISO( vector<int> cuts) { _PhotonDoISO = cuts; }
    void SetCut_PhotonDoEVeto( vector<int> cuts) { _PhotonDoEVeto = cuts; }
    
    void SetCut_phoIDloose( std::map<int, vector<double> > cuts ) { _phoIDloose = cuts; }
    void SetCut_phoIDmedium( std::map<int, vector<double> > cuts ) { _phoIDmedium  = cuts; }
    void SetCut_phoIDtight( std::map<int, vector<double> > cuts ) { _phoIDtight = cuts; }
    void SetCut_phoISOloose( std::map<int, vector<double> > cuts ) { _phoISOloose = cuts; }
    void SetCut_phoISOmedium( std::map<int, vector<double> > cuts ) { _phoISOmedium = cuts; }
    void SetCut_phoISOtight( std::map<int, vector<double> > cuts ) { _phoISOtight = cuts; }
    void SetCut_nhCorr( std::map<int, vector<double> > cuts ) { _nhCorr = cuts; }
    void SetCut_phCorr( std::map<int, vector<double> > cuts ) { _phCorr = cuts; }
    void SetCut_phoWhichID( vector<std::string> cuts) { _phoWhichID = cuts; }
    void SetCut_phoWhichISO( vector<std::string> cuts) { _phoWhichISO = cuts; }
	
	void SetPhotonCR( unsigned int cuts ) { _doPhotonCR = cuts; } 
    
    void SetCut_DoMVAPhotonID(unsigned int cuts){ _DoMVAPhotonID = cuts; } 
    void SetCut_MVAPhotonID(vector<double> cuts){ _MVAPhotonID = cuts; } 
    void SetCut_PhotonMVAEstimator( std::string cuts){ _PhotonMVAEstimator = cuts;}
    		
	//Perform event selection
	bool AnalysisSelection( vector<edm::Ptr<flashgg::DiPhotonCandidate>> diphoCol, 
							JetCollectionVector jetsCol );


	//Get selected objects
	edm::Ptr<flashgg::DiPhotonCandidate> GetSelected_diphoCandidate();
	edm::Ptr<flashgg::Jet> GetSelected_leadingJetCandidate();
	edm::Ptr<flashgg::Jet> GetSelected_subleadingJetCandidate();
	
	//Check cuts
	bool CheckCuts();

private:
	double rho_;
	
	//Selected objects
	edm::Ptr<flashgg::DiPhotonCandidate> diphoCandidate;
	bool hasDiPho;
	edm::Ptr<flashgg::Jet> leadingJetCandidate;
	bool hasLeadJet;
	edm::Ptr<flashgg::Jet> subleadingJetCandidate;
	bool hasSubJet;
	bool _isSignal;
	bool _isPhotonCR;
	
	//Cut values
	unsigned int _doPhotonCR;
	vector<double> _PhotonPtOverDiPhotonMass;
	vector<double> _PhotonEta;
	vector<double> _PhotonR9;
	vector<int> _PhotonElectronVeto;
	vector<int> _PhotonDoElectronVeto;
	vector<double> _DiPhotonPt;
	vector<double> _DiPhotonEta;
	vector<double> _DiPhotonMassWindow;
	unsigned int _DiPhotonOnlyFirst;
	vector<double> _JetPt;
	vector<double> _JetEta;
	vector<double> _JetBDiscriminant;
	vector<double> _JetDrPho;
	vector<int> _JetDoPUID;
	vector<int> _JetDoID;
	unsigned int _n_bJets;
	vector<double> _DiJetPt;
	vector<double> _DiJetEta;
	vector<double> _DiJetMassWindow;
	vector<double> _CandidateMassWindow;
	vector<double> _CandidatePt;
	vector<double> _CandidateEta;
	std::string _bTagType;
	vector<double> _CandidatesDeltaR;
    vector<double> _MVAPhotonID;
    
	vector<int> _PhotonDoID;
	vector<int> _PhotonDoISO;
    vector<int> _PhotonDoEVeto;
    unsigned int _DoMVAPhotonID;
    std::string _PhotonMVAEstimator;
    
    std::map<int, vector<double> > _phoIDloose;
    std::map<int, vector<double> > _phoIDmedium;
    std::map<int, vector<double> > _phoIDtight;
    std::map<int, vector<double> > _phoISOloose;
    std::map<int, vector<double> > _phoISOmedium;
    std::map<int, vector<double> > _phoISOtight;
    std::map<int, vector<double> > _nhCorr;
    std::map<int, vector<double> > _phCorr;    
    
    vector<std::string> _phoWhichID;
    vector<std::string> _phoWhichISO;
    
    
};

#endif
