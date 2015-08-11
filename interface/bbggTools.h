#ifndef FLASHgg_bbggTools_h
#define FLASHgg_bbggTools_h

//FLASHgg files
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
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
	bbggTools() : hasDiPho( 0 ), hasLeadJet( 0 ), hasSubJet( 0 ) { rho_ = -10;}
	~bbggTools() {}
	typedef math::XYZTLorentzVector LorentzVector;
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
	bool isPhoISO(edm::Ptr<flashgg::DiPhotonCandidate> pho, int whichPho, vector<double> cuts);
	void setRho(double rho) {rho_ = rho;}
	
	//Set cuts for selection
	void SetCut_PhotonPtOverDiPhotonMass( vector<double> cuts) { _PhotonPtOverDiPhotonMass = cuts; }
	void SetCut_PhotonEta( vector<double> cuts) { _PhotonEta = cuts; }
	void SetCut_PhotonDoID( vector<int> cuts) { _PhotonDoID = cuts; }
	void SetCut_PhotonDoISO( vector<int> cuts) { _PhotonDoISO = cuts; }
	void SetCut_PhotonHoverE( vector<double> cuts) { _PhotonHoverE = cuts; }
	void SetCut_PhotonSieie( vector<double> cuts) { _PhotonSieie = cuts; }
	void SetCut_PhotonR9( vector<double> cuts) { _PhotonR9 = cuts; }
	void SetCut_PhotonChargedIso( vector<double> cuts) { _PhotonChargedIso = cuts; }
	void SetCut_PhotonNeutralIso( vector<double> cuts) { _PhotonNeutralIso = cuts; }
	void SetCut_PhotonPhotonIso( vector<double> cuts) { _PhotonPhotonIso = cuts; }
	void SetCut_PhotonElectronVeto( vector<int> cuts) { _PhotonElectronVeto = cuts; }
	void SetCut_DiPhotonPt( vector<double> cuts) { _DiPhotonPt = cuts; }
	void SetCut_DiPhotonEta( vector<double> cuts) { _DiPhotonEta = cuts; }
	void SetCut_DiPhotonMassWindow( vector<double> cuts) { _DiPhotonMassWindow = cuts; }
	void SetCut_DiPhotonOnlyFirst( unsigned int cuts) { _DiPhotonOnlyFirst = cuts; }
	void SetCut_JetPt( vector<double> cuts) { _JetPt = cuts; }
	void SetCut_JetEta( vector<double> cuts) { _JetEta = cuts; }
	void SetCut_JetBDiscriminant( vector<double> cuts) { _JetBDiscriminant = cuts; }
	void SetCut_JetDrPho( vector<double> cuts ) { _JetDrPho = cuts; }
	void SetCut_JetDoPUID( vector<int> cuts) { _JetDoPUID = cuts; }
	void SetCut_n_bJets( unsigned int cuts) { _n_bJets = cuts; }
	void SetCut_DiJetPt( vector<double> cuts) { _DiJetPt = cuts; }
	void SetCut_DiJetEta( vector<double> cuts) { _DiJetEta = cuts; }
	void SetCut_DiJetMassWindow( vector<double> cuts) { _DiJetMassWindow = cuts; }
	void SetCut_CandidateMassWindow( vector<double> cuts) { _CandidateMassWindow = cuts; }
	void SetCut_CandidatePt( vector<double> cuts) { _CandidatePt = cuts; }
	void SetCut_CandidateEta( vector<double> cuts) { _CandidateEta = cuts; }
	void SetCut_bTagType( std::string cuts) { _bTagType = cuts; }
	void SetCut_CandidatesDeltaR( vector<double> cuts) { _CandidatesDeltaR = cuts; }
		
	//Perform event selection
	bool AnalysisSelection( edm::Handle<edm::View<flashgg::DiPhotonCandidate> > diphoCol, 
							edm::Handle<edm::View<flashgg::Jet> > jetsCol );
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
	
	//Cut values
	vector<double> _PhotonPtOverDiPhotonMass;
	vector<double> _PhotonEta;
	vector<int> _PhotonDoID;
	vector<int> _PhotonDoISO;
	vector<double> _PhotonHoverE;
	vector<double> _PhotonSieie;
	vector<double> _PhotonR9;
	vector<double> _PhotonChargedIso;
	vector<double> _PhotonNeutralIso;
	vector<double> _PhotonPhotonIso;
	vector<int> _PhotonElectronVeto;
	vector<double> _DiPhotonPt;
	vector<double> _DiPhotonEta;
	vector<double> _DiPhotonMassWindow;
	unsigned int _DiPhotonOnlyFirst;
	vector<double> _JetPt;
	vector<double> _JetEta;
	vector<double> _JetBDiscriminant;
	vector<double> _JetDrPho;
	vector<int> _JetDoPUID;
	unsigned int _n_bJets;
	vector<double> _DiJetPt;
	vector<double> _DiJetEta;
	vector<double> _DiJetMassWindow;
	vector<double> _CandidateMassWindow;
	vector<double> _CandidatePt;
	vector<double> _CandidateEta;
	std::string _bTagType;
	vector<double> _CandidatesDeltaR;
};

#endif
