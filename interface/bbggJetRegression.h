#ifndef FLASHgg_bbggJetRegression_h
#define FLASHgg_bbggJetRegression_h

//FLASHgg files
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
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
#include "TMVA/Reader.h"

using namespace std;

class bbggJetRegression{
public:
	bbggJetRegression() : TMVAReady(0) {}
//	bbggJetRegression(std::string RegType, std::string RegFile);
	void SetupRegression(std::string RegType, std::string RegFile);
	~bbggJetRegression() {}
	typedef math::XYZTLorentzVector LorentzVector;
		
	LorentzVector GetRegression(edm::Ptr<flashgg::Jet> jt, float Rho);
					 
	//Other
	double DeltaR(bbggJetRegression::LorentzVector vec1, bbggJetRegression::LorentzVector vec2);
    	float PtRel(LorentzVector ClosestLep, LorentzVector Jet);
    	edm::Ptr<flashgg::Electron> GetClosestElectron(edm::Ptr<flashgg::Jet> jt, std::vector<edm::Ptr<flashgg::Electron>> Electrons);
    	edm::Ptr<flashgg::Muon> GetClosestMuon(edm::Ptr<flashgg::Jet> jt, std::vector<edm::Ptr<flashgg::Muon>> Muons);
    	float GetLeadingTrackPt(edm::Ptr<flashgg::Jet> jt);

private:
    	TMVA::Reader *RegressionReader;
	bool TMVAReady;
    //Regression variables
    	float Jet_pt, Jet_corr, rho, Jet_eta, Jet_mt, Jet_leadTrackPt, Jet_leptonPtRel, Jet_leptonPt, Jet_leptonDeltaR;
    	float /*Jet_chEmEF, Jet_chHEF,*/ Jet_neHEF, Jet_neEmEF, Jet_chMult, Jet_vtxPt, Jet_vtxMass, Jet_vtx3dL, Jet_vtxNtrk, Jet_vtx3deL;
};

#endif
