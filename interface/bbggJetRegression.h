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
	bbggJetRegression();
	~bbggJetRegression() {}
	typedef math::XYZTLorentzVector LorentzVector;
		
	LorentzVector GetRegression(edm::Ptr<flashgg::Jet> jt, std::vector<edm::Ptr<flashgg::Electron>> Electrons, std::vector<edm::Ptr<flashgg::Muon>> Muons, float Rho);
    void SetupRegression(std::string RegType, std::string RegFile);
					 
	//Other
	double DeltaR(bbggJetRegression::LorentzVector vec1, bbggJetRegression::LorentzVector vec2);
    float PtRel(LorentzVector ClosestLep, LorentzVector Jet);
    edm::Ptr<flashgg::Electron> GetClosestElectron(edm::Ptr<flashgg::Jet> jt, std::vector<edm::Ptr<flashgg::Electron>> Electrons);
    edm::Ptr<flashgg::Muon> GetClosestMuon(edm::Ptr<flashgg::Jet> jt, std::vector<edm::Ptr<flashgg::Muon>> Muons);
    float GetLeadingTrackPt(edm::Ptr<flashgg::Jet> jt);

private:
    TMVA::Reader RegressionReader;
    //Regression variables
    float Jet_pt, Jet_rawPt, rho, Jet_eta, Jet_mt, Jet_leadTrackPt, Jet_leptonPtRel, Jet_leptonPt, Jet_leptonDeltaR;
    float Jet_chEmEF, Jet_chHEF, Jet_neHEF, Jet_neEmEF, Jet_chMult, Jet_vtxPt, Jet_vtxMass, Jet_vtx3dL, Jet_vtxNtrk, Jet_vtx3deL;
    bool TMVAReady;
};

#endif

#ifdef bbggJetRegression_cxx
bbggJetRegression::bbggJetRegression(){
    TMVAReady = 0;
}

bbggJetRegression::SetupRegression(std::string RegType, std::string RegFile){
    RegressionReader = TMVA::Reader();
    RegressionReader.AddVariable("Jet_pt", &Jet_pt);
    RegressionReader.AddVariable("Jet_rawPt", &Jet_rawPt);
    RegressionReader.AddVariable("rho", &rho);
    RegressionReader.AddVariable("Jet_eta", &Jet_eta);
    RegressionReader.AddVariable("Jet_mt", &Jet_mt);
    RegressionReader.AddVariable("Jet_leadTrackPt", &Jet_leadTrackPt);
    RegressionReader.AddVariable("Jet_leptonPtRel", &Jet_leptonPtRel);
    RegressionReader.AddVariable("Jet_leptonPt", &Jet_leptonPt);
    RegressionReader.AddVariable("Jet_leptonDeltaR;", &Jet_leptonDeltaR;);
    RegressionReader.AddVariable("Jet_chEmEF", &Jet_chEmEF);
    RegressionReader.AddVariable("Jet_chHEF", &Jet_chHEF);
    RegressionReader.AddVariable("Jet_neHEF", &Jet_neHEF);
    RegressionReader.AddVariable("Jet_neEmEF", &Jet_neEmEF);
    RegressionReader.AddVariable("Jet_chMult", &Jet_chMult);
    RegressionReader.AddVariable("Jet_vtxPt", &Jet_vtxPt);
    RegressionReader.AddVariable("Jet_vtxMass", &Jet_vtxMass);
    RegressionReader.AddVariable("Jet_vtx3dL", &Jet_vtx3dL);
    RegressionReader.AddVariable("Jet_vtxNtrk", &Jet_vtxNtrk);
    RegressionReader.AddVariable("Jet_vtx3deL", &Jet_vtx3deL);
    RegressionReader.BookMVA(RegType, RegFile);
    TMVAReady = 1;
}

#endif // #ifdef bbggJetRegression_cxx
