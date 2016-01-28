#ifndef FLASHgg_bbggKinFit_h
#define FLASHgg_bbggKinFit_h

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
#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"
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
#include "TMatrixD.h"

using namespace std;

class bbggKinFit{
public:
	bbggKinFit();
	~bbggKinFit() {}
	typedef math::XYZTLorentzVector LorentzVector;
    
    void SetupKinFitter();
	
    void SetJetResolutionParameters(vector<double> etaIntervals, vector<double> pt, vector<double> eta, vector<double> phi);
    float GetPtResolution(bbggKinFit::LorentzVector j);
    float GetEtaResolution(bbggKinFit::LorentzVector j);
    float GetPhiResolution(bbggKinFit::LorentzVector j);
    
    void SetHMass(float hmass) {_hMass = hmass;}
    
    void KinematicFit(bbggKinFit::LorentzVector j1, bbggKinFit::LorentzVector j2);
    
    bbggKinFit::LorentzVector GetJet1() {return Jet1;}
    bbggKinFit::LorentzVector GetJet2() {return Jet2;}
    
	//Other
	double DeltaR(bbggKinFit::LorentzVector vec1, bbggKinFit::LorentzVector vec2);

private:
    float _hMass;
    bbggKinFit::LorentzVector Jet1;
    bbggKinFit::LorentzVector Jet2;
//    TMatrixD m1(3,3);
//    TMatrixD m2(3,3);
    vector<double> ptRes;
    vector<double> etaRes;
    vector<double> phiRes;
    vector<double> etaInts;
            
    bool hasptrel;
    bool hasetarel;
    bool hasphirel;
    bool issetup;
};

#endif

