#ifndef FLASHgg_bbggJetSystematics_h
#define FLASHgg_bbggJetSystematics_h

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
#include "TString.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"


using namespace std;

class bbggJetSystematics{
public:
	bbggJetSystematics() : isSmearSet(0) {}
	~bbggJetSystematics() {}
	typedef math::XYZTLorentzVector LorentzVector;

	double DeltaR(bbggJetSystematics::LorentzVector vec1, bbggJetSystematics::LorentzVector vec2);
    
    void SetupSmear(std::string resFile, std::string sfFile);
    void SmearJets(std::vector<flashgg::Jet> & Jets, float rho, int variation, std::string randomLabel);
    
    void SetupScale(std::string scaleFile);
    void ScaleJets(std::vector<flashgg::Jet> & Jets, int variation);


private:
    bool isSmearSet;
    JME::JetResolution resolution;
    JME::JetResolutionScaleFactor res_sf;
    JetCorrectionUncertainty *jecUnc;
};

#endif
