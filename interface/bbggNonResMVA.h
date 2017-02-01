#ifndef FLASHgg_bbggNonResMVA_h
#define FLASHgg_bbggNonResMVA_h

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

class bbggNonResMVA{
public:
	bbggNonResMVA() : TMVAReady(0) {}
	void SetupNonResMVA(std::string LMFile, std::string HMFile, std::vector<std::string> inVars);
	~bbggNonResMVA() {}
	typedef math::XYZTLorentzVector LorentzVector;

        std::vector<float> mvaDiscriminants(std::map<std::string,float> params);

	//Other
	double DeltaR(bbggNonResMVA::LorentzVector vec1, bbggNonResMVA::LorentzVector vec2);

private:
    	TMVA::Reader *LowMassReader;
    	TMVA::Reader *HighMassReader;
	bool TMVAReady;
    //MVA variables
        std::map<std::string, float> mvaVars;
        std::vector<std::string> orderedVars;

};

#endif
