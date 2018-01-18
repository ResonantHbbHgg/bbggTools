#ifndef FLASHgg_bbggNonResMVA2017_h
#define FLASHgg_bbggNonResMVA2017_h

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

class bbggNonResMVA2017{
public:
	bbggNonResMVA2017() : TMVAReady(0) {}
	void SetupNonResMVA(std::string MVAFile, std::vector<std::string> inVars);
	~bbggNonResMVA2017() {}
	typedef math::XYZTLorentzVector LorentzVector;

        std::vector<float> mvaDiscriminants(std::map<std::string,float> params);

	//Other
	double DeltaR(bbggNonResMVA2017::LorentzVector vec1, bbggNonResMVA2017::LorentzVector vec2);

private:
    	TMVA::Reader *MVAReader;
	bool TMVAReady;
    //MVA variables
        std::map<std::string, float> mvaVars;
        std::vector<std::string> orderedVars;

};

#endif
