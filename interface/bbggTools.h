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
	bbggTools() { rho_ = -10;}
	~bbggTools() {}
	typedef math::XYZTLorentzVector LorentzVector;
	double getCHisoToCutValue(edm::Ptr<flashgg::DiPhotonCandidate> dipho, int whichPho);
    double getNHisoToCutValue(const flashgg::Photon* pho);
	double getPHisoToCutValue(const flashgg::Photon* pho);
	double getEA( float eta, int whichEA);
	double DeltaR( bbggTools::LorentzVector vec1, bbggTools::LorentzVector vec2);
	bool isPhoID(edm::Ptr<flashgg::Photon> pho, vector<double> cuts);
	bool isPhoISO(edm::Ptr<flashgg::DiPhotonCandidate> pho, int whichPho, vector<double> cuts);
	void setRho(double rho) {rho_ = rho;}
private:
	double rho_;
};

#endif
