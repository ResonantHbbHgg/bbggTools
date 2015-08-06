#ifndef FLASHgg_bbggMC_h
#define FLASHgg_bbggMC_h

//FLASHgg files
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
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

class bbggMC{
public:
	bbggMC() : hasDiPho( 0 ), hasLeadJet( 0 ), hasSubJet( 0 ) { }
	~bbggMC() {}
	typedef math::XYZTLorentzVector LorentzVector;
	
	//Set cuts for selection
		
	//Perform event selection
	bool MatchTruth( edm::Handle<edm::View<flashgg::DiPhotonCandidate> > diphoCol, 
							edm::Handle<edm::View<flashgg::Jet> > jetsCol,
							edm::Handle<edm::View<reco::GenParticle> > genCol );
	//Get selected objects
	edm::Ptr<flashgg::DiPhotonCandidate> GetSelected_diphoCandidate();
	edm::Ptr<flashgg::Jet> GetSelected_leadingJetCandidate();
	edm::Ptr<flashgg::Jet> GetSelected_subleadingJetCandidate();
	
	//Other
	double DeltaR(bbggMC::LorentzVector vec1, bbggMC::LorentzVector vec2);

private:
	double rho_;
	
	//Selected objects
	edm::Ptr<flashgg::DiPhotonCandidate> diphoCandidate;
	bool hasDiPho;
	edm::Ptr<flashgg::Jet> leadingJetCandidate;
	bool hasLeadJet;
	edm::Ptr<flashgg::Jet> subleadingJetCandidate;
	bool hasSubJet;
};

#endif
