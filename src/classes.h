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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//FLASHgg files
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"

//Local
#include "flashgg/bbggTools/interface/bbggJetRegression.h"
#include "flashgg/bbggTools/interface/bbggJetSystematics.h"
#include "flashgg/bbggTools/interface/bbggKinFit.h"
#include "flashgg/bbggTools/interface/bbggMC.h"
#include "flashgg/bbggTools/interface/bbggPhotonCorrector.h"
#include "flashgg/bbggTools/interface/bbggTools.h"

namespace {
	struct dictionary {
           std::pair<flashgg::Jet, flashgg::Jet> dummy_jetpair;
           std::vector<std::pair<flashgg::Jet, flashgg::Jet>> dummy_jetpair_vec;
//		bbggTools dummy_tools;
//       bbggLTMaker dummy_ltree;
//        bbggMC dummy_mc;
//        bbgg2DFitter dummy_fit;
	};
}
