#ifndef FLASHgg_bbggPhotonCorrector_h
#define FLASHgg_bbggPhotonCorrector_h

//FLASHgg files
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
// Other DataFormats
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

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
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrection_class.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"


using namespace std;

class bbggPhotonCorrector{
public:
    bbggPhotonCorrector() : isSet(0), variation_(0) {}
    ~bbggPhotonCorrector() {}
    typedef math::XYZTLorentzVector LorentzVector;
 
    void setRandomLabel(std::string randomLabel) { randomLabel_ = randomLabel;};
    void setVariation(float variation) {variation_ = variation;};

    void SetupCorrector(std::string CorrectionFile);

    void SmearPhotonsInDiPhotons(std::vector<flashgg::DiPhotonCandidate> & diPhos, int runNumber);
    void SmearPhoton(flashgg::Photon & photon);
    
    void ScalePhotonsInDiPhotons(std::vector<flashgg::DiPhotonCandidate> & diPhos, int runNumber);
    void ScalePhoton(flashgg::Photon & photon);

    void ExtraScalePhotonsInDiPhotons(std::vector<flashgg::DiPhotonCandidate> & diPhos);//, EcalRecHitCollection _ebRecHits);
    void ExtraScalePhoton(flashgg::Photon & photon);//, EcalRecHitCollection _ebRecHits);

    void SetCustomPhotonIDMVA(std::vector<flashgg::DiPhotonCandidate> & diPhos, edm::Handle<edm::ValueMap<float> > mvaValues);

private:
    bool isSet;
    float variation_;
    int runNumber_;
    std::string randomLabel_;
    EnergyScaleCorrection_class scaler_;
};

#endif
