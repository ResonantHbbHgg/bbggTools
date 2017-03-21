#include "flashgg/bbggTools/interface/bbggPhotonCorrector.h"

#include "TVector3.h"
#include "TLorentzVector.h"


using namespace std;

const int DEBUG = 0;

//bool DEBUG = false;


void bbggPhotonCorrector::SetupCorrector(std::string CorrectionFile)
{
    scaler_ = EnergyScaleCorrection_class(CorrectionFile);
    isSet = 1;
}

void bbggPhotonCorrector::SmearPhotonsInDiPhotons(std::vector<flashgg::DiPhotonCandidate> & diPhos, int runNumber)
{
    if(DEBUG) std::cout << "[bbggPhotonCorrector::SmearPhotonsInDiPhotons] Begin" <<std::endl;
    runNumber_ = runNumber;
    scaler_.doSmearings = true;
    scaler_.doScale = false;
    for ( unsigned int dp = 0; dp < diPhos.size(); dp++) {
        float beforeLeadingCorr = diPhos[dp].getLeadingPhoton().energy();
        SmearPhoton(diPhos[dp].getLeadingPhoton());
        SmearPhoton(diPhos[dp].getSubLeadingPhoton());
        float afterLeadingCorr = diPhos[dp].getLeadingPhoton().energy();
        if (DEBUG) std::cout << "[bbggPhotonCorrector::SmearPhotonsInDiPhotons] Before smear: " << beforeLeadingCorr << " - after smear: " << afterLeadingCorr << std::endl;
    }
}

void bbggPhotonCorrector::ScalePhotonsInDiPhotons(std::vector<flashgg::DiPhotonCandidate> & diPhos, int runNumber)
{
    if(DEBUG) std::cout << "[bbggPhotonCorrector::ScalePhotonsInDiPhotons] Begin" <<std::endl;
    runNumber_ = runNumber;
    scaler_.doSmearings = false;
    scaler_.doScale = true;
    for ( unsigned int dp = 0; dp < diPhos.size(); dp++) {
        float beforeLeadingCorr = diPhos[dp].getLeadingPhoton().energy();
        ScalePhoton(diPhos[dp].getLeadingPhoton());
        ScalePhoton(diPhos[dp].getSubLeadingPhoton());
        float afterLeadingCorr = diPhos[dp].getLeadingPhoton().energy();
        if (DEBUG) std::cout << "[bbggPhotonCorrector::ScalePhotonsInDiPhotons] Before scale: " << beforeLeadingCorr << " - after scale: " << afterLeadingCorr << std::endl;
    }
}

void bbggPhotonCorrector::ExtraScalePhotonsInDiPhotons(std::vector<flashgg::DiPhotonCandidate> & diPhos)//, EcalRecHitCollection _ebRecHits)
{
    for ( unsigned int dp = 0; dp < diPhos.size(); dp++) {
        float beforeLeadingCorr = diPhos[dp].getLeadingPhoton().energy();
        ExtraScalePhoton(diPhos[dp].getLeadingPhoton());//, _ebRecHits);
        ExtraScalePhoton(diPhos[dp].getSubLeadingPhoton());//, _ebRecHits);
        float afterLeadingCorr = diPhos[dp].getLeadingPhoton().energy();
        if (DEBUG) std::cout << "[bbggPhotonCorrector::ScalePhotonsInDiPhotons] Before scale: " << beforeLeadingCorr << " - after scale: " << afterLeadingCorr << std::endl;
    }
}

void bbggPhotonCorrector::SmearPhoton(flashgg::Photon & photon)
{
    unsigned int gain=12;
    if(photon.hasSwitchToGain1()) gain=1;
    if(photon.hasSwitchToGain6()) gain=6;
    auto sigma = scaler_.getSmearingSigma(runNumber_, photon.isEB(), photon.full5x5_r9(), photon.superCluster()->eta(), photon.et(), gain, variation_,variation_);
    float rnd = photon.userFloat(randomLabel_);
    float smearing = (1. + rnd * sigma);
    if( DEBUG ) {
        std::cout << "  " << ": Photon has energy= " << photon.energy() << " eta=" << photon.eta()
        << " and we apply a smearing with sigma " << ( 100 * sigma ) << "% to get new energy=" << smearing*photon.energy() << std::endl;
    }
    bbggPhotonCorrector::LorentzVector thisp4 = photon.p4();
    photon.setP4( smearing*thisp4);
}

void bbggPhotonCorrector::ScalePhoton(flashgg::Photon & photon)
{
    unsigned int gain=12;
    if(photon.hasSwitchToGain1()) gain=1;
    if(photon.hasSwitchToGain6()) gain=6;

    auto shift_val = scaler_.ScaleCorrection(runNumber_, photon.isEB(), photon.full5x5_r9(), photon.superCluster()->eta(), photon.et());

    std::bitset<EnergyScaleCorrection_class::scAll> source_all(std::string("111"));
//    std::bitset<EnergyScaleCorrection_class::scAll> stat(std::string("001"));
//    std::bitset<EnergyScaleCorrection_class::scAll> syst(std::string("010"));
//    std::bitset<EnergyScaleCorrection_class::scAll> gain(std::string("100"));

    auto shift_err = scaler_.ScaleCorrectionUncertainty(runNumber_, photon.isEB(), photon.full5x5_r9(), photon.superCluster()->eta(),
							photon.et(), gain, source_all);
    if (abs(variation_)) shift_val = 1.;
    float scale = shift_val + variation_ * shift_err;

    if( DEBUG ) {
        std::cout << "  " << ": Photon has energy= " << photon.energy() << " eta=" << photon.eta()
        << " and we apply a scale factor of " << ( 100 * scale ) << "% to get new energy=" << scale*photon.energy() << std::endl;
    }
    bbggPhotonCorrector::LorentzVector thisp4 = photon.p4();
    photon.setP4( scale*thisp4);
}

void bbggPhotonCorrector::ExtraScalePhoton(flashgg::Photon & photon)//, EcalRecHitCollection _ebrechits)
{
/*
    float seedEnergy = photon.seedEnergy();
    bool isEB = (fabs(photon.superCluster()->eta()) < 1.5) ? 0 : 1;
    double Ecorr = 1;
    if (!isEB) Ecorr = 1;
    else {
        if(seedEnergy > 200 && seedEnergy < 300)  Ecorr=1.0199 + variation_*0.0008;
        else if(seedEnergy > 300 && seedEnergy < 400) Ecorr=  1.052 + variation_*0.003;
        else if(seedEnergy > 400 && seedEnergy < 500) Ecorr = 1.015 + variation_*0.006;
    }
*/
/* 
    DetId detid = photon.superCluster()->seed()->seed();
    const EcalRecHit * rh = NULL;
    if (detid.subdetId() == EcalBarrel) {
        auto rh_i =  _ebrechits.find(detid);
        if( rh_i != _ebrechits.end()) rh =  &(*rh_i);
        else rh = NULL;
    } 

    if(rh==NULL) Ecorr=1;
    else{
        if(rh->energy() > 200 && rh->energy()<300)  Ecorr=1.0199 + variation_*0.0008;
        else if(rh->energy()>300 && rh->energy()<400) Ecorr=  1.052 + variation_*0.003;
        else if(rh->energy()>400 && rh->energy()<500) Ecorr = 1.015 + variation_*0.006;
    }
*/
    double Ecorr = 1; //This is not used anymore!
    bbggPhotonCorrector::LorentzVector thisp4 = photon.p4();
    photon.setP4( Ecorr*thisp4);
}

void bbggPhotonCorrector::SetCustomPhotonIDMVA( std::vector<flashgg::DiPhotonCandidate> & diPhos, edm::Handle<edm::ValueMap<float> > mvaValues)
{
    for (unsigned int iy = 0; iy < diPhos.size(); iy++)
    {
        float leadingMVA = (*mvaValues)[diPhos[iy].leadingView()->originalPhoton()];
        diPhos[iy].getLeadingPhoton().addUserFloat("EGMMVAID", leadingMVA);
        float subleadingMVA = (*mvaValues)[diPhos[iy].subLeadingView()->originalPhoton()];
        diPhos[iy].getSubLeadingPhoton().addUserFloat("EGMMVAID", subleadingMVA);
        if(DEBUG) {
            std::cout << "[SetCustomPhotonIDMVA] Leading photon EGM MVA: " << diPhos[iy].leadingPhoton()->userFloat("EGMMVAID") << std::endl;
            std::cout << "[SetCustomPhotonIDMVA] SubLeading photon EGM MVA: " << diPhos[iy].subLeadingPhoton()->userFloat("EGMMVAID") << std::endl;
        }
    }
}
