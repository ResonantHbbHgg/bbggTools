#include "flashgg/bbggTools/interface/bbggJetSystematics.h"
//FLASHgg files
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TVector3.h"
#include "TLorentzVector.h"

#include "JetMETCorrections/Modules/interface/JetResolution.h"


using namespace std;

const int DEBUG = 0;

//bool DEBUG = false;

double bbggJetSystematics::DeltaR(bbggJetSystematics::LorentzVector vec1, bbggJetSystematics::LorentzVector vec2)
{
//	double R2 = (vec1.phi() - vec2.phi())*(vec1.phi() - vec2.phi()) + (vec1.eta() - vec2.eta())*(vec1.eta() - vec2.eta());
//	return sqrt(R2);
	return deltaR(vec1, vec2);
}

void bbggJetSystematics::SetupSmear(std::string resFile, std::string sfFile)
{
    resolution = JME::JetResolution(resFile);
    res_sf = JME::JetResolutionScaleFactor(sfFile);
    isSmearSet = 1;
}

void bbggJetSystematics::SmearJets(std::vector<flashgg::Jet> & Jets, float rho, int variation, std::string randomLabel)
{
    for ( unsigned int jt = 0; jt < Jets.size(); jt++){
        JME::JetParameters jet_param;
        jet_param.setJetPt( Jets[jt].pt() );
        jet_param.setJetEta( Jets[jt].eta() );
        jet_param.setRho(rho);
        
        float res = resolution.getResolution(jet_param);
        float scalefactor = -999;
        if(variation == 1) scalefactor = res_sf.getScaleFactor({{JME::Binning::JetEta, Jets[jt].eta()}}, Variation::UP);
        if(variation == -1) scalefactor = res_sf.getScaleFactor({{JME::Binning::JetEta, Jets[jt].eta()}}, Variation::DOWN);
        
        float recpt = Jets[jt].pt();
        auto genjet = Jets[jt].genJet();
        if(genjet != nullptr) {
            float genpt = genjet->pt();
            float newpt = max<float>(0., genpt + scalefactor*( recpt - genpt ));
            bbggJetSystematics::LorentzVector oldP4 = Jets[jt].p4();
            Jets[jt].setP4( (newpt/recpt)*oldP4 );
        } else {
            if(!Jets[jt].hasUserFloat(randomLabel)) {
                std::cout << "[bbggJetSystematics::SmearJets] Jet doesn't have gen or random label and is not gen matched!" << std::endl;
                continue;
            }
            
            float rnd = Jets[jt].userFloat(randomLabel);
            float extra_smear_width = std::sqrt(scalefactor*scalefactor - 1) * res;
            
            if (extra_smear_width<0) {
                std::cout << "[bbggJetSystematics::SmearJets] Calculated smear is negative!" << std::endl;
                continue;
            }
            
            float escale = (1. + rnd * extra_smear_width);
            bbggJetSystematics::LorentzVector oldP4 = Jets[jt].p4();
            Jets[jt].setP4(escale*oldP4);
        }
    }
}

void bbggJetSystematics::SetupScale(std::string scaleFile)
{ 
    if(DEBUG) std::cout << "[bbggJetSystematics::SetupScale] Setting up scale file " << scaleFile << std::endl;
    jecUnc = new JetCorrectionUncertainty(scaleFile);
    if(DEBUG) std::cout << "[bbggJetSystematics::SetupScale] Scale file set up!" << std::endl;
}

void bbggJetSystematics::ScaleJets(std::vector<flashgg::Jet> & Jets, int variation)
{
    if(DEBUG) std::cout << "[bbggJetSystematics::ScaleJets] Looping over jets..." << std::endl;
    for ( unsigned int jt = 0; jt < Jets.size(); jt++){
        if(DEBUG) std::cout << "[bbggJetSystematics::ScaleJets] Jet precorr pt: "<< Jets[jt].pt() << " eta: " << Jets[jt].eta() << std::endl;
        jecUnc->setJetEta(Jets[jt].eta());
        jecUnc->setJetPt(Jets[jt].pt());
        double unc = jecUnc->getUncertainty(true);
        bbggJetSystematics::LorentzVector corrP4 = Jets[jt].p4()*(1 + variation*unc);
        Jets[jt].setP4( corrP4 );
        if(DEBUG) std::cout << "[bbggJetSystematics::ScaleJets] Jet post corr pt: "<< Jets[jt].pt() << std::endl;
    }
}
