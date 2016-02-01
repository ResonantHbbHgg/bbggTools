#include "flashgg/bbggTools/interface/bbggJetRegression.h"
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
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

const int DEBUG = 0;

//bool DEBUG = false;

double bbggJetRegression::DeltaR(bbggJetRegression::LorentzVector vec1, bbggJetRegression::LorentzVector vec2)
{
	double R2 = (vec1.phi() - vec2.phi())*(vec1.phi() - vec2.phi()) + (vec1.eta() - vec2.eta())*(vec1.eta() - vec2.eta());
	return sqrt(R2);
}

bbggJetRegression::LorentzVector bbggJetRegression::GetRegression(edm::Ptr<flashgg::Jet> jt, std::vector<edm::Ptr<flashgg::Electron>> Electrons, std::vector<edm::Ptr<flashgg::Muon>> Muons, float Rho)
{
    bbggJetRegression::LorentzVector NewJet(-1,-1,-1,-1);
    
    if(TMVAReady == 0){
        std::cout << "[bbggJetRegression::GetRegression] You haven't yet setup the regression reader! First do bbggJetRegression::SetupRegression(std::string RegType, std::string RegFile)" << std::endl;
        return NewJet;
    }
    
    //Assigning regression variables:
    Jet_pt = jt->pt();
    Jet_rawPt = jt->correctedJet("Uncorrected").pt();
    rho = Rho;
    Jet_eta = jt->eta();
    Jet_mt = jt->mt();
    
    Jet_leadTrackPt = bbggJetRegression::GetLeadingTrackPt(jt);
        
    float softLeptonPt = jt->softLeptonPt();
    float softLeptonDR = jt->softLeptonDeltaR();
    
    Jet_leptonPtRel = softLeptonPt/Jet_pt; //bbggJetRegression::PtRel(ClosestLepton, jt->p4());
    Jet_leptonPt = softLeptonPt; //ClosestLepton.pt();
    Jet_leptonDeltaR = softLeptonDR; //bbggJetRegression::DeltaR(ClosestLepton, jt->p4());
        
    Jet_chEmEF = jt->chargedEmEnergyFraction();
    Jet_chHEF = jt->chargedHadronEnergyFraction();
    Jet_neHEF = jt->neutralEmEnergyFraction();
    Jet_neEmEF = jt->neutralHadronEnergyFraction();
    Jet_chMult = jt->chargedHadronMultiplicity();
    Jet_vtxPt = sqrt( jt->userFloat("vtxPx")*jt->userFloat("vtxPx") + jt->userFloat("vtxPy")*jt->userFloat("vtxPy"));
    Jet_vtxMass = jt->userFloat("vtxMass");
    Jet_vtx3dL = jt->userFloat("vtx3DVal");
    Jet_vtxNtrk = jt->userFloat("vtxNtracks");
    Jet_vtx3deL = jt->userFloat("vtx3DSig");
    
    const std::vector< Float_t > & RegressedValues = RegressionReader.EvaluateRegression("pt_reg");
    float newPt = RegressedValues[0];
    float newEta = RegressedValues[3];
    
    TLorentzVector TempJet;
    TempJet.SetPtEtaPhiM(newPt, newEta, jt->phi(), jt->mass());
    NewJet = bbggJetRegression::LorentzVector(TempJet.Px(), TempJet.Py(), TempJet.Pz(), TempJet.E());
    return NewJet;
}

float bbggJetRegression::PtRel(LorentzVector ClosestLep, LorentzVector Jet){
    TVector3 A(Jet.Vect().X(), Jet.Vect().Y(), Jet.Vect().Z());
    TLorentzVector O(ClosestLep.Px(),ClosestLep.Py(),ClosestLep.Pz(),ClosestLep.E());
    return O.Perp(A);
}

float bbggJetRegression::GetLeadingTrackPt(edm::Ptr<flashgg::Jet> jt)
{
    reco::TrackRefVector JetTracks = jt->associatedTracks();
    unsigned int nTracks = jt->associatedTracks().size();
    float MaxPt = 0;
    for(unsigned int tr = 0; tr < nTracks; tr++){
        float TrackPt = jt->associatedTracks()[tr]->pt();
        if(TrackPt > MaxPt)
            MaxPt = TrackPt;
    }
    return MaxPt;
}