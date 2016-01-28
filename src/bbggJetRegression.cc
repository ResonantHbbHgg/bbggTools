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
    Jet_rawPt = jt->correctedJet("raw").pt();
    rho = Rho;
    Jet_eta = jt->eta();
    Jet_mt = jt->mt();
    
    Jet_leadTrackPt = bbggJetRegression::GetLeadingTrackPt(jt);
        
    std::vector<LorentzVector> ClosestLeptonVec;
    
    if(Electrons.size() > 0){
        edm::Ptr<flashgg::Electron> closestEl = bbggJetRegression::GetClosestElectron(jt, Electrons);
        ClosestLeptonVec.push_back(closestEl->p4());
    }
    
    if(Muons.size() > 0){
        edm::Ptr<flashgg::Muon> closestMu = bbggJetRegression::GetClosestMuon(jt, Muons);
        ClosestLeptonVec.push_back(closestMu->p4());
    }

    LorentzVector ClosestLepton;
    if(ClosestLeptonVec.size() == 0){
        Jet_leptonPtRel = -99;
        Jet_leptonPt = -99;
        Jet_leptonDeltaR = -99;
    }
    if(ClosestLeptonVec.size() == 1)
        ClosestLepton = ClosestLeptonVec[0];
    if(ClosestLeptonVec.size() == 2){
        float DR0 = bbggJetRegression::DeltaR(ClosestLeptonVec[0], jt->p4());
        float DR1 = bbggJetRegression::DeltaR(ClosestLeptonVec[1], jt->p4());
        ClosestLepton = (DR0 < DR1) ? (ClosestLeptonVec[0]) : (ClosestLeptonVec[1]);
    }
    
    Jet_leptonPtRel = bbggJetRegression::PtRel(ClosestLepton, jt->p4());
    Jet_leptonPt = ClosestLepton.pt();
    Jet_leptonDeltaR = bbggJetRegression::DeltaR(ClosestLepton, jt->p4());
        
    Jet_chEmEF = jt->chargedEmEnergyFraction();
    Jet_chHEF = jt->chargedHadronEnergyFraction();
    Jet_neHEF = jt->neutralEmEnergyFraction();
    Jet_neEmEF = jt->neutralHadronEnergyFraction();
    Jet_chMult = jt->chargedHadronMultiplicity();
    Jet_vtxPt = sqrt( jt->userFloat("vtxPx")*jt->userFloat("vtxPx") + jt->userFloat("vtxPy")*jt->userFloat("vtxPy"));
    Jet_vtxMass = jt->userFloat("vtxMass");
    Jet_vtx3dL = jt->userFloat("vtx3dL");
    Jet_vtxNtrk = jt->userFloat("vtxNtrk");
    Jet_vtx3deL = jt->userFloat("vtx3deL");
    
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

edm::Ptr<flashgg::Electron> bbggJetRegression::GetClosestElectron(edm::Ptr<flashgg::Jet> jt, std::vector<edm::Ptr<flashgg::Electron>> Electrons)
{
    std::vector<edm::Ptr<flashgg::Electron>> closestEl;
    float MinDr = 100;
    for(unsigned int el = 0; el < Electrons.size(); el++){
        edm::Ptr<flashgg::Electron> thisEl = Electrons[el];
        if(bbggJetRegression::DeltaR(thisEl->p4(), jt->p4()) < MinDr ){
            MinDr = bbggJetRegression::DeltaR(thisEl->p4(), jt->p4());
            closestEl.clear();
            closestEl.push_back(thisEl);
        }
    }
    return closestEl[0]; 
}

edm::Ptr<flashgg::Muon> bbggJetRegression::GetClosestMuon(edm::Ptr<flashgg::Jet> jt, std::vector<edm::Ptr<flashgg::Muon>> Muons)
{
    std::vector<edm::Ptr<flashgg::Muon>> closestMu;
    float MinDr = 100;
    for(unsigned int mu = 0; mu < Muons.size(); mu++){
        edm::Ptr<flashgg::Muon> thisMu = Muons[mu];
        if(bbggJetRegression::DeltaR(thisMu->p4(), jt->p4()) < MinDr ){
            MinDr = bbggJetRegression::DeltaR(thisMu->p4(), jt->p4());
            closestMu.clear();
            closestMu.push_back(thisMu);
        }
    }
    return closestMu[0]; 
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