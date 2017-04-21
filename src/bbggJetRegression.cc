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
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

const int DEBUG = 0;

//bool DEBUG = false;


void bbggJetRegression::SetupRegression(std::string RegType, std::string RegFile){
    TMVAReady = 1;
    RegressionReader = new TMVA::Reader();
    RegressionReader->AddVariable("Jet_pt", &Jet_pt);
    RegressionReader->AddVariable("Jet_corr", &Jet_corr);
    RegressionReader->AddVariable("Jet_eta", &Jet_eta);
    RegressionReader->AddVariable("Jet_mt", &Jet_mt);
    RegressionReader->AddVariable("Jet_leadTrackPt", &Jet_leadTrackPt);
    RegressionReader->AddVariable("Jet_leptonPtRel", &Jet_leptonPtRel);
    RegressionReader->AddVariable("Jet_leptonPt", &Jet_leptonPt);
    RegressionReader->AddVariable("Jet_leptonDeltaR", &Jet_leptonDeltaR);
//    RegressionReader->AddVariable("Jet_chHEF+Jet_neHEF", &Jet_chHEF_neHEF);
    RegressionReader->AddVariable("Jet_totHEF", &Jet_totHEF);
    RegressionReader->AddVariable("Jet_neEmEF", &Jet_neEmEF);
//    RegressionReader->AddVariable("Jet_chMult", &Jet_chMult);
//    RegressionReader->AddVariable("Jet_chEmEF", &Jet_chEmEF);
//    RegressionReader->AddVariable("Jet_chHEF", &Jet_chHEF);
    RegressionReader->AddVariable("Jet_vtxPt", &Jet_vtxPt);
    RegressionReader->AddVariable("Jet_vtxMass", &Jet_vtxMass);
    RegressionReader->AddVariable("Jet_vtx3dL", &Jet_vtx3dL);
    RegressionReader->AddVariable("Jet_vtxNtrk", &Jet_vtxNtrk);
    RegressionReader->AddVariable("Jet_vtx3deL", &Jet_vtx3deL);
    RegressionReader->AddVariable("nPVs", &nPVs);
    RegressionReader->AddVariable("Jet_PFMET", &Jet_PFMET);
    RegressionReader->AddVariable("Jet_METDPhi", &Jet_METDPhi);
//    RegressionReader->AddVariable("Jet_JetDR", &Jet_JetDR);
    RegressionReader->BookMVA(RegType, RegFile);
}

void bbggJetRegression::SetupRegression(std::string RegType, std::string RegFileLeading, std::string RegFileSubLeading){
    TMVAReady = 1;
    RegressionReaderLeading = new TMVA::Reader();
    RegressionReaderLeading->AddVariable("Jet_pt", &Jet_pt);
//    RegressionReaderLeading->AddVariable("Jet_corr", &Jet_corr);
    RegressionReaderLeading->AddVariable("Jet_eta", &Jet_eta);
    RegressionReaderLeading->AddVariable("Jet_mt", &Jet_mt);
    RegressionReaderLeading->AddVariable("Jet_leadTrackPt", &Jet_leadTrackPt);
    RegressionReaderLeading->AddVariable("Jet_leptonPtRel", &Jet_leptonPtRel);
    RegressionReaderLeading->AddVariable("Jet_leptonPt", &Jet_leptonPt);
    RegressionReaderLeading->AddVariable("Jet_leptonDeltaR", &Jet_leptonDeltaR);
//    RegressionReaderLeading->AddVariable("Jet_totHEF", &Jet_totHEF);
    RegressionReaderLeading->AddVariable("Jet_neHEF", &Jet_neHEF);
    RegressionReaderLeading->AddVariable("Jet_neEmEF", &Jet_neEmEF);
    RegressionReaderLeading->AddVariable("Jet_vtxPt", &Jet_vtxPt);
    RegressionReaderLeading->AddVariable("Jet_vtxMass", &Jet_vtxMass);
    RegressionReaderLeading->AddVariable("Jet_vtx3dL", &Jet_vtx3dL);
    RegressionReaderLeading->AddVariable("Jet_vtxNtrk", &Jet_vtxNtrk);
    RegressionReaderLeading->AddVariable("Jet_vtx3deL", &Jet_vtx3deL);
    RegressionReaderLeading->AddVariable("nGoodPVs", &nPVs);
    RegressionReaderLeading->AddVariable("Jet_PFMET", &Jet_PFMET);
    RegressionReaderLeading->AddVariable("Jet_METDPhi", &Jet_METDPhi);
    RegressionReaderLeading->AddVariable("Jet_JetDR", &Jet_JetDR);
    RegressionReaderLeading->BookMVA(RegType, RegFileLeading);

    RegressionReaderSubLeading = new TMVA::Reader();
    RegressionReaderSubLeading->AddVariable("Jet_pt", &Jet_pt);
//    RegressionReaderSubLeading->AddVariable("Jet_corr", &Jet_corr);
    RegressionReaderSubLeading->AddVariable("Jet_eta", &Jet_eta);
    RegressionReaderSubLeading->AddVariable("Jet_mt", &Jet_mt);
    RegressionReaderSubLeading->AddVariable("Jet_leadTrackPt", &Jet_leadTrackPt);
    RegressionReaderSubLeading->AddVariable("Jet_leptonPtRel", &Jet_leptonPtRel);
    RegressionReaderSubLeading->AddVariable("Jet_leptonPt", &Jet_leptonPt);
    RegressionReaderSubLeading->AddVariable("Jet_leptonDeltaR", &Jet_leptonDeltaR);
//    RegressionReaderSubLeading->AddVariable("Jet_totHEF", &Jet_totHEF);
    RegressionReaderSubLeading->AddVariable("Jet_neHEF", &Jet_neHEF);
    RegressionReaderSubLeading->AddVariable("Jet_neEmEF", &Jet_neEmEF);
    RegressionReaderSubLeading->AddVariable("Jet_vtxPt", &Jet_vtxPt);
    RegressionReaderSubLeading->AddVariable("Jet_vtxMass", &Jet_vtxMass);
    RegressionReaderSubLeading->AddVariable("Jet_vtx3dL", &Jet_vtx3dL);
    RegressionReaderSubLeading->AddVariable("Jet_vtxNtrk", &Jet_vtxNtrk);
    RegressionReaderSubLeading->AddVariable("Jet_vtx3deL", &Jet_vtx3deL);
    RegressionReaderSubLeading->AddVariable("nGoodPVs", &nPVs);
    RegressionReaderSubLeading->AddVariable("Jet_PFMET", &Jet_PFMET);
    RegressionReaderSubLeading->AddVariable("Jet_METDPhi", &Jet_METDPhi);
    RegressionReaderSubLeading->AddVariable("Jet_JetDR", &Jet_JetDR);
    RegressionReaderSubLeading->BookMVA(RegType, RegFileSubLeading);
}

double bbggJetRegression::DeltaR(bbggJetRegression::LorentzVector vec1, bbggJetRegression::LorentzVector vec2)
{
//	double R2 = (vec1.phi() - vec2.phi())*(vec1.phi() - vec2.phi()) + (vec1.eta() - vec2.eta())*(vec1.eta() - vec2.eta());
//	return sqrt(R2);
	return deltaR(vec1, vec2);
}

void bbggJetRegression::RegressedJets(std::vector<flashgg::Jet> & jets, int npvs, bbggJetRegression::LorentzVector theMET)
{
   if(TMVAReady == 0){
       std::cout << "[bbggJetRegression::GetRegression] You haven't yet setup the regression reader! First do bbggJetRegression::SetupRegression(std::string RegType, std::string RegFile)" << std::endl;
        return;
   }

   for ( unsigned int jt = 0; jt < jets.size(); jt++) { bbggJetRegression::GetRegressedJet( jets[jt], npvs, theMET); }

}

void bbggJetRegression::RegressedJets(std::vector<std::pair<flashgg::Jet, flashgg::Jet>> & jets, int npvs, bbggJetRegression::LorentzVector theMET)
{
   if(TMVAReady == 0){
       std::cout << "[bbggJetRegression::GetRegression] You haven't yet setup the regression reader! First do bbggJetRegression::SetupRegression(std::string RegType, std::string RegFile)" << std::endl;
        return;
   }

   for ( unsigned int jt = 0; jt < jets.size(); jt++) { 
      bbggJetRegression::GetRegressedJet( jets[jt].first, npvs, theMET, jets[jt].second, 1); 
      bbggJetRegression::GetRegressedJet( jets[jt].second, npvs, theMET, jets[jt].first, 0); 
   }

}

void bbggJetRegression::GetRegressedJet(flashgg::Jet & jt, int npvs, bbggJetRegression::LorentzVector theMET, flashgg::Jet oJet, bool isLeading)
{

   if(TMVAReady == 0){
       std::cout << "[bbggJetRegression::GetRegression] You haven't yet setup the regression reader! First do bbggJetRegression::SetupRegression(std::string RegType, std::string RegFile)" << std::endl;
        return;
   }

   Jet_pt = jt.pt();
//   Jet_corr = jt.jecFactor("Uncorrected");
   float JetRawpt = jt.correctedJet("Uncorrected").pt();
   Jet_corr = Jet_pt/JetRawpt;
   if(DEBUG) std::cout << "JETCORR:: " << jt.jecFactor("Uncorrected") << " \t " << Jet_corr << std::endl;
   nPVs = npvs;
   Jet_eta = jt.eta();
   Jet_mt = jt.mt();
   Jet_leadTrackPt = jt.userFloat("leadTrackPt");
   if(Jet_leadTrackPt==-999) Jet_leadTrackPt=0;
   Jet_leptonPtRel = jt.userFloat("softLepRatio");
   if(Jet_leptonPtRel==-999) Jet_leptonPtRel=0;
   Jet_leptonPt = jt.userFloat("softLepPt");
   if(Jet_leptonPt==-999) Jet_leptonPt=0;
   Jet_leptonDeltaR = jt.userFloat("softLepDr");
   if(Jet_leptonDeltaR==-999) Jet_leptonDeltaR=0;
   Jet_totHEF = jt.chargedHadronEnergyFraction() + jt.neutralHadronEnergyFraction();
   Jet_neHEF = jt.neutralHadronEnergyFraction();
   Jet_neEmEF = jt.neutralEmEnergyFraction();
   Jet_vtxPt = sqrt( jt.userFloat("vtxPx")*jt.userFloat("vtxPx") + jt.userFloat("vtxPy")*jt.userFloat("vtxPy"));
   if( jt.userFloat("vtxPx") == -999 || jt.userFloat("vtxPy") == -999) Jet_vtxPt = 0;
   Jet_vtxMass = jt.userFloat("vtxMass");
   if(Jet_vtxMass==-999) Jet_vtxMass = 0;
   Jet_vtx3dL = jt.userFloat("vtx3DVal");
   if(Jet_vtx3dL==-999) Jet_vtx3dL = 0;
   Jet_vtxNtrk = jt.userFloat("vtxNTracks");
   if(Jet_vtxNtrk==-999) Jet_vtxNtrk = 0;
   if(Jet_vtx3deL==-999 || Jet_vtx3dL == 0) Jet_vtx3deL = 0;
   else Jet_vtx3deL = Jet_vtx3dL/jt.userFloat("vtx3DSig");
//   if(Jet_vtx3deL==-999) Jet_vtx3deL = 0;
   //MET stuff
   Jet_PFMET = theMET.pt();
   Jet_METDPhi = deltaPhi( jt.phi(), theMET.phi() );
   //DeltaR
   Jet_JetDR = bbggJetRegression::DeltaR(jt.p4(), oJet.p4());
   

   const std::vector< Float_t > & RegressedValuesL = RegressionReaderLeading->EvaluateRegression("BDTG method");
   const std::vector< Float_t > & RegressedValuesS = RegressionReaderSubLeading->EvaluateRegression("BDTG method");

   float corrFactor = 1.;
   if(isLeading)
     corrFactor = RegressedValuesL[0];
   else
     corrFactor = RegressedValuesS[0];
     
   float oldEnergy = jt.energy();
   float oldPhi = jt.phi();
   TLorentzVector tempVec;
   tempVec.SetPtEtaPhiE(corrFactor*Jet_pt, Jet_eta, oldPhi, corrFactor*oldEnergy);
   bbggJetRegression::LorentzVector temp2;
   temp2.SetPxPyPzE( tempVec.Px(), tempVec.Py(), tempVec.Pz(), tempVec.E());
   jt.setP4(temp2);

   if(DEBUG) {
      std::cout << "[bbggJetRegression::GetRegressedJet] Old pt : " << Jet_pt << " - pt from reg : " << corrFactor*Jet_pt << " - new pt : " << jt.pt() << std::endl;
   }
}

void bbggJetRegression::GetRegressedJet(flashgg::Jet & jt, int npvs, bbggJetRegression::LorentzVector theMET)
{

   if(TMVAReady == 0){
       std::cout << "[bbggJetRegression::GetRegression] You haven't yet setup the regression reader! First do bbggJetRegression::SetupRegression(std::string RegType, std::string RegFile)" << std::endl;
        return;
   }

   Jet_pt = jt.pt();
//   Jet_corr = jt.jecFactor("Uncorrected");
   float JetRawpt = jt.correctedJet("Uncorrected").pt();
   Jet_corr = Jet_pt/JetRawpt;
   if(DEBUG) std::cout << "JETCORR:: " << jt.jecFactor("Uncorrected") << " \t " << Jet_corr << std::endl;
   nPVs = npvs;
   Jet_eta = jt.eta();
   Jet_mt = jt.mt();
   Jet_leadTrackPt = jt.userFloat("leadTrackPt");
   if(Jet_leadTrackPt==-999) Jet_leadTrackPt=0;
   Jet_leptonPtRel = jt.userFloat("softLepRatio");
   if(Jet_leptonPtRel==-999) Jet_leptonPtRel=0;
   Jet_leptonPt = jt.userFloat("softLepPt");
   if(Jet_leptonPt==-999) Jet_leptonPt=0;
   Jet_leptonDeltaR = jt.userFloat("softLepDr");
   if(Jet_leptonDeltaR==-999) Jet_leptonDeltaR=0;
   Jet_totHEF = jt.chargedHadronEnergyFraction() + jt.neutralHadronEnergyFraction();
   Jet_neEmEF = jt.neutralEmEnergyFraction();
   Jet_vtxPt = sqrt( jt.userFloat("vtxPx")*jt.userFloat("vtxPx") + jt.userFloat("vtxPy")*jt.userFloat("vtxPy"));
   if( jt.userFloat("vtxPx") == -999 || jt.userFloat("vtxPy") == -999) Jet_vtxPt = 0;
   Jet_vtxMass = jt.userFloat("vtxMass");
   if(Jet_vtxMass==-999) Jet_vtxMass = 0;
   Jet_vtx3dL = jt.userFloat("vtx3DVal");
   if(Jet_vtx3dL==-999) Jet_vtx3dL = 0;
   Jet_vtxNtrk = jt.userFloat("vtxNTracks");
   if(Jet_vtxNtrk==-999) Jet_vtxNtrk = 0;
   if(Jet_vtx3deL==-999 || Jet_vtx3dL == 0) Jet_vtx3deL = 0;
   else Jet_vtx3deL = Jet_vtx3dL/jt.userFloat("vtx3DSig");
//   if(Jet_vtx3deL==-999) Jet_vtx3deL = 0;
   //MET stuff
   Jet_PFMET = theMET.pt();
   Jet_METDPhi = deltaPhi( jt.phi(), theMET.phi() );
   

   const std::vector< Float_t > & RegressedValues = RegressionReader->EvaluateRegression("BDTG method");
   float corrFactor = RegressedValues[0];

   float oldEnergy = jt.energy();
   float oldPhi = jt.phi();
   TLorentzVector tempVec;
   tempVec.SetPtEtaPhiE(corrFactor*Jet_pt, Jet_eta, oldPhi, corrFactor*oldEnergy);
   bbggJetRegression::LorentzVector temp2;
   temp2.SetPxPyPzE( tempVec.Px(), tempVec.Py(), tempVec.Pz(), tempVec.E());
   jt.setP4(temp2);

   if(DEBUG) {
      std::cout << "[bbggJetRegression::GetRegressedJet] Old pt : " << Jet_pt << " - pt from reg : " << corrFactor*Jet_pt << " - new pt : " << jt.pt() << std::endl;
   }
}

/*
bbggJetRegression::LorentzVector bbggJetRegression::GetRegression(edm::Ptr<flashgg::Jet> jt, float Rho)
{
    bbggJetRegression::LorentzVector NewJet(-1,-1,-1,-1);
    
    if(TMVAReady == 0){
        std::cout << "[bbggJetRegression::GetRegression] You haven't yet setup the regression reader! First do bbggJetRegression::SetupRegression(std::string RegType, std::string RegFile)" << std::endl;
        return NewJet;
    }
    
    //Assigning regression variables:
    Jet_pt = jt->pt();
    Jet_corr = jt->jecFactor("Uncorrected");//jt->correctedJet("Uncorrected").pt();
    rho = Rho;
    Jet_eta = jt->eta();
    Jet_mt = jt->mt();
    
    Jet_leadTrackPt = jt->userFloat("leadTrackPt");//bbggJetRegression::GetLeadingTrackPt(jt);
        
//    float softLeptonPt = jt->softLeptonPt();
//    float softLeptonDR = jt->softLeptonDeltaR();
    
    Jet_leptonPtRel = jt->userFloat("softLepRatio");//softLeptonPt/Jet_pt; //bbggJetRegression::PtRel(ClosestLepton, jt->p4());
    Jet_leptonPt = jt->userFloat("softLepPt");//softLeptonPt; //ClosestLepton.pt();
    Jet_leptonDeltaR = jt->userFloat("softLepDr");//softLeptonDR; //bbggJetRegression::DeltaR(ClosestLepton, jt->p4());
        
//    Jet_chEmEF = jt->chargedEmEnergyFraction();
//    Jet_chHEF = jt->chargedHadronEnergyFraction();
    Jet_neHEF = jt->neutralEmEnergyFraction();
    Jet_neEmEF = jt->neutralHadronEnergyFraction();
    Jet_chMult = jt->chargedHadronMultiplicity();
    Jet_vtxPt = sqrt( jt->userFloat("vtxPx")*jt->userFloat("vtxPx") + jt->userFloat("vtxPy")*jt->userFloat("vtxPy"));
    Jet_vtxMass = jt->userFloat("vtxMass");
    Jet_vtx3dL = jt->userFloat("vtx3DVal");
    Jet_vtxNtrk = jt->userFloat("vtxNTracks");
    Jet_vtx3deL = jt->userFloat("vtx3DSig");
    
    const std::vector< Float_t > & RegressedValues = RegressionReader->EvaluateRegression("BDTG method");
    float newPt = RegressedValues[0];
//    float newEta = RegressedValues[3];
    
    TLorentzVector TempJet;
    TempJet.SetPtEtaPhiM(newPt, jt->eta(), jt->phi(), jt->mass());
    NewJet = bbggJetRegression::LorentzVector(TempJet.Px(), TempJet.Py(), TempJet.Pz(), TempJet.E());
    return NewJet;
}

bbggJetRegression::LorentzVector bbggJetRegression::GetRegression(flashgg::Jet jt, float Rho){
{
    bbggJetRegression::LorentzVector NewJet(-1,-1,-1,-1);
    
    if(TMVAReady == 0){
        std::cout << "[bbggJetRegression::GetRegression] You haven't yet setup the regression reader! First do bbggJetRegression::SetupRegression(std::string RegType, std::string RegFile)" << std::endl;
        return NewJet;
    }
    
    //Assigning regression variables:
    Jet_pt = jt.pt();
    Jet_corr = jt.jecFactor("Uncorrected");//jt.correctedJet("Uncorrected").pt();
    rho = Rho;
    Jet_eta = jt.eta();
    Jet_mt = jt.mt();
    
    Jet_leadTrackPt = jt.userFloat("leadTrackPt");//bbggJetRegression::GetLeadingTrackPt(jt);
        
//    float softLeptonPt = jt.softLeptonPt();
//    float softLeptonDR = jt.softLeptonDeltaR();
    
    Jet_leptonPtRel = jt.userFloat("softLepRatio");//softLeptonPt/Jet_pt; //bbggJetRegression::PtRel(ClosestLepton, jt.p4());
    Jet_leptonPt = jt.userFloat("softLepPt");//softLeptonPt; //ClosestLepton.pt();
    Jet_leptonDeltaR = jt.userFloat("softLepDr");//softLeptonDR; //bbggJetRegression::DeltaR(ClosestLepton, jt.p4());
        
//    Jet_chEmEF = jt.chargedEmEnergyFraction();
//    Jet_chHEF = jt.chargedHadronEnergyFraction();
    Jet_neHEF = jt.neutralEmEnergyFraction();
    Jet_neEmEF = jt.neutralHadronEnergyFraction();
    Jet_chMult = jt.chargedHadronMultiplicity();
    Jet_vtxPt = sqrt( jt.userFloat("vtxPx")*jt.userFloat("vtxPx") + jt.userFloat("vtxPy")*jt.userFloat("vtxPy"));
    Jet_vtxMass = jt.userFloat("vtxMass");
    Jet_vtx3dL = jt.userFloat("vtx3DVal");
    Jet_vtxNtrk = jt.userFloat("vtxNTracks");
    Jet_vtx3deL = jt.userFloat("vtx3DSig");
    
    const std::vector< Float_t > & RegressedValues = RegressionReader->EvaluateRegression("BDTG method");
    float newPt = RegressedValues[0];
//    float newEta = RegressedValues[3];
    
    TLorentzVector TempJet;
    TempJet.SetPtEtaPhiM(newPt, jt.eta(), jt.phi(), jt.mass());
    NewJet = bbggJetRegression::LorentzVector(TempJet.Px(), TempJet.Py(), TempJet.Pz(), TempJet.E());
    return NewJet;
}
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
*/
