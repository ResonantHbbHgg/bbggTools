#include "flashgg/bbggTools/interface/bbggKinFit.h"
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
#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"

using namespace std;

const int DEBUG = 0;

//bool DEBUG = false;

float relFactor = 1.0000000;

bbggKinFit::bbggKinFit(){
    if(DEBUG) cout << "[bbggKinFit::bbggKinFit] Initializing KinFit!" << endl;
    hasptrel = 0;
    hasetarel = 0;
    hasphirel = 0;
    issetup = 0;
//    m1 = TMatrixD(3,3);
//    m2 = TMatrixD(3,3);
//    m1.Zero();
//    m2.Zero();
    _hMass = 125.;
}

double bbggKinFit::DeltaR(bbggKinFit::LorentzVector vec1, bbggKinFit::LorentzVector vec2)
{
	double R2 = (vec1.phi() - vec2.phi())*(vec1.phi() - vec2.phi()) + (vec1.eta() - vec2.eta())*(vec1.eta() - vec2.eta());
	return sqrt(R2);
}

void bbggKinFit::SetJetResolutionParameters(vector<double> etaIntervals, vector<double> pt, vector<double> eta, vector<double> phi)
{
    etaInts = etaIntervals;
    ptRes = pt;
    hasptrel = 1;
    etaRes = eta;
    hasetarel = 1;
    phiRes = phi;
    hasphirel = 1;
    
    // int NVars = eta.size()/etaIntervals.size();
    // if(DEBUG) cout << "[bbggKinFit::SetJetResolutionParameters] Getting resolution parameters! EtaBins: " << etaIntervals.size() << " pt size: " << pt.size() << " eta size: " << eta.size() << " phi size: " << phi.size() << endl;
    // if(pt.size() != (etaIntervals.size()-1)*NVars || eta.size() != (etaIntervals.size()-1)*NVars || phi.size() != (etaIntervals.size()-1)*NVars){
    //     std::cout << "[bbggKinFit] The resolution vectors must have the same size as the intervals in Eta!" << std::endl;
    //     return;
    // } else {
    //     etaInts = etaIntervals;
    //     ptRes = pt;
    //     hasptrel = 1;
    //     etaRes = eta;
    //     hasetarel = 1;
    //     phiRes = phi;
    //     hasphirel = 1;
    // }
}

float bbggKinFit::GetPtResolution(bbggKinFit::LorentzVector j)
{
    if(DEBUG) cout << "[bbggKinFit::GetPtResolution] Jet pt: " << j.pt() << endl;
    float Et = j.pt();
    float eta = fabs(j.eta());
    float a=-1, b=-1, c=-1, d=-1;
    for(unsigned int i = 1; i < etaInts.size(); i++){
	if(DEBUG) std::cout << "eta: " << eta << " between: " << etaInts[i-1] << " and " << etaInts[i] << endl;
        if(eta > etaInts[i-1] && eta < etaInts[i]){
            unsigned int I = (i-1)*3;
            a = ptRes[I];//0, 3, 6
            b = ptRes[I+1];//1, 4, 7
            c = ptRes[I+2];//2, 5, 8
            d = ptRes[I+3];//2, 5, 8
        }
    }
    if(a == -1) {
        std::cout << "[bbggKinFit::GetPtResolution] Unset Kin Fit parameter!" << std::endl;
        return -1;
    }
    
    //Kin fit example function
    //return (a * a) + (b * b) * Et + (c * c) * Et * Et;
    //RunI function
    float ptResolution = ( a * fabs(a) ) + ( b * b ) * pow(Et,d+1) + ( c * c )* Et * Et;
    return ptResolution*relFactor;
}

float bbggKinFit::GetEtaResolution(bbggKinFit::LorentzVector j)
{
    if(DEBUG) cout << "[bbggKinFit::GetEtaResolution] Jet pt: " << j.pt() << endl;
    float Et = j.pt();
    float eta = fabs(j.eta());
    float a=-1, b=-1, c=-1, d=-1, e=-1;
    int NVars = etaRes.size()/(etaInts.size()-1);
    for(unsigned int i = 1; i < etaInts.size(); i++){
        if(eta > etaInts[i-1] && eta < etaInts[i]){
            unsigned int I = (i-1)*NVars;
            a = etaRes[I];//0, 3, 6
            b = etaRes[I+1];//1, 4, 7
            c = etaRes[I+2];//2, 5, 8
            d = etaRes[I+3];//2, 5, 8
            e = etaRes[I+4];//2, 5, 8
        }
    }
    if(a == -1) {
        std::cout << "[bbggKinFit::GetEtaResolution] Unset Kin Fit parameter!" << std::endl;
        return -1;
    }
    
    //Kin fit example function
    //return a/(Et * Et) + b/Et + c;
    //RunI function
    float etaResolution = sqrt(a*a/(Et * Et) + b*b/Et + c*c) + d/Et + e/pow(Et,1.5);
    return (etaResolution*etaResolution)*relFactor;
}

float bbggKinFit::GetPhiResolution(bbggKinFit::LorentzVector j)
{
    if(DEBUG) cout << "[bbggKinFit::GetPhiResolution] Jet pt: " << j.pt() << endl;
    float Et = j.pt();
    float eta = fabs(j.eta());
    float a=-1, b=-1, c=-1, d=-1, e=-1;
    for(unsigned int i = 1; i < etaInts.size(); i++){
        if(eta > etaInts[i-1] && eta < etaInts[i]){
            unsigned int I = (i-1)*3;
            a = phiRes[I];//0, 3, 6
            b = phiRes[I+1];//1, 4, 7
            c = phiRes[I+2];//2, 5, 8
            d = phiRes[I+3];//2, 5, 8
            e = phiRes[I+4];//2, 5, 8
        }
    }
    if(a == -1) {
        std::cout << "[bbggKinFit::GetPhiResolution] Unset Kin Fit parameter!" << std::endl;
        return -1;
    }
    
    //Kin fit example function
    //return a/(Et * Et) + b/Et + c;
    //RunI function
    float phiResolution = sqrt(a*a/(Et * Et) + b*b/Et + c*c) + d/Et + e/pow(Et,1.5);
    return (phiResolution*phiResolution)*relFactor;
}

void bbggKinFit::KinematicFit(bbggKinFit::LorentzVector j1, bbggKinFit::LorentzVector j2)
{
    if(DEBUG) cout << "[bbggKinFit::KinematicFit] Starting fit 1!" << endl;
    TLorentzVector jt1(j1.px(), j1.py(), j1.pz(), j1.energy());
    TLorentzVector jt2(j2.px(), j2.py(), j2.pz(), j2.energy());
    if(DEBUG) cout << "[bbggKinFit::KinematicFit] Starting fit 2!" << endl;
    //Setup covariance matrices:
    TMatrixD m1(3,3);
    m1.Zero();
    TMatrixD m2(3,3);    
    m2.Zero();

    m1(0,0) = bbggKinFit::GetPtResolution( j1);
    m1(1,1) = bbggKinFit::GetEtaResolution(j1);
    m1(2,2) = bbggKinFit::GetPhiResolution(j1);

    m2(0,0) = bbggKinFit::GetPtResolution( j2);
    m2(1,1) = bbggKinFit::GetEtaResolution(j2);
    m2(2,2) = bbggKinFit::GetPhiResolution(j2);

    if(DEBUG) cout << "[bbggKinFit::KinematicFit] Starting fit 3!" << endl;
    TFitParticleEtEtaPhi *jet1 = new TFitParticleEtEtaPhi( "Jet1", "Jet1", &jt1, &m1 );
    TFitParticleEtEtaPhi *jet2 = new TFitParticleEtEtaPhi( "Jet2", "Jet2", &jt2, &m2 );
    
    TFitConstraintM *mCons1 = new TFitConstraintM( "HMassConstraint", "HMass-Constraint", 0, 0 , _hMass);
    mCons1->addParticles1( jet1, jet2 );
    
    if(DEBUG) cout << "[bbggKinFit::KinematicFit] Setting up!" << endl;
    TKinFitter* fitter = new TKinFitter("fitter", "fitter");
    fitter->addMeasParticle( jet1 );
    fitter->addMeasParticle( jet2 );
    fitter->addConstraint( mCons1 );
    
    //Set convergence criteria
    fitter->setMaxNbIter( 30 );
    fitter->setMaxDeltaS( 1e-2 );
    fitter->setMaxF( 1e-1 );
    fitter->setVerbosity(1);
    
    //fit
    if(DEBUG) cout << "[bbggKinFit::KinematicFit] Start fit!" << endl;
    fitter->fit();
    TLorentzVector newJet1(*jet1->getCurr4Vec());
    TLorentzVector newJet2(*jet2->getCurr4Vec());
    Jet1 = bbggKinFit::LorentzVector( newJet1.Px(), newJet1.Py(), newJet1.Pz(), newJet1.E());
    Jet2 = bbggKinFit::LorentzVector( newJet2.Px(), newJet2.Py(), newJet2.Pz(), newJet2.E());
    if(DEBUG) cout << "[bbggKinFit::KinematicFit] Ended fit!" << endl;
    
    //Memory!
    delete jet1;
    delete jet2;
    delete mCons1;
    delete fitter;
}
