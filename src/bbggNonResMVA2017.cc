#include "flashgg/bbggTools/interface/bbggNonResMVA2017.h"
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


void bbggNonResMVA2017::SetupNonResMVA(std::string MVAFile, std::vector<std::string> inVars){
  TMVAReady = 1;
  MVAReader = new TMVA::Reader();

  for (unsigned int iv = 0; iv < inVars.size(); iv++)
    mvaVars[inVars[iv]] = -10;

  orderedVars = inVars;

  std::cout << "[bbggNonResMVA2017::SetupNonResMVA] NonResMVA set with the following variables (attention, the order is important): " << std::endl;
  //    for ( std::map<std::string,float>::iterator it = mvaVars.begin(); it != mvaVars.end(); it++ ) {
  for ( unsigned int vv = 0; vv < orderedVars.size(); vv++) {
    std::map<std::string,float>::iterator it = mvaVars.find(orderedVars[vv]);
    //      MVAReader->AddVariable(it->first, mvaVars[it->first]);
    MVAReader->AddVariable(it->first, &it->second);
    std::cout << "\t" << it->first << std::endl;
  }
  MVAReader->BookMVA("BDT", MVAFile); 

}

double bbggNonResMVA2017::DeltaR(bbggNonResMVA2017::LorentzVector vec1, bbggNonResMVA2017::LorentzVector vec2)
{
  //double R2 = (vec1.phi() - vec2.phi())*(vec1.phi() - vec2.phi()) + (vec1.eta() - vec2.eta())*(vec1.eta() - vec2.eta());
  //return sqrt(R2);
  return deltaR(vec1, vec2);
}

std::vector<float> bbggNonResMVA2017::mvaDiscriminants(std::map<std::string,float> params)
{

  std::vector<float> mvaDis;

  if(TMVAReady == 0){
    std::cout << "[bbggNonResMVA2017::mvaDiscriminants] You haven't yet setup the NonResMVA reader!" << std::endl;
    return mvaDis;
  }

  if( params.size() != mvaVars.size()) {
    std::cout << "[bbggNonResMVA2017::mvaDiscriminants] Your input list does not match the variables list in the MVA method!" << std::endl;
    return mvaDis;
  }

  for ( std::map<std::string,float>::iterator it = mvaVars.begin(); it != mvaVars.end(); it++ ) {
    std::map<std::string,float>::iterator toFind = params.find(it->first);
    if ( toFind == mvaVars.end() ) {
      std::cout << "[bbggNonResMVA2017::mvaDiscriminants] Variable that is in the training ## " << it->first << " not found in the parameters map provided!!" <<std::endl;
      return mvaDis;
    } else {
      mvaVars[it->first] = params[it->first];
    }
  }

  float MVAdis = MVAReader->EvaluateMVA("BDT");

  mvaDis.push_back(MVAdis);
   
  return mvaDis;
}

