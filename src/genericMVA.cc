#include "flashgg/bbggTools/interface/genericMVA.h"
//FLASHgg files

using namespace std;

const int DEBUG = 0;

//bool DEBUG = false;


void genericMVA::setupMVA(std::string File, std::vector<std::string> inVars){
    TMVAReady = 1;
    Reader = new TMVA::Reader();

    for (unsigned int iv = 0; iv < inVars.size(); iv++)
      mvaVars[inVars[iv]] = -10;

    orderedVars = inVars;

    std::cout << "[MVA::SetupMVA] MVA set with the following variables (attention, the order is important): " << std::endl;
//    for ( std::map<std::string,float>::iterator it = mvaVars.begin(); it != mvaVars.end(); it++ ) {
    for ( unsigned int vv = 0; vv < orderedVars.size(); vv++) {
      std::map<std::string,float>::iterator it = mvaVars.find(orderedVars[vv]);
      Reader->AddVariable(it->first, &it->second);
      std::cout << "\t" << it->first << std::endl;
    }
    Reader->BookMVA("BDT", File); 

}

float genericMVA::mvaDiscriminants(std::map<std::string,float> params)
{

  float mvaDis = -10.;

   if(TMVAReady == 0){
       std::cout << "[bbggNonResMVA::mvaDiscriminants] You haven't yet setup the NonResMVA reader!" << std::endl;
        return mvaDis;
   }

   if( params.size() != mvaVars.size()) {
       std::cout << "[bbggNonResMVA::mvaDiscriminants] Your input list does not match the variables list in the MVA method!" << std::endl;
       return mvaDis;
   }

   for ( std::map<std::string,float>::iterator it = mvaVars.begin(); it != mvaVars.end(); it++ ) {
      std::map<std::string,float>::iterator toFind = params.find(it->first);
      if ( toFind == mvaVars.end() ) {
        std::cout << "[bbggNonResMVA::mvaDiscriminants] Variable that is in the training ## " << it->first << " not found in the parameters map provided!!" <<std::endl;
        return mvaDis;
      } else {
        mvaVars[it->first] = params[it->first];
      }
   }

   float dis = Reader->EvaluateMVA("BDT");
   
   return dis;
}

