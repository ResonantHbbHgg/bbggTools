#ifndef FLASHgg_genericMVA_h
#define FLASHgg_genericMVA_h

// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <string>
#include <utility>
#include <cmath>
//root include files
#include "TFile.h"
#include "TMVA/Reader.h"

using namespace std;

class genericMVA{
public:
	genericMVA() : TMVAReady(0) {}
	void setupMVA(std::string File, std::vector<std::string> inVars);
	~genericMVA() {}

        float mvaDiscriminants(std::map<std::string,float>) ;

private:
    	TMVA::Reader *Reader;
	bool TMVAReady;
    //MVA variables
        std::map<std::string, float> mvaVars;
        std::vector<std::string> orderedVars;

};

#endif
