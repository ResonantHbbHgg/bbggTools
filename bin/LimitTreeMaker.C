#include "flashgg/bbggTools/interface/bbggLTMaker.h"
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace std;

int Process(string file, string outFile, string mtotMin, string mtotMax,string scale, string photonCR){
    TFile* iFile = new TFile(TString(file), "READ");
    TTree* iTree = (TTree*) iFile->Get("bbggSelectionTree");
    cout << "[LimitTreeMaker:Process] Processing tree with " << iTree->GetEntries() << " entries." << endl;
    string option = outFile;
    option.append(";");
    option.append(mtotMin);
    option.append(";");
    option.append(mtotMax);
    option.append(";");
    option.append(scale);
    option.append(";");
    option.append(photonCR);
    cout << "[LimitTreeMaker:Process] Option parsing: " << option << endl;
    iTree->Process("flashgg/bbggTools/src/bbggLTMaker.cc+", TString(option)); 
    delete iFile;
    return 1;
}

int main(int argc, char *argv[]) {
    string outputFile = "";
    string inputFile = "";
    string outputLocation = "";
    string mtotMax = "10000";
    string mtotMin = "0";
    string scale = "1";
    string photonCR = "0";

    for( int i = 1; i < argc; i++){
	if ( std::string(argv[i]) == "-i"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		inputFile = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-o"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		outputLocation = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-max"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		mtotMax = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-min"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		mtotMin = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-scale"){
		if ( (i+1) == argc) {
			std::cout << "Invalid number of arguments!" << std::endl;
			break;
		}
		scale = string(std::string(argv[i+1]));
		i++;
	}
	else if ( std::string(argv[i]) =="-photonCR"){
		photonCR = "1";
	}
	else {
		cout << "Usage: LimitTreeMaker -i <input list of files> -o <output location> [optional: -min <min mtot> -max <max mtot> -scale <scale factor> -photonCR (do photon control region" << endl;
		return -1;
	}	
    }

    if(inputFile == "" || outputLocation == "") {
	cout << "Usage: LimitTreeMaker -i <input list of files> -o <output location> [optional: -min <min mtot> -max <max mtot> -scale <scale factor> -photonCR (do photon control region" << endl;
	return 0;
    }
/*
    if(argc < 6) {
        cout << "[LimitTreeMaker] Usage: LimitTreeMaker <inputtextfile> <outputLocation> <mtotMax> <mtotMin> <scale> <isPhotonCR>" << endl;
        return 0;
    } else {
        inputFile = argv[1];
        outputLocation = argv[2];
	mtotMin = argv[3];
	mtotMax = argv[4];
        cout << "Opening file: " << argv[1] << endl;
        cout << "Output location: " << argv[2] << endl;
	cout << "Minimum 4-body mass: " << argv[3] << endl;
	cout << "Maximum 4-body mass: " << argv[4] << endl;
	if(argc > 5) {
		scale = argv[5];
		cout << "Scale factor: " << argv[5] << endl;
	}
    }
*/
    ifstream infile(inputFile);
    string line;
    while(getline(infile,line)){
        string rootFileName = "";
        cout << "Processing file: " << line << endl;
        string token;
        stringstream ss(line);
        while(getline(ss, token, '/')){
            if (token.find(".root")!=std::string::npos) rootFileName = token;
        }
        if(rootFileName.find(".root")==std::string::npos) {
            cout << "[LimitTreeMaker] The following line in file is not a root file: " << line << endl;
            continue;
        }
        size_t sLoc = rootFileName.find("bbggSelectionTree_");
        size_t sLen = std::string("bbggSelectionTree_").length();
        rootFileName.replace(sLoc, sLen, "");
        string outF = outputLocation;
        outF.append("/LT_");
        outF.append(rootFileName);
        cout << "Output file: " << outF << endl;
        Process(line, outF, mtotMin, mtotMax, scale, photonCR);
    }
    
    return 0;
}

