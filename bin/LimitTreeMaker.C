#include "flashgg/bbggTools/interface/bbggLTMaker.h"
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace std;

int Process(string file, string outFile, string mtotMin, string mtotMax,string cat = "-1"){
    TFile* iFile = new TFile(TString(file), "READ");
    TTree* iTree = (TTree*) iFile->Get("bbggSelectionTree");
    cout << "[LimitTreeMaker:Process] Processing tree with " << iTree->GetEntries() << " entries." << endl;
    string option = outFile;
    option.append(";");
    option.append(mtotMin);
    option.append(";");
    option.append(mtotMax);
    option.append(";");
    option.append(cat);
    cout << "[LimitTreeMaker:Process] Option parsing: " << option << endl;
    iTree->Process("flashgg/bbggTools/src/bbggLTMaker.cc+", TString(option)); 
    delete iFile;
    return 1;
}

int main(int argc, char *argv[]) {
    string outputFile = "";
    string inputFile = "";
    string outputLocation = "";
    string mtotMax = "";
    string mtotMin = "";
    if(argc < 5) {
        cout << "[LimitTreeMaker] Usage: LimitTreeMaker <inputtextfile> <outputLocation> <mtotMax> <mtotMin>" << endl;

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
    }
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
        Process(line, outF, mtotMin, mtotMax);
    }
    
    return 0;
}

