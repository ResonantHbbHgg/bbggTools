#include "flashgg/bbggTools/interface/bbggLTMaker.h"
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace std;

int Process(string file, string outFile, string cat){
    TFile* iFile = new TFile(TString(file), "READ");
    TTree* iTree = (TTree*) iFile->Get("bbggSelectionTree");
    cout << "[LimitTreeMaker:Process] Processing tree with " << iTree->GetEntries() << " entries." << endl;
    string option = outFile;
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
    string category = "";
    string outputLocation = "";
    if(argc < 4) {
        cout << "[LimitTreeMaker] Usage: LimitTreeMaker <inputtextfile> <category: 0(no btag)/1(loose btag)/2(tight btag)> <outputLocation>" << endl;
        return 0;
    } else {
        inputFile = argv[1];
        category = argv[2];
        outputLocation = argv[3];
        cout << "Opening file: " << argv[1] << endl;
        cout << "Category used: " << argv[2] << endl;
        cout << "Output location: " << argv[3] << endl;
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
        outF.append(category);
        outF.append("bTag_");
        outF.append(rootFileName);
        cout << "Output file: " << outF << endl;
        Process(line, outF, category);
    }
    
    return 0;
}

