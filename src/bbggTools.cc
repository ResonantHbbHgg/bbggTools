#include "flashgg/bbggTools/interface/bbggTools.h"
//FLASHgg files
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"


using namespace std;

double bbggTools::getCHisoToCutValue(edm::Ptr<flashgg::DiPhotonCandidate> dipho, int whichPho)
{
	if(rho_ == -10 ){
		cout << "[bbggTools::getCHisoToCutValue] You have to do tools->setRho(rho)!" << endl;
		return -1;
	}
	if( whichPho > 1 ) {
		std::cout << "[bbggTools::getCHisoToCutValue] You chose the wrong photon!" << std::endl;
		return -1;
	}
	double PFIso = -1, eta = -99;
	if(whichPho == 0) {
		PFIso = dipho->leadingView().pfChIso03WrtChosenVtx();
		eta = dipho->leadingPhoton()->superCluster()->eta();
	}
	if(whichPho == 1) {
		PFIso = dipho->subLeadingView().pfChIso03WrtChosenVtx();
		eta = dipho->subLeadingPhoton()->superCluster()->eta();
	}
	
	double EA = bbggTools::getEA(eta, 0);
	double finalValue = fmax(PFIso - rho_*EA, 0.);
	return finalValue;
}

double bbggTools::getNHisoToCutValue(const flashgg::Photon* pho)
{
	if(rho_ == -10 ){
		cout << "[bbggTools::getNHisoToCutValue] You have to do tools->setRho(rho)!" << endl;
		return -1;
	}
	double PFIso = pho->pfNeutIso03();
	double eta = pho->superCluster()->eta();
	double EA = bbggTools::getEA(eta, 1);
	double extraFactor = 0;
	if(fabs(eta) < 1.479) extraFactor = exp(0.0028*pho->pt()+0.5408);
	if(fabs(eta) > 1.479) extraFactor = 0.01725*pho->pt();
	double finalValue = fmax(PFIso - rho_*EA, 0.) - extraFactor;
	return finalValue;
}

double bbggTools::getPHisoToCutValue(const flashgg::Photon* pho)
{
	if(rho_ == -10 ){
		cout << "[bbggTools::getPHisoToCutValue] You have to do tools->setRho(rho)!" << endl;
		return -1;
	}
	double PFIso = pho->pfPhoIso03();
	double eta = pho->superCluster()->eta();
	double EA = bbggTools::getEA(eta, 2);
	double extraFactor = 0;
	if(fabs(eta) < 1.479) extraFactor = 0.0014*pho->pt();
	if(fabs(eta) > 1.479) extraFactor = 0.0091*pho->pt();
	double finalValue = fmax(PFIso - rho_*EA, 0.) - extraFactor;
	return finalValue;
}

double bbggTools::getEA( float eta, int whichEA){
        if(whichEA < 0 || whichEA > 2){
                std::cout << "WRONG EA TYPE" << std::endl;
                return -1;
        }

        float EA[7][3];

        EA[0][0] = 0.0234; EA[0][1] = 0.0053; EA[0][2] = 0.078;
        EA[1][0] = 0.0189; EA[1][1] = 0.0103; EA[1][2] = 0.0629;
        EA[2][0] = 0.0171; EA[2][1] = 0.0057; EA[2][2] = 0.0264;
        EA[3][0] = 0.0129; EA[3][1] = 0.0070; EA[3][2] = 0.0462;
        EA[4][0] = 0.0110; EA[4][1] = 0.0152; EA[4][2] = 0.0740;
        EA[5][0] = 0.0074; EA[5][1] = 0.0232; EA[5][2] = 0.0924;
        EA[6][0] = 0.0035; EA[6][1] = 0.1709; EA[6][2] = 0.1484;

        float feta = fabs(eta);

        if(feta > 0.000 && feta < 1.000 ) return EA[0][whichEA];
        if(feta > 1.000 && feta < 1.479 ) return EA[1][whichEA];
        if(feta > 1.479 && feta < 2.000 ) return EA[2][whichEA];
        if(feta > 2.000 && feta < 2.200 ) return EA[3][whichEA];
        if(feta > 2.200 && feta < 2.300 ) return EA[4][whichEA];
        if(feta > 2.300 && feta < 2.400 ) return EA[5][whichEA];
        if(feta > 2.400 && feta < 10.00 ) return EA[6][whichEA];

        return -1;
}

double bbggTools::DeltaR(bbggTools::LorentzVector vec1, bbggTools::LorentzVector vec2)
{
	double R2 = (vec1.phi() - vec2.phi())*(vec1.phi() - vec2.phi()) + (vec1.eta() - vec2.eta())*(vec1.eta() - vec2.eta());
	return sqrt(R2);
}

bool bbggTools::isPhoID(edm::Ptr<flashgg::Photon> pho, vector<double> cuts)
{
	if(rho_ == -10 ){
		cout << "[bbggTools::isPhoID] You have to do tools->setRho(rho)!" << endl;
		return -1;
	}
	if(cuts.size() != 3){
		cout << "[bbggTools::isPhoID] ERROR: the input cuts vector must have size three (sieie/hoe/el-veto)" << endl;
		return 0;
	}
	bool isid = true;
	double hoe = pho->hadronicOverEm();
	double sieie = pho->full5x5_sigmaIetaIeta();
	int elveto = pho->passElectronVeto();
	
	if( hoe > cuts[0] ) 	isid = false;
  	if( sieie > cuts[1] ) isid = false;
	if( elveto != (int) cuts[2] ) isid = false;
	
	return isid;
}

bool bbggTools::isPhoISO(edm::Ptr<flashgg::DiPhotonCandidate> dipho, int whichPho, vector<double> cuts)
{
	if(rho_ == -10 ){
		cout << "[bbggTools::isPhoISO] You have to do tools->setRho(rho)!" << endl;
		return -1;
	}
	if(cuts.size() != 3){
		cout << "[bbggTools::isPhoISO] ERROR: the input cuts vector must have size three (ch/nh/ph)" << endl;
		return 0;
	}
	if(whichPho != 0 && whichPho != 1){
		cout << "[bbggTools::isPhoISO] ERROR: whichPho should be 0 (leading photon) or 1 (subleading photon)" << endl;
		return 0;
	}
	
  	double chiso = 0, nhiso = 0, phiso = 0;
	chiso = bbggTools::getCHisoToCutValue( dipho, whichPho);
  	if(whichPho == 0) nhiso = bbggTools::getNHisoToCutValue( dipho->leadingPhoton() );
  	if(whichPho == 1) nhiso = bbggTools::getNHisoToCutValue( dipho->subLeadingPhoton() );
  	if(whichPho == 0) phiso = bbggTools::getPHisoToCutValue( dipho->leadingPhoton() );
	if(whichPho == 1) phiso = bbggTools::getPHisoToCutValue( dipho->subLeadingPhoton() );
	
	bool isiso = true;
	if(chiso > cuts[0]) isiso = false;
	if(nhiso > cuts[1]) isiso = false;
	if(phiso > cuts[2]) isiso = false;
	
	return isiso;
	
}
