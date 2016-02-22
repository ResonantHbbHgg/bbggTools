#include "flashgg/bbggTools/interface/bbggTools.h"
//FLASHgg files
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"


using namespace std;

bool DEBUG = 0;

std::vector<int> bbggTools::TriggerSelection(std::vector<std::string> myTriggers, const edm::TriggerNames &names, edm::Handle<edm::TriggerResults> triggerBits)
{
    std::vector<int> triggerResults;
    for(unsigned int j = 0; j < myTriggers.size(); j++)
    {
        int accepted = 0;
        for ( unsigned int i = 0; i < triggerBits->size(); i++)
        {
            if(DEBUG) std::cout << "[bbggTools::TriggerSelection] Trigger name: " << names.triggerName(i) << " \t Decision: " << triggerBits->accept(i) << std::endl;
            if((names.triggerName(i)).find(myTriggers[j]) != std::string::npos )
            {
                if(triggerBits->accept(i) == 1){
                    accepted = 1;
                }
            }
        }
        triggerResults.push_back(accepted);
    }
    return triggerResults;
}


bool bbggTools::isJetID(edm::Ptr<flashgg::Jet> jet)
{
    double NHF = jet->neutralHadronEnergyFraction();
    double NEMF = jet->neutralEmEnergyFraction();
    double NumConst = jet->chargedMultiplicity()+jet->neutralMultiplicity();
    double CHF = jet->chargedHadronEnergyFraction();
    double CHM = jet->chargedMultiplicity();
    double CEMF = jet->chargedEmEnergyFraction();
    double NNP = jet->neutralMultiplicity();
    
    if( fabs(jet->eta()) < 3.0 )
    {
        if(NHF > 0.99) return 0;
        if(NEMF > 0.99) return 0;
        if(NumConst < 2) return 0;
        if( fabs(jet->eta()) < 2.4 ){
            if(CHF < 0) return 0;
            if(CHM < 0) return 0;
            if(CEMF > 0.99) return 0;
        }    
    }
    if ( fabs(jet->eta()) > 3.0 )
    {
        if(NEMF > 0.90) return 0;
        if(NNP < 10) return 0;
    }
    
    return 1;
}

std::map<int, vector<double> > bbggTools::getWhichID (std::string wpoint)
{
    if( wpoint.find("loose") != std::string::npos) {
        return _phoIDloose;
    }
    if( wpoint.find("medium") != std::string::npos) {
        return _phoIDmedium;
    }
    if( wpoint.find("tight") != std::string::npos) {
        return _phoIDtight;
    }
    std::cout << "[bbggTools::getWhichID] " << wpoint << " has to be either loose, medium or tight!" << std::endl;
    std::map<int, vector<double> > empty;
    return empty;
}
std::map<int, vector<double> > bbggTools::getWhichISO (std::string wpoint)
{
    if( wpoint.find("loose") != std::string::npos) {
        return _phoISOloose;
    }
    if( wpoint.find("medium") != std::string::npos) {
        return _phoISOmedium;
    }
    if( wpoint.find("tight") != std::string::npos) {
        return _phoISOtight;
    }
    std::cout << "[bbggTools::getWhichISO] " << wpoint << " has to be either loose, medium or tight!" << std::endl;
    std::map<int, vector<double> > empty;
    return empty;
}

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
		PFIso = dipho->leadingView()->pfChIso03WrtChosenVtx();
		eta = dipho->leadingPhoton()->superCluster()->eta();
	}
	if(whichPho == 1) {
		PFIso = dipho->subLeadingView()->pfChIso03WrtChosenVtx();
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

double bbggTools::getNHisoToCutValue(const flashgg::Photon* pho, vector<double> nhCorr)
{
	if(rho_ == -10 ){
		cout << "[bbggTools::getNHisoToCutValue] You have to do tools->setRho(rho)!" << endl;
		return -1;
	}
	if( nhCorr.size() < 2 ) {
		cout << "[bbggTools::getNHisoToCutValue] nhCorr vector not initialized correctly!" << endl;
		return -1;
	}
	double PFIso = pho->pfNeutIso03();
	double eta = pho->superCluster()->eta();
	double EA = bbggTools::getEA(eta, 1);
//	double extraFactor = exp(nhCorr[0]*pho->pt()+nhCorr[1]);
	double extraFactor = nhCorr[0]*pho->pt()+nhCorr[1]*(pho->pt())*(pho->pt());
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

double bbggTools::getPHisoToCutValue(const flashgg::Photon* pho, vector<double> phCorr)
{
	if(rho_ == -10 ){
		cout << "[bbggTools::getPHisoToCutValue] You have to do tools->setRho(rho)!" << endl;
		return -1;
	}
	if( phCorr.size() < 2 ) {
		cout << "[bbggTools::getPHisoToCutValue] phCorr vector not initialized correctly!" << endl;
		return -1;
	}
	double PFIso = pho->pfPhoIso03();
	double eta = pho->superCluster()->eta();
	double EA = bbggTools::getEA(eta, 2);
	double extraFactor = phCorr[0]*pho->pt() + phCorr[1];
	double finalValue = fmax(PFIso - rho_*EA, 0.) - extraFactor;
	return finalValue;
}

double bbggTools::getEA( float eta, int whichEA){
        if(whichEA < 0 || whichEA > 2){
                std::cout << "WRONG EA TYPE" << std::endl;
                return -1;
        }

        float EA[7][3];

        //Spring 15 50ns values (25ns values still not available)
//	EA[0][0] = 0.0157; EA[0][1] = 0.0143; EA[0][2] = 0.0725;
//      EA[1][0] = 0.0143; EA[1][1] = 0.0210; EA[1][2] = 0.0604;
//      EA[2][0] = 0.0115; EA[2][1] = 0.0147; EA[2][2] = 0.0320;
//      EA[3][0] = 0.0094; EA[3][1] = 0.0082; EA[3][2] = 0.0512;
//      EA[4][0] = 0.0095; EA[4][1] = 0.0124; EA[4][2] = 0.0766;
//      EA[5][0] = 0.0068; EA[5][1] = 0.0186; EA[5][2] = 0.0949;
//      EA[6][0] = 0.0053; EA[6][1] = 0.0320; EA[6][2] = 0.1160;
		
        //Spring 15 25ns values ("Please note that the Charged Hadron Isolation is NOT corrected with Effective Areas in the SP15 25ns ID")
	EA[0][0] = 0.0; EA[0][1] = 0.0599; EA[0][2] = 0.1271;
        EA[1][0] = 0.0; EA[1][1] = 0.0819; EA[1][2] = 0.1101;
        EA[2][0] = 0.0; EA[2][1] = 0.0696; EA[2][2] = 0.0756;
        EA[3][0] = 0.0; EA[3][1] = 0.0360; EA[3][2] = 0.1175;
        EA[4][0] = 0.0; EA[4][1] = 0.0360; EA[4][2] = 0.1498;
        EA[5][0] = 0.0; EA[5][1] = 0.0462; EA[5][2] = 0.1857;
        EA[6][0] = 0.0; EA[6][1] = 0.0656; EA[6][2] = 0.2183;

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
//	int elveto = pho->passElectronVeto();
	
	if( hoe > cuts[0] ) 	isid = false;
  	if( sieie > cuts[1] ) isid = false;
//	if( elveto != (int) cuts[2] ) isid = false;
	
	return isid;
}

bool bbggTools::isPhoID(const flashgg::Photon* pho, vector<double> cuts)
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
	
	if( hoe > cuts[0] ) 	isid = false;
  	if( sieie > cuts[1] ) isid = false;
	
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
	const flashgg::Photon* pho = (whichPho) ? dipho->subLeadingPhoton() : dipho->leadingPhoton();
	nhiso = bbggTools::getNHisoToCutValue( pho );
	phiso = bbggTools::getPHisoToCutValue( pho );
	bool isiso = true;
	if(chiso > cuts[0]) isiso = false;
	if(nhiso > cuts[1]) isiso = false;
	if(phiso > cuts[2]) isiso = false;
	
	return isiso;
	
}

bool bbggTools::isPhoISO(edm::Ptr<flashgg::DiPhotonCandidate> dipho, int whichPho, vector<double> cuts, vector<double> nhCorr, vector<double> phCorr)
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
	const flashgg::Photon* pho = (whichPho) ? dipho->subLeadingPhoton() : dipho->leadingPhoton();
	nhiso = bbggTools::getNHisoToCutValue( pho, nhCorr );
	phiso = bbggTools::getPHisoToCutValue( pho, phCorr );
	bool isiso = true;
	if(chiso > cuts[0]) isiso = false;
	if(nhiso > cuts[1]) isiso = false;
	if(phiso > cuts[2]) isiso = false;
	
	return isiso;
	
}


bool bbggTools::CheckCuts()
{
	if(	_PhotonPtOverDiPhotonMass.size() == 0 ||
		_PhotonEta.size() == 0 ||
		_PhotonDoID.size() == 0 ||
		_PhotonDoISO.size() == 0 ||
		_PhotonR9.size() == 0 ||
		_PhotonElectronVeto.size() == 0 ||
		_DiPhotonPt.size() == 0 ||
		_DiPhotonEta.size() == 0 ||
		_DiPhotonMassWindow.size() == 0 ||
		_JetPt.size() == 0 ||
		_JetEta.size() == 0 ||
		_JetBDiscriminant.size() == 0 ||
		_JetDrPho.size() == 0 ||
		_JetDoPUID.size() == 0 ||
		_DiJetPt.size() == 0 ||
		_DiJetEta.size() == 0 ||
		_DiJetMassWindow.size() == 0 ||
		_CandidateMassWindow.size() == 0 ||
		_CandidatePt.size() == 0 ||
		_CandidateEta.size() == 0 ||
		_CandidatesDeltaR.size() == 0 ||
        _PhotonDoElectronVeto.size() == 0 ) {
		return 0;
	} else {
		return 1;
	}
}

edm::Ptr<flashgg::DiPhotonCandidate> bbggTools::GetSelected_diphoCandidate()
{
	if(!hasDiPho){
		std::cout << "[bbggTools::GetSelected_diphoCandidate] ERROR: diphoCandidate has not been set yet. Run bbggTools::AnalysisSelection. If you did, then the event failed the selection requirements." << std::endl;
		return diphoCandidate;
	}
	return diphoCandidate;
}

edm::Ptr<flashgg::Jet> bbggTools::GetSelected_leadingJetCandidate()
{
	if(!hasLeadJet){
		std::cout << "[bbggTools::GetSelected_leadingJetCandidate] ERROR: leadingJetCandidate has not been set yet. Run bbggTools::AnalysisSelection. If you did, then the event failed the selection requirements." << std::endl;
		return leadingJetCandidate;
	}
	return leadingJetCandidate;
}

edm::Ptr<flashgg::Jet> bbggTools::GetSelected_subleadingJetCandidate()
{
	if(!hasSubJet){
		std::cout << "[bbggTools::GetSelected_subleadingJetCandidate] ERROR: subleadingJetCandidate has not been set yet. Run bbggTools::AnalysisSelection. If you did, then the event failed the selection requirements." << std::endl;
		return subleadingJetCandidate;
	}
	return subleadingJetCandidate;
}

vector<edm::Ptr<flashgg::DiPhotonCandidate>> bbggTools::DiPhotonKinematicSelection( vector<edm::Ptr<flashgg::DiPhotonCandidate>> diphoCol)
{
    vector<edm::Ptr<flashgg::DiPhotonCandidate>> PreselDiPhotons;
    bool isValidDiPhotonCandidate = false;
    //Begin DiPhoton Loop/Selection -----------------------------------------------------------
    for( unsigned int diphoIndex = 0; diphoIndex < diphoCol.size(); diphoIndex++ )
    {

     if(DEBUG) std::cout << "[bbggTools::AnalysisSelection] Photon loop 1... " << diphoCol.size() << "\t" << diphoIndex << std::endl;

         if(_DiPhotonOnlyFirst && diphoIndex > 0 ) break;

         edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphoCol[ diphoIndex ];

         double dipho_pt = dipho->pt();
         double dipho_eta = dipho->eta();
         double dipho_mass = dipho->mass();

     //Cuts on diphoton object
//         if(dipho_mass < _DiPhotonMassWindow[0] || dipho_mass > _DiPhotonMassWindow[1]) continue;
         if(fabs(dipho_eta) > _DiPhotonEta[0] ) continue;
         if(dipho_pt < _DiPhotonPt[0] ) continue;

     //Cuts on photons
         double pho1_pt = dipho->leadingPhoton()->pt();
         if( pho1_pt < dipho_mass*_PhotonPtOverDiPhotonMass[0] ) continue;

         double pho2_pt = dipho->subLeadingPhoton()->pt();
         if( pho2_pt < dipho->mass()*_PhotonPtOverDiPhotonMass[1] ) continue;

         double pho1_eta = dipho->leadingPhoton()->superCluster()->eta();
         if( fabs(pho1_eta) > _PhotonEta[1] ) continue;

         double pho2_eta = dipho->subLeadingPhoton()->superCluster()->eta();
         if( fabs(pho2_eta) > _PhotonEta[1] ) continue;

	 PreselDiPhotons.push_back(dipho);

     }

     return PreselDiPhotons;

}

vector<edm::Ptr<flashgg::DiPhotonCandidate>> bbggTools::DiPhotonIDSelection( vector<edm::Ptr<flashgg::DiPhotonCandidate>> diphoCol)
{
    vector<edm::Ptr<flashgg::DiPhotonCandidate>> PreselDiPhotons;
    bool isValidDiPhotonCandidate = false;
    //Begin DiPhoton Loop/Selection -----------------------------------------------------------
   for( unsigned int diphoIndex = 0; diphoIndex < diphoCol.size(); diphoIndex++ )
   {
    //
         if(DEBUG) std::cout << "[bbggTools::AnalysisSelection] Photon loop 2..." << std::endl;
         edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphoCol[ diphoIndex ];

         double dipho_pt = dipho->pt();
         double dipho_eta = dipho->eta();
         double dipho_mass = dipho->mass();

         int pho_elv[2];
         pho_elv[0] = dipho->leadingPhoton()->passElectronVeto();
         pho_elv[1] = dipho->subLeadingPhoton()->passElectronVeto();

     //Cuts on photons
 	     double pho1_pt = dipho->leadingPhoton()->pt();
	     double pho2_pt = dipho->subLeadingPhoton()->pt();
 	     double pho1_eta = dipho->leadingPhoton()->superCluster()->eta();
	     double pho2_eta = dipho->subLeadingPhoton()->superCluster()->eta();

         int pho_ids[2];
         pho_ids[0] = 1;
         pho_ids[1] = 1;
         if(DEBUG) std::cout << "[bbggTools::AnalysisSelection] Photon loop 3..." << std::endl; 
         for( int whichPho = 0; whichPho < 2; whichPho++)
         {
             if(DEBUG) cout << "WHICH PHO: " << whichPho << " _PhotonDoElectronVeto[whichPho] " << _PhotonDoElectronVeto[whichPho] << endl; 
             if( _PhotonDoElectronVeto[whichPho] ) {
                 pho_ids[whichPho] = pho_elv[whichPho];
                 if(DEBUG) std::cout << "DOING VETO! pass: " << pho_elv[whichPho] << " " << whichPho << " "<< dipho->leadingPhoton()->passElectronVeto() << " "<< dipho->subLeadingPhoton()->passElectronVeto() << endl; 
             }
             int pho1Index = 0; //Index 0 = barrel, 1 = endcap
             double pho_eta = (whichPho) ? fabs(pho2_eta) : fabs(pho1_eta);
             if( pho_eta > _PhotonEta[0] ) pho1Index = 1;	 
             if( _PhotonDoID[whichPho] )
             { 
                 std::map<int, vector<double> > theIDWP = bbggTools::getWhichID(_phoWhichID[whichPho]);
                 if(theIDWP.size() < 1) {
                     std::cout << "[bbggTools::AnalysisSelection] _phoWhichID[ " << whichPho << "] = " << _phoWhichID[whichPho] << " has to be either loose, medium or tight!" << std::endl;
                     break;
                 }
                 int pho1_id = (whichPho) ? bbggTools::isPhoID(dipho->subLeadingPhoton(), theIDWP[pho1Index]) : bbggTools::isPhoID(dipho->leadingPhoton(), theIDWP[pho1Index]) ;
                 if (!pho1_id) pho_ids[whichPho] = 0;
             }
             if(_PhotonDoISO[whichPho])
             {
                 std::map<int, vector<double> > theISOWP = bbggTools::getWhichISO(_phoWhichISO[whichPho]);
                 if(theISOWP.size() < 1) {
                     std::cout << "[bbggTools::AnalysisSelection] _phoWhichISO[ " << whichPho << "] = " << _phoWhichISO[whichPho] << " has to be either loose, medium or tight!" << std::endl;
                     break;
                 }
                 int pho1_id = bbggTools::isPhoISO(dipho, whichPho, theISOWP[pho1Index], _nhCorr[pho1Index], _phCorr[pho1Index]);
                 if (!pho1_id) pho_ids[whichPho] = 0;
             }//here
         }
         if(DEBUG) std::cout << "[bbggTools::AnalysisSelection] After Photon loop..." << std::endl;

         if(pho_ids[0] == true && pho_ids[1] == true){
             isValidDiPhotonCandidate = true;
             //CandVtx = dipho->vtx();
             diphoCandidate = dipho;
             hasDiPho = true;
             _isSignal = 1;
             _isPhotonCR = 0;
             PreselDiPhotons.push_back(dipho);
         }
     }
     return PreselDiPhotons;
}

std::vector<edm::Ptr<flashgg::Jet>> bbggTools::DiJetSelection(std::vector<edm::Ptr<flashgg::Jet>> Jets)
{


    edm::Ptr<flashgg::Jet> jet1, jet2;
    std::vector<edm::Ptr<flashgg::Jet>> SelDijet;
    bbggTools::LorentzVector DiJet(0,0,0,0);
//    double dijetPt_ref = 0;
    double sumbtag_ref = 0;
    double sumpt_ref = 0;
    bool hasDiJet = false;
    sumbtag_ref = sumbtag_ref;
    sumpt_ref = sumpt_ref;

	if(DEBUG) std::cout << "Jet sorting... " << std::endl;
    for(unsigned int iJet = 0; iJet < Jets.size(); iJet++)
    {
 		unsigned int isbjet = 0;
 		if( Jets[iJet]->bDiscriminator(_bTagType) > _JetBDiscriminant[1] ) isbjet = 1;
 		for(unsigned int jJet = iJet+1; jJet < Jets.size(); jJet++)
 		{
 	  	 	unsigned int isbjet2 = 0;
 	  		if( Jets[jJet]->bDiscriminator(_bTagType) > _JetBDiscriminant[1] ) isbjet2 = 1;
 	  		unsigned int totalbjet = isbjet + isbjet2;
 	  		if(_n_bJets && totalbjet != _n_bJets) continue;

 	  		bbggTools::LorentzVector dijet = Jets[iJet]->p4() + Jets[jJet]->p4();
//	  		if(dijet.mass() < _DiJetMassWindow[0] || dijet.mass() > _DiJetMassWindow[1]) continue;
//			if( bbggTools::DeltaR( dijet, diphoCandidate->p4() ) < _CandidatesDeltaR[0] ) continue;
						
			double sumbtag = Jets[iJet]->bDiscriminator(_bTagType) + Jets[jJet]->bDiscriminator(_bTagType);

			//TESTCUT
//			double prodbtag = (Jets[iJet]->bDiscriminator(_bTagType))*(Jets[jJet]->bDiscriminator(_bTagType));
//			if(prodbtag < 0.3)
//				continue;
			if( bbggTools::DeltaR(Jets[iJet]->p4(), Jets[jJet]->p4()) < 0.5)
				continue;

			double sumpt = Jets[iJet]->pt() + Jets[jJet]->pt();
// 	  		if(dijet.pt() > dijetPt_ref && dijet.pt() > _DiJetPt[0] && fabs(dijet.Eta()) < _DiJetEta[0] )
 	  		if(sumbtag > sumbtag_ref && dijet.pt() > _DiJetPt[0] && fabs(dijet.Eta()) < _DiJetEta[0] )
// 	  		if(sumpt > sumpt_ref && dijet.pt() > _DiJetPt[0] && fabs(dijet.Eta()) < _DiJetEta[0] )
 	  		{
				hasDiJet = true;
//				dijetPt_ref = dijet.pt();
				sumbtag_ref = sumbtag;
				sumpt_ref = sumpt;
				if( Jets[iJet]->pt() > Jets[jJet]->pt() ) {
					jet1 = Jets.at(iJet);
					jet2 = Jets.at(jJet);
				} else {
					jet2 = Jets.at(iJet);
					jet1 = Jets.at(jJet);
				} 
			}
		}
	}

	if(hasDiJet){
		SelDijet.push_back(jet1);
		SelDijet.push_back(jet2);
	}
	return SelDijet;

}

std::vector<edm::Ptr<flashgg::Jet>> bbggTools::JetPreSelection(JetCollectionVector jetsCol, edm::Ptr<flashgg::DiPhotonCandidate> diphoCandidate)
{

    //Begin Jets Loop/Selection ------------------------------------------------------------
    std::vector<edm::Ptr<flashgg::Jet>> Jets;
    if(DEBUG) std::cout << "Begin Jet selection..." << std::endl;
    unsigned int jetCollectionIndex = diphoCandidate->jetCollectionIndex();
//    std::cout << "[JetPreSelection] Number of input jets for preselection: " << jetsCol[jetCollectionIndex]->size() << std::endl;
    for( unsigned int jetIndex = 0; jetIndex < jetsCol[jetCollectionIndex]->size(); jetIndex++ )
    {
    	edm::Ptr<flashgg::Jet> jet = jetsCol[jetCollectionIndex]->ptrAt( jetIndex );
        
    	bool isJet1 = true, isJet2 = true;
        
	//if(jet->partonFlavour() != -5 && jet->partonFlavour() != 5) continue;

        if(_JetDoID[0] && !(bbggTools::isJetID(jet)))
            isJet1 = false;
        if(_JetDoID[1] && !(bbggTools::isJetID(jet)))
            isJet2 = false;
//    	if(jet->pt() < _JetPt[0])
  //          isJet1 = false;
    	if(fabs(jet->eta()) > _JetEta[0] )
            isJet1 = false;
    	if( _JetDoPUID[0]  )
           isJet1 = false;
    	if(jet->pt() < _JetPt[1])
            isJet2 = false;
    	if(fabs(jet->eta()) > _JetEta[1] )
            isJet2 = false;
//    	if( _JetDoPUID[1] )
//            isJet2 = false;
 	if(jet->bDiscriminator(_bTagType) < _JetBDiscriminant[0])
            continue;
 	if( !isJet1 && !isJet2 )
            continue;
 	if( bbggTools::DeltaR(jet->p4(), diphoCandidate->leadingPhoton()->p4()) < _JetDrPho[0] 
             || bbggTools::DeltaR(jet->p4(), diphoCandidate->subLeadingPhoton()->p4()) < _JetDrPho[0] ) continue;

 	Jets.push_back(jet);
     }
//	std::cout << "After Jet selection... 1, Jets.size(): " << Jets.size() << std::endl;

return Jets;
}




bool bbggTools::AnalysisSelection( vector<edm::Ptr<flashgg::DiPhotonCandidate>> diphoCol,
                                    JetCollectionVector jetsCol)
{
	bool cutsChecked = bbggTools::CheckCuts();
	if(!cutsChecked) {
		std::cout << "You haven't filled all the cuts correctly!" << std::endl;
		return 0;
	}

    _isSignal = 0;
    _isPhotonCR = 0;
    hasLeadJet = 0;
    hasSubJet = 0;
	
    bool isValidDiPhotonCandidate = false;
    if(DEBUG) std::cout << "Being Analysis Selection..." << std::endl;

    edm::Ptr<reco::Vertex> CandVtx;

    //Begin DiPhoton Loop/Selection -----------------------------------------------------------
    for( unsigned int diphoIndex = 0; diphoIndex < diphoCol.size(); diphoIndex++ )
    {
        
     if(DEBUG) std::cout << "[bbggTools::AnalysisSelection] Photon loop 1... " << diphoCol.size() << "\t" << diphoIndex << std::endl;
     
 	 if(_DiPhotonOnlyFirst && diphoIndex > 0 ) break;

 	 edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphoCol[ diphoIndex ];
	 	 
 	 double dipho_pt = dipho->pt();
	 double dipho_eta = dipho->eta();
	 double dipho_mass = dipho->mass();
     
     //Cuts on diphoton object
 	 if(dipho_mass < _DiPhotonMassWindow[0] || dipho_mass > _DiPhotonMassWindow[1]) continue;
 	 if(fabs(dipho_eta) > _DiPhotonEta[0] ) continue;
 	 if(dipho_pt < _DiPhotonPt[0] ) continue;
		 
     //Cuts on photons
 	 double pho1_pt = dipho->leadingPhoton()->pt();
 	 if( pho1_pt < dipho_mass*_PhotonPtOverDiPhotonMass[0] ) continue;
     
	 double pho2_pt = dipho->subLeadingPhoton()->pt();
 	 if( pho2_pt < dipho->mass()*_PhotonPtOverDiPhotonMass[1] ) continue;
     
 	 double pho1_eta = dipho->leadingPhoton()->superCluster()->eta();
 	 if( fabs(pho1_eta) > _PhotonEta[1] ) continue;
     
	 double pho2_eta = dipho->subLeadingPhoton()->superCluster()->eta();
 	 if( fabs(pho2_eta) > _PhotonEta[1] ) continue;
     
     if(DEBUG) std::cout << "[bbggTools::AnalysisSelection] Photon loop 2..." << std::endl;
     
     int pho_elv[2];
 	 pho_elv[0] = dipho->leadingPhoton()->passElectronVeto();
	 pho_elv[1] = dipho->subLeadingPhoton()->passElectronVeto();

	 int pho_ids[2];
     pho_ids[0] = 1;
     pho_ids[1] = 1;
     if(DEBUG) std::cout << "[bbggTools::AnalysisSelection] Photon loop 3..." << std::endl; 
     for( int whichPho = 0; whichPho < 2; whichPho++)
     {
	 if(DEBUG) cout << "WHICH PHO: " << whichPho << " _PhotonDoElectronVeto[whichPho] " << _PhotonDoElectronVeto[whichPho] << endl; 
         if( _PhotonDoElectronVeto[whichPho] ) {
		pho_ids[whichPho] = pho_elv[whichPho];
		if(DEBUG) std::cout << "DOING VETO! pass: " << pho_elv[whichPho] << " " << whichPho << " "<< dipho->leadingPhoton()->passElectronVeto() << " "<< dipho->subLeadingPhoton()->passElectronVeto() << endl; 
	 }
         int pho1Index = 0; //Index 0 = barrel, 1 = endcap
         double pho_eta = (whichPho) ? fabs(pho2_eta) : fabs(pho1_eta);
         if( pho_eta > _PhotonEta[0] ) pho1Index = 1;	 
         if( _PhotonDoID[whichPho] )
         { 
             std::map<int, vector<double> > theIDWP = bbggTools::getWhichID(_phoWhichID[whichPho]);
             if(theIDWP.size() < 1) {
                 std::cout << "[bbggTools::AnalysisSelection] _phoWhichID[ " << whichPho << "] = " << _phoWhichID[whichPho] << " has to be either loose, medium or tight!" << std::endl;
                 return 0;
             }
             int pho1_id = (whichPho) ? bbggTools::isPhoID(dipho->subLeadingPhoton(), theIDWP[pho1Index]) : bbggTools::isPhoID(dipho->leadingPhoton(), theIDWP[pho1Index]) ;
             if (!pho1_id) pho_ids[whichPho] = 0;
             
         }
         if(_PhotonDoISO[whichPho])
         {
             std::map<int, vector<double> > theISOWP = bbggTools::getWhichISO(_phoWhichISO[whichPho]);
             if(theISOWP.size() < 1) {
                 std::cout << "[bbggTools::AnalysisSelection] _phoWhichISO[ " << whichPho << "] = " << _phoWhichISO[whichPho] << " has to be either loose, medium or tight!" << std::endl;
                 return 0;
             }
             int pho1_id = bbggTools::isPhoISO(dipho, whichPho, theISOWP[pho1Index], _nhCorr[pho1Index], _phCorr[pho1Index]);
             if (!pho1_id) pho_ids[whichPho] = 0;
         }//here
     }
     if(DEBUG) std::cout << "[bbggTools::AnalysisSelection] After Photon loop..." << std::endl;

     if(pho_ids[0] == true && pho_ids[1] == true){
 	  isValidDiPhotonCandidate = true;
   	  CandVtx = dipho->vtx();
 	  diphoCandidate = dipho;
	  hasDiPho = true;
 	  _isSignal = 1;
	  _isPhotonCR = 0;
	  break;
     }
     if(_isPhotonCR) continue;
     if(pho_ids[0] == true && pho_ids[1] == false){
 	  isValidDiPhotonCandidate = true;
   	  CandVtx = dipho->vtx();
 	  diphoCandidate = dipho;
	  hasDiPho = true;
 	  _isSignal = 0;
	  _isPhotonCR = 1;
     }
     if(pho_ids[0] == false && pho_ids[1] == true){
 	  isValidDiPhotonCandidate = true;
   	  CandVtx = dipho->vtx();
 	  diphoCandidate = dipho;
	  hasDiPho = true;
 	  _isSignal = 0;
	  _isPhotonCR = 1;
     }
     

    }
    if( isValidDiPhotonCandidate == false ) return 0;
    if(DEBUG) std::cout << "Passed diphoton selection..." << std::endl;
    //End DiPhoton Loop/Selection -----------------------------------------------------------

    //Begin Jets Loop/Selection ------------------------------------------------------------
    std::vector<edm::Ptr<flashgg::Jet>> Jets;
    if(DEBUG) std::cout << "Begin Jet selection..." << std::endl;
    unsigned int jetCollectionIndex = diphoCandidate->jetCollectionIndex();
    for( unsigned int jetIndex = 0; jetIndex < jetsCol[jetCollectionIndex]->size(); jetIndex++ )
    {
    	edm::Ptr<flashgg::Jet> jet = jetsCol[jetCollectionIndex]->ptrAt( jetIndex );
        
    	bool isJet1 = true, isJet2 = true;
        
	//if(jet->partonFlavour() != -5 && jet->partonFlavour() != 5) continue;

        if(_JetDoID[0] && !(bbggTools::isJetID(jet)))
            isJet1 = false;
        if(_JetDoID[1] && !(bbggTools::isJetID(jet)))
            isJet2 = false;
    	if(jet->pt() < _JetPt[0])
            isJet1 = false;
    	if(fabs(jet->eta()) > _JetEta[0] )
            isJet1 = false;
    	if( _JetDoPUID[0] && jet->passesPuJetId(CandVtx) == 0 )
           isJet1 = false;
    	if(jet->pt() < _JetPt[1])
            isJet2 = false;
    	if(fabs(jet->eta()) > _JetEta[1] )
            isJet2 = false;
    	if( _JetDoPUID[1] && jet->passesPuJetId(CandVtx) == 0 )
            isJet2 = false;
 		if(jet->bDiscriminator(_bTagType) < _JetBDiscriminant[0])
            continue;
 		if( !isJet1 && !isJet2 )
            continue;
 		if( bbggTools::DeltaR(jet->p4(), diphoCandidate->leadingPhoton()->p4()) < _JetDrPho[0] 
             || bbggTools::DeltaR(jet->p4(), diphoCandidate->subLeadingPhoton()->p4()) < _JetDrPho[0] ) continue;
 		Jets.push_back(jet);
	}
	if(DEBUG) std::cout << "After Jet selection... 1, Jets.size(): " << Jets.size() << std::endl;

    if(Jets.size() < 2 ) return 0;
	
	if(DEBUG) std::cout << "After Jet selection... 2" << std::endl;

    edm::Ptr<flashgg::Jet> jet1, jet2;
    bbggTools::LorentzVector DiJet(0,0,0,0);
//    double dijetPt_ref = 0;
    double sumbtag_ref = 0;
    double sumpt_ref = 0;
    bool hasDiJet = false;
    sumbtag_ref = sumbtag_ref;
    sumpt_ref = sumpt_ref;

	if(DEBUG) std::cout << "Jet sorting... " << std::endl;
    for(unsigned int iJet = 0; iJet < Jets.size(); iJet++)
    {
 		unsigned int isbjet = 0;
 		if( Jets[iJet]->bDiscriminator(_bTagType) > _JetBDiscriminant[1] ) isbjet = 1;
 		for(unsigned int jJet = iJet+1; jJet < Jets.size(); jJet++)
 		{
 	  	 	unsigned int isbjet2 = 0;
 	  		if( Jets[jJet]->bDiscriminator(_bTagType) > _JetBDiscriminant[1] ) isbjet2 = 1;
 	  		unsigned int totalbjet = isbjet + isbjet2;
 	  		if(_n_bJets && totalbjet != _n_bJets) continue;

 	  		bbggTools::LorentzVector dijet = Jets[iJet]->p4() + Jets[jJet]->p4();
	  		if(dijet.mass() < _DiJetMassWindow[0] || dijet.mass() > _DiJetMassWindow[1]) continue;
			if( bbggTools::DeltaR( dijet, diphoCandidate->p4() ) < _CandidatesDeltaR[0] ) continue;
						
			double sumbtag = Jets[iJet]->bDiscriminator(_bTagType) + Jets[jJet]->bDiscriminator(_bTagType);

			//TESTCUT
//			double prodbtag = (Jets[iJet]->bDiscriminator(_bTagType))*(Jets[jJet]->bDiscriminator(_bTagType));
//			if(prodbtag < 0.3)
//				continue;
			if( bbggTools::DeltaR(Jets[iJet]->p4(), Jets[jJet]->p4()) < 0.5)
				continue;

			double sumpt = Jets[iJet]->pt() + Jets[jJet]->pt();
// 	  		if(dijet.pt() > dijetPt_ref && dijet.pt() > _DiJetPt[0] && fabs(dijet.Eta()) < _DiJetEta[0] )
 	  		if(sumbtag > sumbtag_ref && dijet.pt() > _DiJetPt[0] && fabs(dijet.Eta()) < _DiJetEta[0] )
// 	  		if(sumpt > sumpt_ref && dijet.pt() > _DiJetPt[0] && fabs(dijet.Eta()) < _DiJetEta[0] )
 	  		{
				hasDiJet = true;
//				dijetPt_ref = dijet.pt();
				sumbtag_ref = sumbtag;
				sumpt_ref = sumpt;
				if( Jets[iJet]->pt() > Jets[jJet]->pt() ) {
					jet1 = Jets.at(iJet);
					jet2 = Jets.at(jJet);
				} else {
					jet2 = Jets.at(iJet);
					jet1 = Jets.at(jJet);
				} 
			}
		}
	}

    if( hasDiJet == false ) return 0;
   
	leadingJetCandidate = jet1;
	hasLeadJet = true;
	subleadingJetCandidate = jet2;
	hasSubJet = true;

	if(DEBUG) std::cout << "Passed analysis selection!..." << std::endl;
	return 1;
}
