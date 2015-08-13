#include "flashgg/bbggTools/interface/bbggTools.h"
//FLASHgg files
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"


using namespace std;

bool DEBUG = false;

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
	double extraFactor = exp(nhCorr[0]*pho->pt()+nhCorr[1]);
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
	const flashgg::Photon* pho = (whichPho) ? dipho->subLeadingPhoton() : dipho->leadingPhoton();
	nhiso = bbggTools::getNHisoToCutValue( pho );
	phiso = bbggTools::getPHisoToCutValue( pho );
	/*
  	if(whichPho == 0) nhiso = bbggTools::getNHisoToCutValue( dipho->leadingPhoton() );
  	if(whichPho == 1) nhiso = bbggTools::getNHisoToCutValue( dipho->subLeadingPhoton() );
  	if(whichPho == 0) phiso = bbggTools::getPHisoToCutValue( dipho->leadingPhoton() );
	if(whichPho == 1) phiso = bbggTools::getPHisoToCutValue( dipho->subLeadingPhoton() );
	*/
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
	/*
  	if(whichPho == 0) nhiso = bbggTools::getNHisoToCutValue( dipho->leadingPhoton() );
  	if(whichPho == 1) nhiso = bbggTools::getNHisoToCutValue( dipho->subLeadingPhoton() );
  	if(whichPho == 0) phiso = bbggTools::getPHisoToCutValue( dipho->leadingPhoton() );
	if(whichPho == 1) phiso = bbggTools::getPHisoToCutValue( dipho->subLeadingPhoton() );
	*/	
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
		_PhotonHoverE.size() == 0 ||
		_PhotonSieie.size() == 0 ||
		_PhotonR9.size() == 0 ||
		_PhotonChargedIso.size() == 0 ||
		_PhotonNeutralIso.size() == 0 ||
		_PhotonPhotonIso.size() == 0 ||
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
		_CandidatesDeltaR.size() == 0 ) {
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

bool bbggTools::AnalysisSelection( edm::Handle<edm::View<flashgg::DiPhotonCandidate> > diphoCol, 
						edm::Handle<edm::View<flashgg::Jet> > jetsCol )
{
    bool isValidDiPhotonCandidate = false;
    if(DEBUG) std::cout << "Being AnalysisSelection..." << std::endl;

    edm::Ptr<reco::Vertex> CandVtx;

    //Begin DiPhoton Loop/Selection -----------------------------------------------------------
    for( unsigned int diphoIndex = 0; diphoIndex < diphoCol->size(); diphoIndex++ )
    {
 	 if(_DiPhotonOnlyFirst && diphoIndex > 0 ) break;

 	 edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphoCol->ptrAt( diphoIndex );
	 
 	 double dipho_pt = dipho->pt();
	 double dipho_eta = dipho->eta();
	 double dipho_mass = dipho->mass();
		 
 	 if(dipho_mass < _DiPhotonMassWindow[0] || dipho_mass > _DiPhotonMassWindow[1]) continue;
 	 if(fabs(dipho_eta) > _DiPhotonEta[0] ) continue;
 	 if(dipho_pt < _DiPhotonPt[0] ) continue;
		 
 	 double pho1_pt = dipho->leadingPhoton()->pt();
	 double pho2_pt = dipho->subLeadingPhoton()->pt();
 	 double pho1_eta = dipho->leadingPhoton()->superCluster()->eta();
	 double pho2_eta = dipho->subLeadingPhoton()->superCluster()->eta();
 	 double pho1_hoe = dipho->leadingPhoton()->hadronicOverEm();
	 double pho2_hoe = dipho->subLeadingPhoton()->hadronicOverEm();
 	 double pho1_sieie = dipho->leadingPhoton()->full5x5_sigmaIetaIeta();
	 double pho2_sieie = dipho->subLeadingPhoton()->full5x5_sigmaIetaIeta();
 	 int pho1_elveto = dipho->leadingPhoton()->passElectronVeto();
	 int pho2_elveto = dipho->subLeadingPhoton()->passElectronVeto();

 	 double pho1_chiso = bbggTools::getCHisoToCutValue( dipho, 0);
 	 double pho2_chiso = bbggTools::getCHisoToCutValue( dipho, 1);
 	 double pho1_nhiso = bbggTools::getNHisoToCutValue( dipho->leadingPhoton() );
 	 double pho2_nhiso = bbggTools::getNHisoToCutValue( dipho->subLeadingPhoton() );
 	 double pho1_phiso = bbggTools::getPHisoToCutValue( dipho->leadingPhoton() );
 	 double pho2_phiso = bbggTools::getPHisoToCutValue( dipho->subLeadingPhoton() ); 
	 		 
 	 if( pho1_pt < dipho_mass*_PhotonPtOverDiPhotonMass[0] ) continue;
 	 if( fabs(pho1_eta) > _PhotonEta[1] ) continue;
 	 if( pho2_pt < dipho->mass()*_PhotonPtOverDiPhotonMass[1] ) continue;
 	 if( fabs(pho2_eta) > _PhotonEta[1] ) continue;
		 
 	 bool pho1_id = true, pho2_id = true;
 	 if( _PhotonDoID[0] )
 	 {
 	   int pho1Index = 0;
 	   if( fabs(pho1_eta) > _PhotonEta[0] ) pho1Index = 1;
			 
 	   if( pho1_hoe > _PhotonHoverE[pho1Index] ) 	pho1_id = false;
 	   if( pho1_sieie > _PhotonSieie[pho1Index] )	 pho1_id = false;
 	   if( pho1_elveto != _PhotonElectronVeto[0] )    pho1_id = false;
 	}
 	if(_PhotonDoISO[0])
 	{
 	   int pho1Index = 0;
            if( fabs(pho1_eta) > _PhotonEta[0] ) pho1Index = 1;

 	   if( pho1_chiso > _PhotonChargedIso[pho1Index] ) pho1_id = false;
 	   if( pho1_nhiso > _PhotonNeutralIso[pho1Index] ) pho1_id = false;
 	   if( pho1_phiso > _PhotonPhotonIso[pho1Index] ) pho1_id = false;
 	}
 	if( _PhotonDoID[1] )
 	{
 	   int pho2Index = 2;
 	   if( fabs(pho2_eta) > _PhotonEta[0] ) pho2Index = 3;

            if( pho2_hoe > _PhotonHoverE[pho2Index] )   pho2_id = false;
            if( pho2_sieie > _PhotonSieie[pho2Index] ) pho2_id = false;
 	   if( pho2_elveto != _PhotonElectronVeto[1] )    pho2_id = false;
 	}
 	if(_PhotonDoISO[1])
 	{
            int pho2Index = 2;
            if( fabs(pho2_eta) > _PhotonEta[0] ) pho2Index = 3;	 

 	   if( pho2_chiso > _PhotonChargedIso[pho2Index] ) pho2_id = false;
 	   if( pho2_nhiso > _PhotonNeutralIso[pho2Index] ) pho2_id = false;
 	   if( pho2_phiso > _PhotonPhotonIso[pho2Index] ) pho2_id = false;
 	}

 	if(pho1_id == true && pho2_id == true){
 	  isValidDiPhotonCandidate = true;
   	  CandVtx = dipho->vtx();
 	  diphoCandidate = dipho;
	  hasDiPho = true;
 	  break;
 	}

    }
    if( isValidDiPhotonCandidate == false ) return 0;
    if(DEBUG) std::cout << "Passed diphoton selection..." << std::endl;
    //End DiPhoton Loop/Selection -----------------------------------------------------------

    //Begin Jets Loop/Selection ------------------------------------------------------------
    std::vector<edm::Ptr<flashgg::Jet>> Jets;
    int nJet1 = 0, nJet2 = 0;
    if(DEBUG) std::cout << "Begin Jet selection..." << std::endl;
    for( unsigned int jetIndex = 0; jetIndex < jetsCol->size(); jetIndex++ )
    {
    	edm::Ptr<flashgg::Jet> jet = jetsCol->ptrAt( jetIndex );
    	bool isJet1 = true, isJet2 = true;
    	if(jet->pt() < _JetPt[0]) isJet1 = false;
    	if(fabs(jet->eta()) > _JetEta[0] ) isJet1 = false;
    	if( _JetDoPUID[0] && jet->passesPuJetId(CandVtx) == 0 ) isJet1 = false;
    	if(jet->pt() < _JetPt[1]) isJet2 = false;
    	if(fabs(jet->eta()) > _JetEta[1] ) isJet2 = false;
    	if( _JetDoPUID[1] && jet->passesPuJetId(CandVtx) == 0 ) isJet2 = false;
    	if(isJet1) nJet1++;
    	if(isJet1 == false && isJet2) nJet2++;
 		if(jet->bDiscriminator(_bTagType) < _JetBDiscriminant[0]) continue;
 		if( !isJet1 && !isJet2 ) continue;
 		if( bbggTools::DeltaR(jet->p4(), diphoCandidate->leadingPhoton()->p4()) < _JetDrPho[0] 
             || bbggTools::DeltaR(jet->p4(), diphoCandidate->subLeadingPhoton()->p4()) < _JetDrPho[0] ) continue;
 		Jets.push_back(jet);
	}


    if(Jets.size() < 2 ) return 0;

    edm::Ptr<flashgg::Jet> jet1, jet2;
    bbggTools::LorentzVector DiJet(0,0,0,0);
//    double dijetPt_ref = 0;
//	double sumbtag_ref = 0;
	double sumpt_ref = 0;
    bool hasDiJet = false;

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
 //	  if(totalbjet != n_bJets) continue;

 	  		bbggTools::LorentzVector dijet = Jets[iJet]->p4() + Jets[jJet]->p4();
	  		if(dijet.mass() < _DiJetMassWindow[0] || dijet.mass() > _DiJetMassWindow[1]) continue;
			
			if( bbggTools::DeltaR( dijet, diphoCandidate->p4() ) < _CandidatesDeltaR[0] ) continue;
			
//			if( Jets[iJet]->pt() < 0.3*dijet.mass() + 15 ) continue;
//			if( Jets[jJet]->pt() < 0.3*dijet.mass() + 15) continue;
//			if( Jets[iJet]->pt() <  25. + 0.01*exp( 0.05*dijet.mass()) ) continue;
//			if( Jets[jJet]->pt() <  25. + 0.01*exp( 0.05*dijet.mass()) ) continue;
//			if( Jets[jJet]->pt() < 10) continue;
						
//			double sumbtag = Jets[iJet]->bDiscriminator(_bTagType) + Jets[jJet]->bDiscriminator(_bTagType);
			double sumpt = Jets[iJet]->pt() + Jets[jJet]->pt();

// 	  		if(dijet.pt() > dijetPt_ref && dijet.pt() > _DiJetPt[0] && fabs(dijet.Eta()) < _DiJetEta[0] )
// 	  		if(sumbtag > sumbtag_ref && dijet.pt() > _DiJetPt[0] && fabs(dijet.Eta()) < _DiJetEta[0] )
 	  		if(sumpt > sumpt_ref && dijet.pt() > _DiJetPt[0] && fabs(dijet.Eta()) < _DiJetEta[0] )
 	  		{
				hasDiJet = true;
//				dijetPt_ref = dijet.pt();
//				sumbtag_ref = sumbtag;
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

	return 1;
}
