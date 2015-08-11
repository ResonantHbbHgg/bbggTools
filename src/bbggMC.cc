#include "flashgg/bbggTools/interface/bbggMC.h"
//FLASHgg files
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"


using namespace std;

//bool DEBUG = false;

double bbggMC::DeltaR(bbggMC::LorentzVector vec1, bbggMC::LorentzVector vec2)
{
	double R2 = (vec1.phi() - vec2.phi())*(vec1.phi() - vec2.phi()) + (vec1.eta() - vec2.eta())*(vec1.eta() - vec2.eta());
	return sqrt(R2);
}

edm::Ptr<flashgg::DiPhotonCandidate> bbggMC::GetSelected_diphoCandidate()
{
	if(!hasDiPho){
		std::cout << "[bbggMC::GetSelected_diphoCandidate] ERROR: diphoCandidate has not been set yet. Run bbggMC::AnalysisSelection. If you did, then the event failed the selection requirements." << std::endl;
		return diphoCandidate;
	}
	return diphoCandidate;
}

edm::Ptr<flashgg::Jet> bbggMC::GetSelected_leadingJetCandidate()
{
	if(!hasLeadJet){
		std::cout << "[bbggMC::GetSelected_leadingJetCandidate] ERROR: leadingJetCandidate has not been set yet. Run bbggMC::AnalysisSelection. If you did, then the event failed the selection requirements." << std::endl;
		return leadingJetCandidate;
	}
	return leadingJetCandidate;
}

edm::Ptr<flashgg::Jet> bbggMC::GetSelected_subleadingJetCandidate()
{
	if(!hasSubJet){
		std::cout << "[bbggMC::GetSelected_subleadingJetCandidate] ERROR: subleadingJetCandidate has not been set yet. Run bbggMC::AnalysisSelection. If you did, then the event failed the selection requirements." << std::endl;
		return subleadingJetCandidate;
	}
	return subleadingJetCandidate;
}

bool bbggMC::MatchTruth( edm::Handle<edm::View<flashgg::DiPhotonCandidate> > diphoCol, 
						edm::Handle<edm::View<flashgg::Jet> > jetsCol,
						edm::Handle<edm::View<reco::GenParticle> > genCol )
{
	vector<LorentzVector> genPhos;
	vector<LorentzVector> genJets;
	for( unsigned int iGen = 0; iGen < genCol->size(); iGen++)
	{
		edm::Ptr<reco::GenParticle> genPar = genCol->ptrAt(iGen);
		
		int id = genPar->pdgId();
		int motherId = (genPar->mother()) ? genPar->mother()->pdgId() : -1;
		int status = genPar->status();
		if( motherId == 25 )
		{
			if (id == 22 && status == 1) genPhos.push_back(genPar->p4());
			if (id == 5 || id == -5) genJets.push_back(genPar->p4());
		}
	}
	if(genPhos.size() != 2 || genJets.size() != 2) {
		std::cout << "[bbggMC::MatchTruth] size of genPhos col: " << genPhos.size() << std::endl;
		std::cout << "[bbggMC::MatchTruth] size of genJets col: " << genJets.size() << std::endl;
		return 0;	
	}
	
	bool isPhomatched = false;
	for (unsigned int diPho = 0; diPho < diphoCol->size(); diPho++)
	{
		edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphoCol->ptrAt( diPho );
		bool areMatched = true;
		LorentzVector PHO1 = dipho->leadingPhoton()->p4();
		LorentzVector PHO2 = dipho->subLeadingPhoton()->p4();
		for (unsigned int GPho = 0; GPho < genPhos.size(); GPho++){
			if(DeltaR(PHO1, genPhos[GPho]) < 0.3 && DeltaR(PHO2, genPhos[GPho]) < 0.3)
				areMatched = false;
		}
		if(!areMatched) continue;
		if(areMatched){
			diphoCandidate = dipho;
			isPhomatched = true;
			hasDiPho = true;
			break;
		}
	}
	
	if(!isPhomatched) {
		std::cout << "[bbggMC::MatchTruth] No matched diphoton :(" << std::endl;
		return 0;
	}
	
	vector<edm::Ptr<flashgg::Jet>> matchedJets;
	for (unsigned int GJet = 0; GJet < genJets.size(); GJet++)
	{
		float minDR = 10;
		unsigned int iMatched = 1000;
		for (unsigned int iJet = 0; iJet < jetsCol->size(); iJet++)
		{
			edm::Ptr<flashgg::Jet> jet = jetsCol->ptrAt( iJet );
			LorentzVector lJet = jet->p4(); //Jet->p4();
//			if(jet->partonFlavour() != 5 || jet->partonFlavour() != -5) continue;
			if(!jet->passesPuJetId(diphoCandidate->vtx())) continue;
			if( DeltaR(lJet, diphoCandidate->leadingPhoton()->p4()) < 0.4) continue;
			if( DeltaR(lJet, diphoCandidate->subLeadingPhoton()->p4()) < 0.4) continue;
			if(DeltaR(lJet, genJets[GJet]) < minDR) {
				minDR = DeltaR(lJet, genJets[GJet]);
				iMatched = iJet;
			}
		}
		if(minDR < 0.4 && iMatched < 1000){
			edm::Ptr<flashgg::Jet> Jet = jetsCol->ptrAt( iMatched );
			matchedJets.push_back(Jet);
		}
	}
	
	
	if(matchedJets.size() < 2) {
		std::cout << "[bbggMC::MatchTruth] No matched dijet :(" << std::endl;
		return 0;
	}
	
	if(matchedJets.size() > 2) {
		std::cout << "[bbggMC::MatchTruth] More than 2 jets matched, check the code! :(" << std::endl;
		return 0;
	}
	
	if(matchedJets[0]->pt() > matchedJets[1]->pt()) {
		leadingJetCandidate = matchedJets[0];
		subleadingJetCandidate = matchedJets[1];
	} else {
		leadingJetCandidate = matchedJets[1];
		subleadingJetCandidate = matchedJets[0];		
	}
	hasLeadJet = true;
	hasSubJet = true;
	return 1;
}
