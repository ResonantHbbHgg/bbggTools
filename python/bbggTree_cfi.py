import FWCore.ParameterSet.Config as cms
import flashgg.bbggTools.parameters as param

### For more information on each parameter, see parameters.py

bbggtree = cms.EDAnalyzer('bbggTree',
	DiPhotonTag=param._DiPhotonTag,
	JetTag=param._JetTag,
	rhoFixedGridCollection=param._rhoFixedGridCollection,
	PhotonPtOverDiPhotonMass=param._PhotonPtOverDiPhotonMass,
	PhotonEta=param._PhotonEta,
	PhotonHoverE=param._PhotonHoverE,
	PhotonSieie=param._PhotonSieie,
	PhotonR9=param._PhotonR9,
	PhotonChargedIso=param._PhotonChargedIso,
	PhotonNeutralIso=param._PhotonNeutralIso,
	PhotonPhotonIso=param._PhotonPhotonIso,
	PhotonElectronVeto=param._PhotonElectronVeto,
	PhotonDoID=param._PhotonDoID,
	PhotonDoISO=param._PhotonDoISO,
	DiPhotonPt=param._DiPhotonPt,
	DiPhotonEta=param._DiPhotonEta,
	DiPhotonMassWindow=param._DiPhotonMassWindow,
	DiPhotonOnlyFirst=param._DiPhotonOnlyFirst,
	JetPtOverDiJetMass=param._JetPtOverDiJetMass,
	JetEta=param._JetEta,
	JetBDiscriminant=param._JetBDiscriminant,
	JetDoPUID=param._JetDoPUID,
	n_bJets=param._n_bJets,
	DiJetPt=param._DiJetPt,
	DiJetEta=param._DiJetEta,
	DiJetMassWindow=param._DiJetMassWindow,
	CandidateMassWindow=param._CandidateMassWindow,
	CandidatePt=param._CandidatePt,
	CandidateEta=param._CandidateEta,
	bTagType=param._bTagType
)
