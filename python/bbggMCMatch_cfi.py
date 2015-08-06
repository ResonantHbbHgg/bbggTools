import FWCore.ParameterSet.Config as cms
import flashgg.bbggTools.parameters as param

#### Comments refer to the structure
#### of the vector below the comment line

bbggmcmatch = cms.EDAnalyzer('bbggMCMatch',
	DiPhotonTag=param._DiPhotonTag,
	JetTag=param._JetTag,
	GenTag=param._GenTag,
	bTagType=param._bTagType
)