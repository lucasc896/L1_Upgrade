import FWCore.ParameterSet.Config as cms

from UserCode.L1TriggerDPG.l1Ntuple_cfg import *

# global tag
process.GlobalTag.globaltag = 'GR_R_44_V14::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.l1NtupleProducer.gtEvmSource 			= cms.InputTag("none")
process.l1NtupleProducer.rctSource            	= cms.InputTag("none")
process.l1NtupleProducer.dttfSource           	= cms.InputTag("none")
process.l1NtupleProducer.csctfTrkSource       	= cms.InputTag("none")
process.l1NtupleProducer.csctfLCTSource       	= cms.InputTag("none")
process.l1NtupleProducer.csctfStatusSource    	= cms.InputTag("none")
process.l1NtupleProducer.csctfDTStubsSource   	= cms.InputTag("none")

process.p.remove(process.gtDigis)
process.p.remove(process.gtEvmDigis)
process.p.remove(process.gctDigis)
process.p.remove(process.dttfDigis)
process.p.remove(process.csctfDigis)
process.p.remove(process.l1MuonRecoTreeProducer)


readFiles.extend( [
#	'file:/gpfs_phys/storm/cms/data/Run2011B/L1MuHPF/AOD/PromptReco-v1/000/178/203/C8E95CE1-6EF5-E011-AAF0-BCAEC5329700.root'
	'file:/gpfs_phys/storm/cms/data/Run2011B/ZeroBiasHPF0/RECO/PromptReco-v1/000/178/203/063F2E38-86F5-E011-A90B-001D09F2532F.root'
#	'file:/storage/cl7359/trigger_studies/MinBias_2012_PeriodA.root'
] )

