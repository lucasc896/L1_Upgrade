import FWCore.ParameterSet.Config as cms

from UserCode.L1TriggerDPG.l1Ntuple_cfg import *

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cfi')
from Configuration.AlCa.autoCond import autoCond

#process.load("CondCore.DBCommon.CondDBSetup_cfi")
#process.load("FWCore.MessageLogger.MessageLogger_cfi")

# global tag
process.GlobalTag.globaltag = 'GR_R_52_V7::All'
#process.GlobalTag.globaltag = autoCond['hltonline'].split(',')[0]

####################################################################################

## Trigger bit requirement to use only zero bias triggered events
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt
process.ZeroBiasAve = hlt.triggerResultsFilter.clone()
process.ZeroBiasAve.triggerConditions = cms.vstring('HLT_ZeroBias*',)
process.ZeroBiasAve.hltResults = cms.InputTag( "TriggerResults", "", "HLT" )
process.ZeroBiasAve.l1tResults = cms.InputTag("")
process.ZeroBiasAve.throw = cms.bool( False )

####################################################################################

#process.skimL1Algos = cms.Path(process.ZeroBiasAve)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# get emulator modules
process.load('Configuration.StandardSequences.L1HwVal_cff')

# write out internal jets
process.valGctDigis.writeInternalData = cms.bool(True)

# convert internal jet digis to L1Extra
process.load('EventFilter.GctRawToDigi.gctInternJetProd_cfi')

# clone L1Extra ntuple branches
import L1Trigger.L1ExtraFromDigis.l1extraParticles_cfi
process.l1EmulatorExtraTree = L1Trigger.L1ExtraFromDigis.l1extraParticles_cfi.l1extraParticles.clone()
process.l1EmulatorExtraTree.centralJetSource = cms.InputTag("gctInternJetProd", "Internal")

#valL1extraParticles.muonSource = cms.InputTag("valGmtDigis")
#valL1extraParticles.nonIsolatedEmSource = cms.InputTag("valGctDigis","nonIsoEm")
#valL1extraParticles.isolatedEmSource = cms.InputTag("valGctDigis","isoEm")
#valL1extraParticles.forwardJetSource = cms.InputTag("valGctDigis","forJets")
#valL1extraParticles.centralJetSource = cms.InputTag("valGctDigis","cenJets")
#valL1extraParticles.tauJetSource = cms.InputTag("valGctDigis","tauJets")
#valL1extraParticles.etMissSource = cms.InputTag("valGctDigis")
#valL1extraParticles.htMissSource = cms.InputTag("valGctDigis")
#valL1extraParticles.etTotalSource = cms.InputTag("valGctDigis")
#valL1extraParticles.etHadSource = cms.InputTag("valGctDigis")
#valL1extraParticles.hfRingEtSumsSource = cms.InputTag("valGctDigis")
#valL1extraParticles.hfRingBitCountsSource = cms.InputTag("valGctDigis")

process.l1EmulatorExtraTree.maxL1Extra = cms.uint32(500)

#process.p+=(process.valGctDigis)    # GCT emulator
#process.p+=(process.gctInternJetProducer)  # convert to L1Extra
#process.p+=(process.l1EmulatorExtraTree)   # add an ntuple branch

process.p+=(process.ZeroBiasAve)

process.p.remove(process.gtDigis)
process.p.remove(process.gtEvmDigis)
process.p.remove(process.gctDigis)
process.p.remove(process.dttfDigis)
process.p.remove(process.csctfDigis)
#process.p.remove(process.l1RecoTreeProducer)
#process.p.remove(process.l1MuonRecoTreeProducer)
#process.p.remove(process.l1extraParticles)



readFiles.extend( [
#	'file:/gpfs_phys/storm/cms/data/Run2011B/L1MuHPF/AOD/PromptReco-v1/000/178/203/C8E95CE1-6EF5-E011-AAF0-BCAEC5329700.root'
#	'file:/gpfs_phys/storm/cms/data/Run2011B/ZeroBiasHPF0/RECO/PromptReco-v1/000/178/203/063F2E38-86F5-E011-A90B-001D09F2532F.root'
	'file:/storage/cl7359/trigger_studies/MinBias_2012_PeriodA.root'
] )

