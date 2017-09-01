import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	    #'file:miniaod-00379.root'
	    '/store/data/Run2017B/MuOnia/MINIAOD/PromptReco-v1/000/297/723/00000/9040368C-DE5E-E711-ACFF-02163E0134FF.root'
    #'/store/group/phys_bphys/asanchez/MINIAOD/120EB174-5D92-E611-BD71-001E67E6F931.root'
    )
)

process.upkPatTrigger = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
    patTriggerObjectsStandAlone = cms.InputTag( 'slimmedPatTrigger' ),
    triggerResults              = cms.InputTag( 'TriggerResults::HLT' ),
    unpackFilterLabels          = cms.bool( True )
)

process.demo1 = cms.EDAnalyzer("MiniAODTriggerAnalyzer",
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("unpackedPatTrigger")#"slimmedPatTrigger")#selectedPatTrigger"),
)

process.demo2 = cms.EDAnalyzer("PackedCandAnalyzer",
    electrons = cms.InputTag("slimmedElectrons"),
    muons = cms.InputTag("slimmedMuonsWithTrigger"), #slimmedMuons"),
    jets = cms.InputTag("slimmedJets"),
    pfCands = cms.InputTag("packedPFCandidates"),
)

process.load("Rootupler.Generation.slimmedMuonsTriggerMatcher_cfi")

#process.p = cms.Path(process.slimmedMuonsWithTriggerSequence*process.demo1*process.slimmedMuonsWithTriggerSequence*process.demo2)
process.p = cms.Path(process.slimmedMuonsWithTriggerSequence*process.demo1*process.demo2)
