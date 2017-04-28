import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/group/phys_bphys/asanchez/MINIAOD/120EB174-5D92-E611-BD71-001E67E6F931.root')
)

process.demo1 = cms.EDAnalyzer("MiniAODTriggerAnalyzer",
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
)

process.demo2 = cms.EDAnalyzer("PackedCandAnalyzer",
    electrons = cms.InputTag("slimmedElectrons"),
    muons = cms.InputTag("miniaodPATMuonsWithTrigger"), #slimmedMuons"),
    jets = cms.InputTag("slimmedJets"),
    pfCands = cms.InputTag("packedPFCandidates"),
)

process.load("myAnalyzers.JPsiKsPAT.miniaodTriggerMatcher_cfi")

process.p = cms.Path(process.demo1*process.miniaodPATMuonsWithTriggerSequence*process.demo2)
