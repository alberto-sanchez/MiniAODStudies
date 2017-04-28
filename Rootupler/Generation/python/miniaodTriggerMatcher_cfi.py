import FWCore.ParameterSet.Config as cms

#this is our version of the patMuonsWithTrigger using MINIAOD


### ==== Then perform a match for all HLT triggers of interest
PATmuonTriggerMatchHLT = cms.EDProducer( "PATTriggerMatcherDRDPtLessByR",
    src     = cms.InputTag( "slimmedMuons" ),
    matched = cms.InputTag( "selectedPatTrigger" ),
    matchedCuts = cms.string(""),
    maxDPtRel = cms.double( 0.5 ),
    maxDeltaR = cms.double( 0.5 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( True ) #change with respect to previous tag
)

PATmuonMatchHLTL2   = PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltL2MuonCandidates")'), 
                                                   maxDeltaR = 0.3, maxDPtRel = 10.0)       #maxDeltaR Changed accordingly to Zoltan tuning. It was: 1.2
PATmuonMatchHLTL3   = PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltL3MuonCandidates")'), 
                                                   maxDeltaR = 0.1, maxDPtRel = 10.0)       #maxDeltaR Changed accordingly to Zoltan tuning. It was: 0.5
PATmuonMatchHLTL3T  = PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltGlbTrkMuonCands")'),  
                                                   maxDeltaR = 0.1, maxDPtRel = 10.0)       #maxDeltaR Changed accordingly to Zoltan tuning. It was: 0.5
PATmuonMatchHLTTkMu = PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltHighPtTkMuonCands")'),  
                                                   maxDeltaR = 0.1, maxDPtRel = 10.0)       #maxDeltaR Changed accordingly to Zoltan tuning. It was: 0.5

miniaodPATTriggerMatchers1Mu = cms.Sequence(
      PATmuonMatchHLTL2 +
      PATmuonMatchHLTL3 +
      PATmuonMatchHLTL3T +
      PATmuonMatchHLTTkMu
)

miniaodPATTriggerMatchers1MuInputTags = [
    cms.InputTag('PATmuonMatchHLTL2'),
    cms.InputTag('PATmuonMatchHLTL3'),
    cms.InputTag('PATmuonMatchHLTL3T'),
    cms.InputTag('PATmuonMatchHLTTkMu'),
]

## ==== Embed ====
miniaodPATMuonsWithTrigger = cms.EDProducer( "PATTriggerMatchMuonEmbedder",
    src     = cms.InputTag(  "slimmedMuons" ),
    matches = cms.VInputTag()
)
miniaodPATMuonsWithTrigger.matches += miniaodPATTriggerMatchers1MuInputTags

## ==== Trigger Sequence ====
miniaodPATMuonsWithTriggerSequence = cms.Sequence(
    miniaodPATTriggerMatchers1Mu *
    miniaodPATMuonsWithTrigger
)
