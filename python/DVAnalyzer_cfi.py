import FWCore.ParameterSet.Config as cms

DVAnalyzer = cms.EDAnalyzer('DVAnalyzer',
                      jetTag = cms.untracked.InputTag("slimmedJets"),
                      trackRefTag = cms.untracked.InputTag("prod"),
                      trackRefMapTag = cms.untracked.InputTag("prod"),
                      trackRescaleMapTag = cms.untracked.InputTag("rescale"),
                      beamspotTag = cms.untracked.InputTag("offlineBeamSpot"),
                      #tracksTag = cms.untracked.InputTag("packedPFCandidates"),
                      muonTag = cms.untracked.InputTag("slimmedMuons"),
                      vertexTag = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
                      trackAssoTag = cms.untracked.InputTag("prod"),
                      triggerTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
                      processName = cms.untracked.string("HLT"),
                      triggerNameTag = cms.untracked.string("HLT_PFHT1050_v18"),
                      debug = cms.untracked.bool(False),
                    )
