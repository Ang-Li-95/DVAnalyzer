import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
                                  #'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/StopStopbarTo2Dbar2D_M-800_CTau-1mm_TuneCP2_13TeV-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/623DC9AB-7C01-E911-BFEC-246E96D14D4C.root'
                                  #'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/GluinoGluinoToNeutralinoNeutralinoTo2T2B2S_M-800_CTau-1mm_TuneCP2_13TeV_2018-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/20000/C1239AD3-101D-B040-9466-35C20A621626.root'
                                  #'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/GluinoGluinoToNeutralinoNeutralinoTo2T2B2S_M-800_CTau-1mm_TuneCP2_13TeV_2018-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/80000/E7677EE5-376B-6F4B-B504-BC820873679A.root'
                                  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/80000/00A7C4D5-8881-5D47-8E1F-FADDC4B6FA96.root'
                )
                            )

process.prod = cms.EDProducer('trackRefProd',
                              packed_candidates_src = cms.InputTag('packedPFCandidates'),
                              pvTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              skip_weirdos = cms.bool(True)
)

process.rescale = cms.EDProducer('RescaleTracks',
                                 trackTag = cms.InputTag('prod'),
                                 doRescale = cms.bool(True)
)

#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('myOutputFile.root')
#    #,outputCommands = cms.untracked.vstring('drop *',
#    #  "keep *_generalTracks_*_*",
#    #  "keep *_globalMuons_*_*",
#    #  "keep *_MuonTrackPoints_*_*",
#    #  "keep *_TrackTrackPoints_*_*")
#
#)

process.load('LLPAnalyzer.DVAnalyzer.DVAnalyzer_cfi')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v19'
#process.demo = cms.EDAnalyzer('DemoAnalyzer'
#                              )

#process.p = cms.Path(process.demo)
DVAnalyzer = process.DVAnalyzer

process.TFileService = cms.Service("TFileService", fileName = cms.string("DV_sig_asso.root") )

process.runseq = cms.Sequence()
process.runseq += process.prod
process.runseq += process.rescale
process.runseq += DVAnalyzer
process.path = cms.Path(process.runseq)
process.schedule = cms.Schedule(process.path)
