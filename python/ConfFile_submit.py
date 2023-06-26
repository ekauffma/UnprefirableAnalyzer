import FWCore.ParameterSet.Config as cms

process = cms.Process("L1JetMini")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('/store/data/Run2023C/ZeroBias/MINIAOD/PromptReco-v4/000/368/423/00000/23a45b1e-da4c-442d-bcb8-df5735c99ef7.root')
                            )

process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('368382:37-368382:792')

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('out_hist.root')
                                   )

process.demo = cms.EDAnalyzer('MiniAnalyzer',

                              )


process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.p = cms.Path(process.demo)
