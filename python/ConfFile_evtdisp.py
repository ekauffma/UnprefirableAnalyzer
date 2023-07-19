import FWCore.ParameterSet.Config as cms

process = cms.Process("L1Trigger")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('/store/data/Run2023C/Muon1/MINIAOD/PromptReco-v4/000/368/822/00000/124a3db9-23cb-485f-b700-9182df8673a3.root')
                            )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('testhisto.root')
                                   )

process.demo = cms.EDAnalyzer('PrefiringEventDisplay',

                              )


process.MessageLogger.cerr.FwkReport.reportEvery = 100000

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.p = cms.Path(process.demo)
