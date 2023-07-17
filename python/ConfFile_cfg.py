import FWCore.ParameterSet.Config as cms

process = cms.Process("L1JetMini")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('/store/data/Run2023C/Muon1/MINIAOD/PromptReco-v4/000/367/770/00000/17020b67-65fa-49cd-8e07-5ff930d563cd.root')
                            )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('testhisto.root')
                                   )

process.demo = cms.EDAnalyzer('UnprefirableAnalyzer',

                              )


process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.p = cms.Path(process.demo)
