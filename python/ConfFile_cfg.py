import FWCore.ParameterSet.Config as cms

process = cms.Process("L1Trigger")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('file:/eos/cms/store/data/Run2023C/JetMET1/MINIAOD/22Sep2023_v2-v1/50000/81edc4a1-b25b-45f0-8c2d-6ac30010a78e.root')
                            )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('testhisto.root')
                                   )

process.demo = cms.EDAnalyzer('UnprefirableAnalyzer',

                              )


process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data', '')

process.p = cms.Path(process.demo)
