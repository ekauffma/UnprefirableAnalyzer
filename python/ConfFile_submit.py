import FWCore.ParameterSet.Config as cms

process = cms.Process("L1Trigger")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('file:/eos/cms/store/data/Run2023C/JetMET0/MINIAOD/22Sep2023_v1-v1/30000/56222725-fbae-46ca-96f8-f91acf1dc8ec.root')
                            )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('out_hist.root')
                                   )

process.demo = cms.EDAnalyzer('UnprefirableAnalyzer',

                              )


process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data', '')

process.p = cms.Path(process.demo)
