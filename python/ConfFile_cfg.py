import FWCore.ParameterSet.Config as cms

process = cms.Process("L1JetMini")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('file:/eos/cms/store/data/Run2023C/MinimumBias/MINIAOD/PromptReco-v1/000/367/513/00000/d0f64962-750a-428b-8b9b-b3d9efc8e8de.root')
                            )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('testhisto.root')
                                   )

process.demo = cms.EDAnalyzer('MiniAnalyzer',

                              )

process.p = cms.Path(process.demo)
