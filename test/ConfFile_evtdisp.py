import FWCore.ParameterSet.Config as cms

process = cms.Process("L1Trigger")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('/store/data/Run2023C/Muon1/MINIAOD/PromptReco-v4/000/368/636/00000/116fb49a-58d8-4619-9f58-3634f8104c4f.root',
    '/store/data/Run2023C/Muon1/MINIAOD/PromptReco-v4/000/368/636/00000/2c6d6d2e-2d14-4464-b546-4beb13ac0d7f.root',
    '/store/data/Run2023C/Muon1/MINIAOD/PromptReco-v4/000/368/636/00000/4d024e6e-eb65-4f86-bf2a-9c439fa84cf8.root',
    '/store/data/Run2023C/Muon1/MINIAOD/PromptReco-v4/000/368/636/00000/52bb4f8b-61f3-4789-9379-ddb8f9ea277b.root',
    '/store/data/Run2023C/Muon1/MINIAOD/PromptReco-v4/000/368/636/00000/601ee259-1393-4915-9c19-8c35f3af7320.root',
    '/store/data/Run2023C/Muon1/MINIAOD/PromptReco-v4/000/368/636/00000/8a89af9f-14d2-4b03-9f8a-a0d792a37945.root',
    '/store/data/Run2023C/Muon1/MINIAOD/PromptReco-v4/000/368/636/00000/ab4e3f47-4b07-48bd-aba0-6bcb2c8c5041.root',
    '/store/data/Run2023C/Muon1/MINIAOD/PromptReco-v4/000/368/636/00000/b1c75c5f-7a62-4b8c-9655-95848110f4ab.root',
    '/store/data/Run2023C/Muon1/MINIAOD/PromptReco-v4/000/368/636/00000/cb4292e1-abc0-4efe-b8b6-10dc00248b53.root',
    '/store/data/Run2023C/Muon1/MINIAOD/PromptReco-v4/000/368/636/00000/ccf7b0f0-0c92-42bc-bcd2-958d40813fac.root')
                            )

process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('368636:798')

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
