import sys
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = '231013_JetMET0_Run2023C-PromptReco-v3'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = 'ConfFile_submit.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['out_hist.root']
config.JobType.maxJobRuntimeMin = 600
config.JobType.maxMemoryMB = 2500


config.Data.inputDataset = '/JetMET0/Run2023C-PromptReco-v3/MINIAOD'
#/Muon1/Run2023C-PromptReco-v4/MINIAOD

config.Data.unitsPerJob = 10
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/L1Trigger/ekauffma'

config.JobType.lumiMask = '/eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_367790_Golden.json'

config.Site.storageSite = 'T2_CH_CERN'
config.Site.blacklist = ['T2_US_MIT'] # jobs are failing here often, so blacklisting
