import sys
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Muon1_2023C_manyruns_recoInfo_v4'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = 'ConfFile_submit.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['out_hist.root']

config.Data.inputDataset = '/Muon1/Run2023C-PromptReco-v4/MINIAOD'
config.Data.unitsPerJob = 10
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/L1Trigger/ekauffma'

config.Site.storageSite = 'T2_CH_CERN'
