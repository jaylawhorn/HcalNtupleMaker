from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'HighPtJet80_v3'
config.General.workArea = 'hcalHeavyIon2015'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'hcalLocalReco_cff_py_RAW2DIGI_RECO.py'
config.JobType.outputFiles = ['HCALTree.root']
config.section_("Data")
config.Data.inputDataset = '/HighPtJet80/Run2015E-v1/RAW'
#config.Data.runRange = '258656'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/group/phys_susy/razor/HCALDPG/HcalNtuple' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag = 'Samples'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
