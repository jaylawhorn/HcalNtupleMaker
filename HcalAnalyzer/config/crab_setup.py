from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'data'
config.General.workArea = 'hcalLoneBunchData'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'hcalLocalReco_cff_py_RAW2DIGI_RECO.py'
config.JobType.outputFiles = ['Output.root']
config.section_("Data")
config.Data.inputDataset = '/ZeroBias/Run2015D-v1/RAW'
config.Data.runRange = '258656'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 7 
config.Data.outLFNDirBase = '/store/user/jlawhorn' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag = 'Samples'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
