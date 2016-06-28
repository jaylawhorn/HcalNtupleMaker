from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'FlatQCD_v3'
config.General.workArea = 'hcalFlatQCD2016'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'mc_hcalLocalReco_cff_py_RECOBEFMIX_DIGI_L1_DIGI2RAW_RECO.py'
config.JobType.outputFiles = ['HCALTree.root']
config.section_("Data")
#config.Data.inputDataset = '/HighPtJet80/Run2015E-v1/RAW'
config.Data.inputDataset = '/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISpring16FSPremix-ForSUSJECs_80X_mcRun2_asymptotic_v12-v1/GEN-SIM'
#config.Data.runRange = '258656'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = 200
config.Data.ignoreLocality = True
config.Data.outLFNDirBase = '/store/group/phys_susy/razor/HCALDPG/HcalNtuple' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag = 'Samples'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
