from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'SingleNeutrino_v2'
config.General.workArea = 'SingleNeutrino_80X'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'mc_hcalLocalReco_cff_py_RECOBEFMIX_DIGI_L1_DIGI2RAW_RECO.py'
#config.JobType.psetName = 'hcalLocalReco_cff_py_RAW2DIGI_RECO.py'
config.JobType.outputFiles = ['HCALTree.root']
config.section_("Data")
#config.Data.inputDataset = '/ZeroBias/Run2016B-v2/RAW'
#config.Data.inputDataset = '/HighPtJet80/Run2015E-v1/RAW'
#config.Data.inputDataset = '/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISpring16FSPremix-ForSUSJECs_80X_mcRun2_asymptotic_v12-v1/GEN-SIM'
config.Data.inputDataset = '/SingleNeutrino/RunIISpring16DR80-PUSpring16_RECO_NZS_80X_mcRun2_asymptotic_2016_v3-v1/GEN-SIM-RAW'
#config.Data.runRange = '258656'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
#config.Data.totalUnits = 200
config.Data.ignoreLocality = True
config.Data.outLFNDirBase = '/store/group/dpg_hcal/comm_hcal/RecoAlgos/Summer16Method2Update/HcalNtuple' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag = 'Samples'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
