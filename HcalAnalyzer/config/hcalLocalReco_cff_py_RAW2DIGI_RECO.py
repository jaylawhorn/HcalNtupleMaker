# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/python/hcalLocalReco_cff.py -s RAW2DIGI,RECO --conditions=GR_E_V49::All --filein=foo.root --fileout=bar.root --no_exec --data
import FWCore.ParameterSet.Config as cms

process = cms.Process('RERECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/data/Run2016B/HLTPhysics/RAW/v1/000/272/762/00000/688B49D6-E613-E611-9849-02163E012611.root'),
#    fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/s/sabrandt/public/Samples/pickevents_RAW_merged_50ns.root'),
#    fileNames = cms.untracked.vstring('/store/data/Run2015E/HighPtJet80/RAW/v1/000/261/395/00000/B69F2D67-5F8D-E511-9E1F-02163E013904.root'),
)

process.options = cms.untracked.PSet(
)

from Configuration.StandardSequences.RawToDigi_cff import *
process.RawToDigi = cms.Sequence(hcalDigis)

from RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi import *
process.hcalOOTPileupESProducer = cms.ESProducer('OOTPileupDBCompatibilityESProducer')

from RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_hbhe_cfi import *
#from RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_ho_cfi import *
#from RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_hf_cfi import *
#from RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_zdc_cfi import *

hbheprereco.setNegativeFlags=False

process.hcalLocalRecoSequence = cms.Sequence(hbheprereco)


#import RecoLocalCalo.HcalRecProducers.HcalHitReconstructorHLT_hbhe_cfi

#process.load('RecoLocalCalo.HcalRecProducers.HcalHitReconstructorHLT_hbhe_cfi')
#process.hcalLocalRecoSequenceHLT = cms.Sequence(process.hbheprerecoHLT)

#process.load('RecoLocalCalo.HcalRecProducers.HcalHitReconstructorSinglePulse_hbhe_cfi')
#process.hcalLocalRecoSequenceSingle = cms.Sequence(process.hbheprerecoSingle)

#process.hbheprereco.puCorrMethod = cms.int32(3)

from RecoLocalCalo.HcalRecProducers.HBHEIsolatedNoiseReflagger_cfi import *
process.hcalGlobalRecoSequence = cms.Sequence(hbhereco)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Configuration/python/hcalLocalReco_cff.py nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string("HCALTree.root")
)
process.ExportTree = cms.EDAnalyzer("HcalAnalyzer",
  hbheInput = cms.InputTag('hbheprereco'), 
  IsData = cms.untracked.bool(True),
  TriggerResults = cms.InputTag('TriggerResults','','HLT'),
  TotalChargeThreshold = cms.untracked.double(-9999)
)


# Output definition
#process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
#    dataset = cms.untracked.PSet(
#        dataTier = cms.untracked.string(''),
#        filterName = cms.untracked.string('')
#    ),
#    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
#    fileName = cms.untracked.string('file:foobar_254608.root'),
##    outputCommands = cms.untracked.vstring( #process.RECOSIMEventContent.outputCommands,
##        keep 
##),
#    splitLevel = cms.untracked.int32(0)
#)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v8', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
#process.reco_step1_ecal = cms.Path(process.ecalLocalRecoSequence)
#process.reco_step1 = cms.Path(process.hcalLocalRecoSequence+process.hcalLocalRecoSequenceHLT+process.hcalLocalRecoSequenceSingle)
#process.reco_step1 = cms.Path(process.hcalLocalRecoSequence*process.hcalGlobalRecoSequence)
process.reco_step1 = cms.Path(process.hcalLocalRecoSequence)
process.reco_step2 = cms.Path(process.ExportTree)
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definitio
process.schedule = cms.Schedule(process.raw2digi_step,process.reco_step1,process.reco_step2,process.endjob_step)
#process.schedule = cms.Schedule(process.raw2digi_step,process.reco_step1_ecal,process.reco_step1,process.reco_step2,process.endjob_step)
#process.schedule = cms.Schedule(process.raw2digi_step,process.reco_step1,process.reco_step2,process.endjob_step,process.RECOSIMoutput_step)
#process.schedule = cms.Schedule(process.reco_step1,process.reco_step2,process.endjob_step,process.RECOSIMoutput_step)


