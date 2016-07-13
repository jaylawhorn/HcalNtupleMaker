# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --filein dbs:/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISummer15GS-magnetOn_MCRUN2_71_V1-v1/GEN-SIM --fileout file:JME-RunIISpring16DR80-00006_step1.root --mc --eventcontent RAWSIM --pileup NoPileUp --datatier GEN-SIM-RAW --conditions 80X_mcRun2_asymptotic_2016_v3 --step DIGI,L1,DIGI2RAW --era Run2_25ns --no_exec -n 82
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('HCALNTUPLE',eras.Run2_25ns)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
#    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/QCD_Pt-15To7000_IsoTrkFilter_TuneCUETP8M1_Flat_13TeV-pythia8/GEN-SIM-RECO/NoPU_RECO_80X_mcRun2_asymptotic_2016_v3-v1/00000/02A03576-9E30-E611-8853-44A842CFC9BF.root'),
#    fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/j/jlawhorn/ISOTRK_DAMMIT.root'),
                            fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIISummer15GS/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/GEN-SIM/magnetOn_MCRUN2_71_V1-v1/10000/00572292-615B-E511-BDE8-C81F66B73923.root'),
    inputCommands = cms.untracked.vstring('keep *', 
        'drop *_genParticles_*_*', 
        'drop *_genParticlesForJets_*_*', 
        'drop *_kt4GenJets_*_*', 
        'drop *_kt6GenJets_*_*', 
        'drop *_iterativeCone5GenJets_*_*', 
        'drop *_ak4GenJets_*_*', 
        'drop *_ak7GenJets_*_*', 
        'drop *_ak8GenJets_*_*', 
        'drop *_ak4GenJetsNoNu_*_*', 
        'drop *_ak8GenJetsNoNu_*_*', 
        'drop *_genCandidatesForMET_*_*', 
        'drop *_genParticlesForMETAllVisible_*_*', 
        'drop *_genMetCalo_*_*', 
        'drop *_genMetCaloAndNonPrompt_*_*', 
        'drop *_genMetTrue_*_*', 
        'drop *_genMetIC5GenJs_*_*'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
)

from RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi import *
process.hcalOOTPileupESProducer = cms.ESProducer('OOTPileupDBCompatibilityESProducer')

from RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_hbhe_cfi import *
#from RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_ho_cfi import *
#from RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_hf_cfi import *
#from RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_zdc_cfi import *

hbheprereco.setNegativeFlags=False

process.hcalLocalRecoSequence = cms.Sequence(hbheprereco)

from RecoLocalCalo.HcalRecProducers.HBHEIsolatedNoiseReflagger_cfi import *
process.hcalGlobalRecoSequence = cms.Sequence(hbhereco)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:82'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string("HCALTree.root")
)
process.ExportTree = cms.EDAnalyzer("HcalAnalyzer",
  hbheInput = cms.InputTag('hbheprereco'),
  IsData = cms.untracked.bool(False),
  TriggerResults = cms.InputTag('TriggerResults','','HLT'),
  TotalChargeThreshold = cms.untracked.double(-9999)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_v3', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reco_step1 = cms.Path(process.hcalLocalRecoSequence)
process.reco_step2 = cms.Path(process.ExportTree)
process.endjob_step = cms.EndPath(process.endOfProcess)


# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.raw2digi_step,process.reco_step1,process.reco_step2,process.endjob_step)
