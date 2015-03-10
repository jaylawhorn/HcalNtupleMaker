import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations                                             
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/data/Run2012D/JetHT/RAW/v1/000/208/352/D22C6707-A03B-E211-9F28-003048D37538.root'),
    skipEvents= cms.untracked.uint32(0)
)

process.options = cms.untracked.PSet(
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
   version = cms.untracked.string('$Revision: 1.19 $'),
   annotation = cms.untracked.string('RECO nevts:1'),
   name = cms.untracked.string('Applications')
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_data', '')

# Hcal noise analyzers
process.HBHENoiseFilterResultProducer = cms.EDProducer(
   'HBHENoiseFilterResultProducer',
   noiselabel = cms.InputTag('hcalnoise'),
   minRatio = cms.double(-999),
   maxRatio = cms.double(999),
   minHPDHits = cms.int32(17), 
   minRBXHits = cms.int32(999),
   minHPDNoOtherHits = cms.int32(10),
   minZeros = cms.int32(10),
   minHighEHitTime = cms.double(-9999.0),
   maxHighEHitTime = cms.double(9999.0),
   maxRBXEMF = cms.double(-999.0),
   minNumIsolatedNoiseChannels = cms.int32(10),
   minIsolatedNoiseSumE = cms.double(50.0),
   minIsolatedNoiseSumEt = cms.double(25.0),
   useTS4TS5 = cms.bool(False),
   useRBXRechitR45Loose = cms.bool(False),
   useRBXRechitR45Tight = cms.bool(False),
   IgnoreTS4TS5ifJetInLowBVRegion = cms.bool(True),
   jetlabel = cms.InputTag('ak5PFJets'),
   maxjetindex = cms.int32(0),
   maxNHF = cms.double(0.9)
)
process.TFileService = cms.Service("TFileService",
  fileName = cms.string("HCALTree.root")
)
process.ExportTree = cms.EDAnalyzer("HcalAnalyzer",
  HBHERecHitCollection = cms.untracked.string('hbhereco'),
  IsData = cms.untracked.bool(True),
  TotalChargeThreshold = cms.untracked.double(-9999)
)

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction * process.HBHENoiseFilterResultProducer * process.ExportTree)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.endjob_step)
