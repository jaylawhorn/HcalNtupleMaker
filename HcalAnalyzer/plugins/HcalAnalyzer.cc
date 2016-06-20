// -*- C++ -*-
//
// Package:    Hcal/HcalAnalyzer
// Class:      HcalAnalyzer
// 
/**\class HcalAnalyzer HcalAnalyzer.cc Hcal/HcalAnalyzer/plugins/HcalAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Katharina Bierwagen
//         Created:  Thu, 05 Feb 2015 11:36:53 GMT
//
//


// system include files
#include <memory>
#include <string>
#include <map>
#include <iostream>
using namespace std;

//---------------------------------------------------------------------------   
#include "TTree.h"
#include "TFile.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSimpleRecAlgo.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseContainmentManager.h"
#include "CondFormats/HcalObjects/interface/HcalRecoParams.h"
#include "CondFormats/HcalObjects/interface/HcalRecoParam.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalDbASCIIIO.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseShapes.h"

#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalSimParameterMap.h"

#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/METReco/interface/HcalNoiseRBX.h"
#include "RecoMET/METAlgorithms/interface/HcalHPDRBXMap.h"

// Trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"

#include "FWCore/Common/interface/TriggerNames.h"

const int NUM_TRIGGERS_MAX = 500;

//--------------------------------------------------------------------------- 
class HcalNtuplelizer;
//--------------------------------------------------------------------------- 

class HcalAnalyzer : public edm::EDAnalyzer {
public:
  explicit HcalAnalyzer(const edm::ParameterSet&);
  ~HcalAnalyzer();
  
  
private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  //virtual void endRun(edm::Run const&, edm::EventSetup const&);
  
  
private:
  std::map<int, double> hitEnergySumMap_;
  HcalSimParameterMap simParameterMap_;

  //edm::EDGetTokenT<HBHERecHitCollection> tok_hbhe_;
  edm::EDGetTokenT<HBHERecHitCollection> tok_hbhe_;
  edm::EDGetTokenT<HBHEDigiCollection> tok_hbhe_digi_;

  bool FillHBHE;                  // Whether to store HBHE digi-level information or not                                                                      
  bool IsData_;  
  double TotalChargeThreshold;    // To avoid trees from overweight, only store digis above some threshold                                                    
  edm::Service<TFileService> FileService;

  edm::EDGetTokenT<edm::TriggerResults> triggerResults_;
  std::string triggerPath_;
  edm::InputTag triggerFilter_;

  // Basic event coordinates                                                   
  long long RunNumber;
  long long EventNumber;
  long long LumiSection;
  long long Bunch;
  long long Orbit;
  long long Time;
  
  // HBHE rechits and digis
  int PulseCount;
  double Charge[5184][10];
  double Pedestal[5184][10];
  double Gain[5184][10];
  double SimHitEnergy[5184];
  int IEta[5184];
  int IPhi[5184];
  int Depth[5184];

  // Summary variables for baseline Hcal noise filter                           
  //int HPDHits;
  //int HPDNoOtherHits;
  //int MaxZeros;
  //double MinE2E10;
  //double MaxE2E10;
  //bool HasBadRBXR45;
  //bool HasBadRBXRechitR45Loose;
  //bool HasBadRBXRechitR45Tight;

  // Official decision from the baseline hcal noise filter                      
  //bool OfficialDecision;
  
private:
  TTree *OutputTree;
  
  const CaloGeometry *Geometry;
  
  void ClearVariables();  
};


//
// constructors and destructor
//
HcalAnalyzer::HcalAnalyzer(const edm::ParameterSet& iConfig)  
{
  FillHBHE = iConfig.getUntrackedParameter<bool>("FillHBHE", true);
  TotalChargeThreshold = iConfig.getUntrackedParameter<double>("TotalChargeThreshold", 0);
  
  //tok_hbhe_ = consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>>(iConfig.getParameter<edm::InputTag>("hbheInput"));
  tok_hbhe_ = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hbheInput"));
  tok_hbhe_digi_ = consumes<HBHEDigiCollection>(edm::InputTag("hcalDigis",""));

  IsData_ = iConfig.getUntrackedParameter<bool>("IsData");

  triggerResults_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"));
}


HcalAnalyzer::~HcalAnalyzer()
{
}

// ------------ method called for each event  ------------
void
HcalAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   ClearVariables();

   // get stuff                                                                                                                                                                            
                 
   Handle<HBHERecHitCollection> hRecHits;
   //Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>> hRecHits;
   //iEvent.getByLabel(InputTag(sHBHERecHitCollection, "RERECO"), hRecHits);
   iEvent.getByToken(tok_hbhe_, hRecHits);

   Handle<HBHEDigiCollection> hHBHEDigis;
   //iEvent.getByLabel(InputTag("hcalDigis"), hHBHEDigis);

   iEvent.getByToken(tok_hbhe_digi_, hHBHEDigis);

   Handle<PCaloHitContainer> hSimHits;
   iEvent.getByLabel("g4SimHits","HcalHits",hSimHits);

   ESHandle<HcalDbService> hConditions;
   iSetup.get<HcalDbRecord>().get(hConditions);

   ESHandle<CaloGeometry> hGeometry;
   iSetup.get<CaloGeometryRecord>().get(hGeometry);
   Geometry = hGeometry.product();

   Handle<TriggerResults> hltresults;

   iEvent.getByToken(triggerResults_,hltresults);
   if (!hltresults.isValid()) {
     LogError("HcalAnalyzer") << "invalid collection!" << "\n";
     return;
   }

   const TriggerNames& trigNames = iEvent.triggerNames(*hltresults);
   
   //uint numTriggers = trigNames.size();

   //bool passTrig= false;
   //
   //if (passTrig==false) return;
   
   //loop over triggers
   //for( unsigned int hltIndex=0; hltIndex<numTriggers; ++hltIndex ){
   //  if (hltIndex > NUM_TRIGGERS_MAX) break; //avoid memory leaks from too many triggers
   //
   //  if (hltresults->wasrun(hltIndex)) {
   //    std::cout << hltIndex << " " << trigNames.triggerName(hltIndex) << endl;
   //  }
   //  //myTriggerNames->push_back(trigNames.triggerName(hltIndex));
   //  if (hltresults->wasrun(hltIndex) && hltresults->accept(hltIndex)) std::cout << "passed" << endl;
   //}

   // basic event coordinates                                                                                                                                                                               
   RunNumber = iEvent.id().run();
   EventNumber = iEvent.id().event();
   LumiSection = iEvent.luminosityBlock();
   Bunch = iEvent.bunchCrossing();
   Orbit = iEvent.orbitNumber();
   Time = iEvent.time().value();

   if(!IsData_) {
     // store the energy of each hit in a map
     hitEnergySumMap_.clear();
     PCaloHitContainer::const_iterator hitItr = hSimHits->begin();
     PCaloHitContainer::const_iterator last = hSimHits->end();
     
     for( ; hitItr != last; ++hitItr)
       {
	 HcalDetId hcalDetId(hitItr->id());
	 
	 if(hcalDetId.subdet()== HcalBarrel || hcalDetId.subdet() == HcalEndcap)
	   {
	     int id = hitItr->id();
	     double samplingFactor=1.;
	     if(hcalDetId.subdet()== HcalBarrel)
	       {
		 samplingFactor = simParameterMap_.hbParameters().samplingFactor(DetId(id));
	       }
	     else if(hcalDetId.subdet() == HcalEndcap)
	       {
		 samplingFactor = simParameterMap_.heParameters().samplingFactor(DetId(id));
	       }
	     double energy = hitItr->energy() * samplingFactor;
	     // add it to the map
	     std::map<int, double>::iterator mapItr = hitEnergySumMap_.find(id);
	     if(mapItr == hitEnergySumMap_.end()) {
	       hitEnergySumMap_[id] = energy;
	     }
	     else
	       {
		 hitEnergySumMap_[id] += energy;
	       }
	   }
       }
   }
   
   // HBHE rechit maps - we want to link rechits and digis together                                                                                                                                         
   map<HcalDetId, int> RecHitIndex;
   for(int i = 0; i < (int)hRecHits->size(); i++)
     {
       HcalDetId id = (*hRecHits)[i].id();
       RecHitIndex.insert(pair<HcalDetId, int>(id, i));
     }

   // loop over digis                                                                                                                                                                                       
   for(HBHEDigiCollection::const_iterator iter = hHBHEDigis->begin(); iter != hHBHEDigis->end(); iter++)
     {
       HcalDetId id = iter->id();

       // First let's convert ADC to deposited charge                                                                                                                                                       
       const HcalCalibrations &Calibrations = hConditions->getHcalCalibrations(id);
       const HcalQIECoder *ChannelCoder = hConditions->getHcalCoder(id);
       const HcalQIEShape *Shape = hConditions->getHcalShape(ChannelCoder);
       HcalCoderDb Coder(*ChannelCoder, *Shape);
       CaloSamples Tool;
       Coder.adc2fC(*iter, Tool);

       // Total charge of the digi                                                                                                                                                                           
       double TotalCharge = 0;
       for(int i = 0; i < (int)iter->size(); i++)
         TotalCharge = TotalCharge + Tool[i] - Calibrations.pedestal(iter->sample(i).capid());

       // If total charge is smaller than threshold, don't store this rechit/digi into the tree                                                                                                              
       if(TotalCharge < TotalChargeThreshold)
	 continue;

       // Safety check - there are only 5184 channels in HBHE, but just in case...                                                                                                                           
       if(PulseCount >= 5184)
	 {
	   PulseCount = PulseCount + 1;
	   continue;
	 }

       // Fill things into the tree                                                                                                                                                                          
       for(int i = 0; i < (int)iter->size(); i++)
	 {
	   const HcalQIESample &QIE = iter->sample(i);

	   Gain[PulseCount][i] = Calibrations.respcorrgain(QIE.capid());
	   Pedestal[PulseCount][i] = Calibrations.pedestal(QIE.capid());
	   Charge[PulseCount][i] = Tool[i] - Pedestal[PulseCount][i];
	 }

       if(!IsData_) SimHitEnergy[PulseCount]=hitEnergySumMap_[(*hRecHits)[RecHitIndex[id]].id().rawId()];

       IEta[PulseCount] = id.ieta();
       IPhi[PulseCount] = id.iphi();
       Depth[PulseCount] = id.depth();

       PulseCount = PulseCount + 1;
     }

   // hcal sumamry objects               

   //HPDHits = hSummary->maxHPDHits();
   //HPDNoOtherHits = hSummary->maxHPDNoOtherHits();
   //MaxZeros = hSummary->maxZeros();
   //MinE2E10 = hSummary->minE2Over10TS();
   //MaxE2E10 = hSummary->maxE2Over10TS();
   //HasBadRBXR45 = hSummary->HasBadRBXTS4TS5();
   //HasBadRBXRechitR45Loose = hSummary->HasBadRBXRechitR45Loose();
   //HasBadRBXRechitR45Tight = hSummary->HasBadRBXRechitR45Tight();

   // finally actually fill the tree

   OutputTree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
HcalAnalyzer::beginJob()
{
  // Make branches in the output trees                                                                                                                                                                     
  OutputTree = FileService->make<TTree>("HcalTree", "Hcal tree");

  OutputTree->Branch("RunNumber", &RunNumber, "RunNumber/LL");
  OutputTree->Branch("EventNumber", &EventNumber, "EventNumber/LL");
  OutputTree->Branch("LumiSection", &LumiSection, "LumiSection/LL");
  OutputTree->Branch("Bunch", &Bunch, "Bunch/LL");
  OutputTree->Branch("Orbit", &Orbit, "Orbit/LL");
  OutputTree->Branch("Time", &Time, "Time/LL");

  //OutputTree->Branch("HPDHits", &HPDHits, "HPDHits/I");
  //OutputTree->Branch("HPDNoOtherHits", &HPDNoOtherHits, "HPDNoOtherHits/I");
  //OutputTree->Branch("MaxZeros", &MaxZeros, "MaxZeros/I");
  //OutputTree->Branch("MinE2E10", &MinE2E10, "MinE2E10/D");
  //OutputTree->Branch("MaxE2E10", &MaxE2E10, "MaxE2E10/D");
  //OutputTree->Branch("HasBadRBXR45", &HasBadRBXR45, "HasBadRBXR45/O");
  //OutputTree->Branch("HasBadRBXRechitR45Loose", &HasBadRBXRechitR45Loose, "HasBadRBXRechitR45Loose/O");
  //OutputTree->Branch("HasBadRBXRechitR45Tight", &HasBadRBXRechitR45Tight, "HasBadRBXRechitR45Tight/O");

  //OutputTree->Branch("OfficialDecision", &OfficialDecision, "OfficialDecision/O");

  if(FillHBHE == true)
    {
      OutputTree->Branch("PulseCount", &PulseCount, "PulseCount/I");
      OutputTree->Branch("Charge", &Charge, "Charge[5184][10]/D");
      OutputTree->Branch("Pedestal", &Pedestal, "Pedestal[5184][10]/D");
      OutputTree->Branch("Gain", &Gain, "Gain[5184][10]/D");
      
      if(!IsData_) OutputTree->Branch("SimHitEnergy", &SimHitEnergy, "SimHitEnergy[5184]/D");
      
      OutputTree->Branch("IEta", &IEta, "IEta[5184]/I");
      OutputTree->Branch("IPhi", &IPhi, "IPhi[5184]/I");
      OutputTree->Branch("Depth", &Depth, "Depth[5184]/I");
    }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HcalAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
HcalAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
HcalAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

 //------------------------------------------------------------------------------                                                                                                                            
void HcalAnalyzer::ClearVariables()
{
  RunNumber = 0;
  EventNumber = 0;
  LumiSection = 0;
  Bunch = 0;
  Orbit = 0;
  Time = 0;

  PulseCount = 0;
  for(int i = 0; i < 5184; i++)
    {
      for(int j = 0; j < 10; j++)
	{
	  Gain[i][j] = 0;
	  Charge[i][j] = 0;
	  Pedestal[i][j] = 0;
	}

      if(!IsData_) SimHitEnergy[i] = 0;
      IEta[i] = 0;
      IPhi[i] = 0;
      Depth[i] = 0;

    }

  //HPDHits = 0;
  //HPDNoOtherHits = 0;
  //MaxZeros = 0;
  //MinE2E10 = 0;
  //MaxE2E10 = 0;
  //HasBadRBXR45 = false;
  //HasBadRBXRechitR45Loose = false;
  //HasBadRBXRechitR45Tight = false;

  //OfficialDecision = false;

}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalAnalyzer);
