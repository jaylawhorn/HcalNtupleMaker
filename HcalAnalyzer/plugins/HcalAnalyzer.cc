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
#include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"

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

  edm::EDGetTokenT<HBHERecHitCollection> tok_hbhe_;
  edm::EDGetTokenT<HBHEDigiCollection> tok_hbhe_digi_;
  edm::EDGetTokenT<edm::PCaloHitContainer> tok_hbhe_sim_;
  
  bool FillHBHE;                  // Whether to store HBHE digi-level information or not                                                                      
  bool IsData_;  
  double TotalChargeThreshold;    // To avoid trees from overweight, only store digis above some threshold                                                    
  edm::Service<TFileService> FileService;

  const HcalDDDRecConstants *hcons;

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
  double Method3Energy[5184];
  double Method2Energy[5184];
  int IEta[5184];
  int IPhi[5184];
  int Depth[5184];
  
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
  
  tok_hbhe_ = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hbheInput"));
  tok_hbhe_digi_ = consumes<HBHEDigiCollection>(edm::InputTag("hcalDigis",""));
  tok_hbhe_sim_ = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","HcalHits"));

  IsData_ = iConfig.getUntrackedParameter<bool>("IsData");
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

   edm::ESHandle<HcalDDDRecConstants> pHRNDC;
   iSetup.get<HcalRecNumberingRecord>().get( pHRNDC );
   hcons = &(*pHRNDC);

   Handle<HBHERecHitCollection> hRecHits;
   iEvent.getByToken(tok_hbhe_, hRecHits);

   Handle<HBHEDigiCollection> hHBHEDigis;

   iEvent.getByToken(tok_hbhe_digi_, hHBHEDigis);

   Handle<PCaloHitContainer> hSimHits;
   iEvent.getByToken(tok_hbhe_sim_,hSimHits);

   ESHandle<HcalDbService> hConditions;
   iSetup.get<HcalDbRecord>().get(hConditions);

   ESHandle<CaloGeometry> hGeometry;
   iSetup.get<CaloGeometryRecord>().get(hGeometry);
   Geometry = hGeometry.product();

   // basic event coordinates                                                                                                                                                                               
   RunNumber = iEvent.id().run();
   EventNumber = iEvent.id().event();
   LumiSection = iEvent.luminosityBlock();
   Bunch = iEvent.bunchCrossing();
   Orbit = iEvent.orbitNumber();
   Time = iEvent.time().value();


   // loop over digis
   for(HBHEDigiCollection::const_iterator iter = hHBHEDigis->begin(); iter != hHBHEDigis->end(); iter++)
     {
       HcalDetId id = iter->id();

       //temporary hack to remove SiPM wedge at HEP17
       if (id.ieta()>14 && id.iphi()>=63 && id.iphi()<=66) continue;

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

       for (int x = 0; x<(int)hRecHits->size(); x++) {
	 HcalDetId id2 = (*hRecHits)[x].id();
	 if (id==id2) {
	   Method2Energy[PulseCount] = (*hRecHits)[x].energy();
	   Method3Energy[PulseCount] = (*hRecHits)[x].eaux();
	 }
       }
                             
       for(int i = 0; i < (int)iter->size(); i++)
	 {
	   const HcalQIESample &QIE = iter->sample(i);

	   Gain[PulseCount][i] = Calibrations.respcorrgain(QIE.capid());
	   Pedestal[PulseCount][i] = Calibrations.pedestal(QIE.capid());
	   Charge[PulseCount][i] = Tool[i] - Pedestal[PulseCount][i];
	 }

       if(!IsData_) {

	 double SamplingFactor = 1;
	 if(id.subdet() == HcalBarrel) {
	   SamplingFactor = simParameterMap_.hbParameters().samplingFactor(id);
	 } else if (id.subdet() == HcalEndcap) {
	   SamplingFactor = simParameterMap_.heParameters().samplingFactor(id);
	 }

	 for (int j = 0; j < (int) hSimHits->size(); j++) {

	   HcalDetId simId = HcalHitRelabeller::relabel((*hSimHits)[j].id(), hcons);
	   //if ((*hSimHits)[j].time() < 0 || (*hSimHits)[j].time() > 40) continue;
	   
	   if (simId == id) { // && (simId.iphi() == id.iphi() && abs(id.ieta())<14) {
	     SimHitEnergy[PulseCount]+=SamplingFactor*((*hSimHits)[j].energy());
	     //std::cout << 
	   }
	   /*	   else if (simId.ieta() == id.ieta() && simId.iphi() == id.iphi() && abs(id.ieta())>17) {
		   if (id.depth()==1 && ((HcalDetId)(*hSimHits)[j].id()).depth()<5) {
		   SimHitEnergy[PulseCount]+=SamplingFactor*((*hSimHits)[j].energy());
		   }
		   else if (id.depth()==2 && ((HcalDetId)(*hSimHits)[j].id().depth()>4) {
		   SimHitEnergy[PulseCount]+=SamplingFactor*((*hSimHits)[j].energy());
		   }*/
		   
		   //std::cout << "RecHit: " << id.ieta() << ", "<< id.iphi() << ", " << id.depth() << "; " << std::endl;
		   //std::cout << "SimHit: " << ((HcalDetId)(*hSimHits)[j].id()).ieta() << ", " << ((HcalDetId)(*hSimHits)[j].id()).iphi() << ", " << ((HcalDetId)(*hSimHits)[j].id()).depth() << std::endl;
		   //if ((HcalDetId)(*hSimHits)[j].id() == id) std::cout << "matches" << std::endl;
		   //std::cout << std::endl;
	   //}
	 }
       }
      
       
       IEta[PulseCount] = id.ieta();
       IPhi[PulseCount] = id.iphi();
       Depth[PulseCount] = id.depth();

       PulseCount = PulseCount + 1;
     }

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

  if(FillHBHE == true)
    {
      OutputTree->Branch("PulseCount", &PulseCount, "PulseCount/I");
      OutputTree->Branch("Charge", &Charge, "Charge[5184][10]/D");
      OutputTree->Branch("Pedestal", &Pedestal, "Pedestal[5184][10]/D");
      OutputTree->Branch("Gain", &Gain, "Gain[5184][10]/D");
      OutputTree->Branch("Method2Energy", &Method2Energy, "Method2Energy[5184]/D");
      OutputTree->Branch("Method3Energy", &Method3Energy, "Method3Energy[5184]/D");
      
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

      Method2Energy[i]=0;
      Method3Energy[i]=0;

      if(!IsData_) SimHitEnergy[i] = 0;
      IEta[i] = 0;
      IPhi[i] = 0;
      Depth[i] = 0;

    }

}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalAnalyzer);
