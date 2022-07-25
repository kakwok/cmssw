// -*- C++ -*-
//
// Package:    L1Trigger/CSCTriggerPrimitives
// Class:      SimpleCSCshowerFilter
//
/**\class SimpleCSCshowerFilter SimpleCSCshowerFilter.cc L1Trigger/CSCTriggerPrimitives/plugins/SimpleCSCshowerFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ka Hei Martin Kwok
//         Created:  Thu, 07 Jul 2022 16:53:22 GMT
//
//

// system include files
#include <memory>
#include "TTree.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"
#include "DataFormats/MuonReco/interface/MuonRecHitCluster.h"
#include "DataFormats/CSCDigi/interface/CSCShowerDigiCollection.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"

#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

//
// class declaration
//

class SimpleCSCshowerFilter : public edm::stream::EDFilter<> {
public:
  explicit SimpleCSCshowerFilter(const edm::ParameterSet&);
  ~SimpleCSCshowerFilter();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual bool filter(edm::Event&, const edm::EventSetup&) override;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  
  void reset();
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  //
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> geometryToken_;
  edm::EDGetTokenT<CSCRecHit2DCollection> inputToken_;
  const edm::EDGetTokenT<CSCShowerDigiCollection> dataLCTShower_token_;
  const edm::EDGetTokenT<reco::MuonRecHitClusterCollection> ca4CSCrechitClusters_token_;
  bool AsL1filter;
  bool AsRecofilter;
  TTree* hmtTree;
  int passL1;
  int hasCluster;

#define CSCRECHITARRAYSIZE 100000
  int ncscRechits;
  float cscRechitsPhi[CSCRECHITARRAYSIZE];
  float cscRechitsEta[CSCRECHITARRAYSIZE];
  float cscRechitsX[CSCRECHITARRAYSIZE];
  float cscRechitsY[CSCRECHITARRAYSIZE];
  float cscRechitsZ[CSCRECHITARRAYSIZE];
  float cscRechitsE[CSCRECHITARRAYSIZE];
  float cscRechitsTpeak[CSCRECHITARRAYSIZE];
  float cscRechitsTwire[CSCRECHITARRAYSIZE];
  int   cscRechitsQuality[CSCRECHITARRAYSIZE];
  int   cscRechitsChannels[CSCRECHITARRAYSIZE];
  int   cscRechitsNStrips[CSCRECHITARRAYSIZE];
  int   cscRechitsHitWire[CSCRECHITARRAYSIZE];
  int   cscRechitsWGroupsBX[CSCRECHITARRAYSIZE];
  int   cscRechitsNWireGroups[CSCRECHITARRAYSIZE];
  int   cscRechitsDetId[CSCRECHITARRAYSIZE];
  int   cscRechitsChamber[CSCRECHITARRAYSIZE];
  int   cscRechitsStation[CSCRECHITARRAYSIZE]; 
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SimpleCSCshowerFilter::SimpleCSCshowerFilter(const edm::ParameterSet& iConfig) :
  //now do what ever initialization is needed
  geometryToken_(esConsumes<CSCGeometry, MuonGeometryRecord>()),
  inputToken_(consumes<CSCRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("recHitLabel"))),
  dataLCTShower_token_(consumes(iConfig.getParameter<edm::InputTag>("dataLCTShower"))),
  ca4CSCrechitClusters_token_(consumes(iConfig.getParameter<edm::InputTag>("ca4CSCrechitClusters")))
{
  edm::Service<TFileService> fs;
  hmtTree = fs->make<TTree>("hmt", "HMT tree");
  AsL1filter   = iConfig.getParameter<bool>("AsL1filter");
  AsRecofilter = iConfig.getParameter<bool>("AsRecofilter");
}

SimpleCSCshowerFilter::~SimpleCSCshowerFilter() {
}

//
// member functions
//

// ------------ method called on each new Event  ------------
bool SimpleCSCshowerFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  bool has_LCTshs=false;
  bool has_CSCclusters=false;
  edm::Handle<CSCShowerDigiCollection> dataLCTshs;
  edm::Handle<reco::MuonRecHitClusterCollection> ca4CSCrechitClusters;

  auto const& geo = iSetup.getData(geometryToken_);
  auto const& rechits = iEvent.get(inputToken_);
  iEvent.getByToken(dataLCTShower_token_, dataLCTshs);
  iEvent.getByToken(ca4CSCrechitClusters_token_, ca4CSCrechitClusters);

  reset();

  ncscRechits = rechits.size();
  //std::cout << "ncscRechits  "<<ncscRechits<<std::endl;
  for (auto const& rechit : rechits) {

    LocalPoint recHitLocalPosition = rechit.localPosition();
    auto detid = rechit.cscDetId();;
    auto thischamber = geo.chamber(detid);
    int endcap = CSCDetId::endcap(detid) == 1 ? 1 : -1;
    if (thischamber) {
      GlobalPoint globalPosition = thischamber->toGlobal(recHitLocalPosition);

      cscRechitsX[ncscRechits] = globalPosition.x();
      cscRechitsY[ncscRechits] = globalPosition.y();
      cscRechitsZ[ncscRechits] = globalPosition.z();
      cscRechitsPhi[ncscRechits] = globalPosition.phi();
      cscRechitsEta[ncscRechits] = globalPosition.eta();
      cscRechitsE[ncscRechits] = rechit.energyDepositedInLayer();//not saved
      cscRechitsTpeak[ncscRechits] = rechit.tpeak();
      cscRechitsTwire[ncscRechits] = rechit.wireTime();
      cscRechitsQuality[ncscRechits] = rechit.quality();
      cscRechitsChamber[ncscRechits] = endcap * (CSCDetId::station(detid)*10 + CSCDetId::ring(detid));
      if (CSCDetId::ring(detid) == 4) cscRechitsChamber[ncscRechits] = endcap * (CSCDetId::station(detid)*10 + 1);
    }
  }

  const int min_endcap = CSCDetId::minEndcapId();
  const int max_endcap = CSCDetId::maxEndcapId();
  const int min_station = CSCDetId::minStationId();
  const int max_station = CSCDetId::maxStationId();
  const int min_sector = CSCTriggerNumbering::minTriggerSectorId();
  const int max_sector = CSCTriggerNumbering::maxTriggerSectorId();
  const int min_subsector = CSCTriggerNumbering::minTriggerSubSectorId();
  const int max_subsector = CSCTriggerNumbering::maxTriggerSubSectorId();
  const int min_chamber = CSCTriggerNumbering::minTriggerCscId();
  const int max_chamber = CSCTriggerNumbering::maxTriggerCscId();

  if (ca4CSCrechitClusters->size()>=1) has_CSCclusters=true;
  for (int endc = min_endcap; endc <= max_endcap; endc++) {
    // loop on all stations
    for (int stat = min_station; stat <= max_station; stat++) {
      int numsubs = ((stat == 1) ? max_subsector : 1);
      // loop on sectors and subsectors
      for (int sect = min_sector; sect <= max_sector; sect++) {
        for (int subs = min_subsector; subs <= numsubs; subs++) {
          // loop on all chambers
          for (int cham = min_chamber; cham <= max_chamber; cham++) {
            // extract the ring number
            int ring = CSCTriggerNumbering::ringFromTriggerLabels(stat, cham);

            // actual chamber number =/= trigger chamber number
            int chid = CSCTriggerNumbering::chamberFromTriggerLabels(sect, subs, stat, cham);

            // 0th layer means whole chamber.
            CSCDetId detid(endc, stat, ring, chid, 0);
            
            auto range_dataLCTshs = dataLCTshs->get(detid);
            for (auto dlct = range_dataLCTshs.first; dlct != range_dataLCTshs.second; dlct++) {
               if (dlct->isNominalInTime()) has_LCTshs=true;
            }
          }
        }
      }
    }
  }
  hasCluster = int(has_CSCclusters);
  passL1 = int(has_LCTshs);
  if(AsL1filter){
      if (has_LCTshs){
          hmtTree->Fill();
      }
      return has_LCTshs;
  } else if(AsRecofilter){
      if (has_CSCclusters){
          hmtTree->Fill();
      }
      return has_CSCclusters;
  }else
    return has_LCTshs;
}

void SimpleCSCshowerFilter::reset(){
  ncscRechits = 0;
  for ( int i = 0; i < CSCRECHITARRAYSIZE; i++) {
    cscRechitsPhi[i] = 0.0;
    cscRechitsEta[i] = 0.0;
    cscRechitsX[i] = 0.0;
    cscRechitsY[i] = 0.0;
    cscRechitsZ[i] = 0.0;
    cscRechitsE[i] = 0.0;
    cscRechitsTwire[i] = 0.0;
    cscRechitsTpeak[i] = 0.0;
    cscRechitsQuality[i] = 0.0;
    cscRechitsChannels[i] = 0;
    cscRechitsNStrips[i] = 0;
    cscRechitsHitWire[i] = 0;
    cscRechitsWGroupsBX[i] = 0;
    cscRechitsNWireGroups[i] = 0;
    cscRechitsDetId[i] = 0;
    cscRechitsStation[i] = 0;
    cscRechitsChamber[i] = 0;
  }
}


void
SimpleCSCshowerFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
  hmtTree->Branch("passL1",&passL1,"passL1/I");
  hmtTree->Branch("hasCluster",&hasCluster,"hasCluster/I");
  hmtTree->Branch("ncscRechits",&ncscRechits,"ncscRechits/I");
  hmtTree->Branch("cscRechitsPhi",&cscRechitsPhi,"cscRechitsPhi[ncscRechits]/F");
  hmtTree->Branch("cscRechitsEta",&cscRechitsEta,"cscRechitsEta[ncscRechits]/F");
  hmtTree->Branch("cscRechitsX",&cscRechitsX,"cscRechitsX[ncscRechits]/F");
  hmtTree->Branch("cscRechitsY",&cscRechitsY,"cscRechitsY[ncscRechits]/F");
  hmtTree->Branch("cscRechitsZ",&cscRechitsZ,"cscRechitsZ[ncscRechits]/F");
  hmtTree->Branch("cscRechitsE",&cscRechitsE,"cscRechitsE[ncscRechits]/F");
  hmtTree->Branch("cscRechitsTpeak",&cscRechitsTpeak,"cscRechitsTpeak[ncscRechits]/F");
  hmtTree->Branch("cscRechitsTwire",&cscRechitsTwire,"cscRechitsTwire[ncscRechits]/F");
  hmtTree->Branch("cscRechitsQuality",&cscRechitsQuality,"cscRechitsQuality[ncscRechits]/I");
  hmtTree->Branch("cscRechitsChamber",&cscRechitsChamber,"cscRechitsChamber[ncscRechits]/I");
  hmtTree->Branch("cscRechitsStation",&cscRechitsStation,"cscRechitsStation[ncscRechits]/I");

}

// ------------ method called when ending the processing of a run  ------------
/*
void
SimpleCSCshowerFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
SimpleCSCshowerFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
SimpleCSCshowerFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SimpleCSCshowerFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(SimpleCSCshowerFilter);
