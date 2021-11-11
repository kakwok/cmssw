// -*- C++ -*-
//
// Package:    HLTCSCRechitClusterFilter
// Class:      HLTCSCRechitClusterFilter
//
/**\class HLTCSCRechitClusterFilter HLTCSCRechitClusterFilter.cc filter/HLTCSCRechitClusterFilter/src/HLTCSCRechitClusterFilter.cc

Description:

Implementation:
<Notes on implementation>
*/
//
//
//

// system include files
#include <vector>
#include <map>
#include <iostream>
#include <memory>

// user include files
#include "HLTCSCRechitClusterFilter.h"

//
// constructors and destructor
//
HLTCSCRechitClusterFilter::HLTCSCRechitClusterFilter(const edm::ParameterSet& iConfig)
    : HLTFilter(iConfig),
      m_cscCluster_tag(iConfig.getParameter<edm::InputTag>("cscClusterTag")),
      min_N_(iConfig.getParameter<int>("MinN")),
      min_Size_(iConfig.getParameter<int>("Min_Size")),
      max_nME11_(iConfig.getParameter<int>("Max_nME11")),
      max_nME12_(iConfig.getParameter<int>("Max_nME12")),
      min_nStation10_(iConfig.getParameter<int>("Min_nStation10")),
      min_avgStation10_(iConfig.getParameter<double>("Min_avgStation10")),
      min_Time_(iConfig.getParameter<double>("Min_Time")),
      max_Time_(iConfig.getParameter<double>("Max_Time")),
      min_Eta_(iConfig.getParameter<double>("Min_Eta")),
      max_Eta_(iConfig.getParameter<double>("Max_Eta")),
      max_TimeSpread_(iConfig.getParameter<double>("Max_TimeSpread"))   {
  m_cscCluster_token = consumes<reco::MuonCSCRecHitClusterCollection>(m_cscCluster_tag);
}

HLTCSCRechitClusterFilter::~HLTCSCRechitClusterFilter() = default;

void HLTCSCRechitClusterFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  makeHLTFilterDescription(desc);
  desc.add<edm::InputTag>("cscClusterTag", edm::InputTag("hltcscRechitClusters"));
  desc.add<int>("MinN", 1);
  desc.add<int>("Min_Size", 50);
  desc.add<int>("Max_nME11", 0);
  desc.add<int>("Max_nME12", 0);
  desc.add<int>("Min_nStation10", 0);
  desc.add<double>("Min_avgStation10", 0.0);
  desc.add<double>("Min_Time", -999);
  desc.add<double>("Max_Time", 999);
  desc.add<double>("Min_Eta", -1.0);
  desc.add<double>("Max_Eta", -1.0);
  desc.add<double>("Max_TimeSpread", 999);
  descriptions.add("hltCSCRechitClusterFilter", desc);
}

//
// member functions
//

// ------------ method called on each new Event  ------------
bool HLTCSCRechitClusterFilter::hltFilter(edm::Event& iEvent,
                                     const edm::EventSetup& iSetup,
                                     trigger::TriggerFilterObjectWithRefs& filterproduct) const {
  using namespace edm;
  using namespace std;
  using namespace trigger;

  int nClusterPassed = 0;

  edm::Handle<reco::MuonCSCRecHitClusterCollection> cscRechitClusters;
  iEvent.getByToken(m_cscCluster_token, cscRechitClusters);

  for (auto cluster : *cscRechitClusters) {
    if ( (cluster.size() >= min_Size_) &&
        (cluster.nME11() <= max_nME11_) &&
        (cluster.nME12() <= max_nME12_) &&
        (cluster.nStation10() >= min_nStation10_) &&
        (cluster.avgStation10() >= min_avgStation10_) &&
        ((min_Eta_ < 0.0) || (std::abs(cluster.eta()) >= min_Eta_)) &&
        ((max_Eta_ < 0.0) || (std::abs(cluster.eta()) <= max_Eta_)) &&
        (cluster.time() >min_Time_) && (cluster.time() <=max_Time_) &&
        (cluster.timeSpread() <= max_TimeSpread_)){
        nClusterPassed++;
    }
  }

  return (nClusterPassed >= min_N_);
}

// define as a framework module
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTCSCRechitClusterFilter);
