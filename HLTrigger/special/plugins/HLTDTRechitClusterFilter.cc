// -*- C++ -*-
//
// Package:    HLTDTRechitClusterFilter
// Class:      HLTDTRechitClusterFilter
//
/**\class HLTDTRechitClusterFilter HLTDTRechitClusterFilter.cc filter/HLTDTRechitClusterFilter/src/HLTDTRechitClusterFilter.cc

Description:

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Carlo Battilana
//         Created:  Tue Jan 22 13:55:00 CET 2008
//
//

// system include files
#include <vector>
#include <map>
#include <iostream>
#include <memory>

// user include files
#include "HLTDTRechitClusterFilter.h"

//
// constructors and destructor
//
HLTDTRechitClusterFilter::HLTDTRechitClusterFilter(const edm::ParameterSet& iConfig)
    : HLTFilter(iConfig),
      m_dtCluster_tag(iConfig.getParameter<edm::InputTag>("dtClusterTag")),
      min_N_(iConfig.getParameter<int>("MinN")),
      min_Size_(iConfig.getParameter<int>("Min_Size")),
      max_nMB1_(iConfig.getParameter<int>("Max_nMB1")),
      max_nMB2_(iConfig.getParameter<int>("Max_nMB2")),
      min_nStation10_(iConfig.getParameter<int>("Min_nStation10")),
      min_avgStation10_(iConfig.getParameter<double>("Min_avgStation10")),
      min_Eta_(iConfig.getParameter<double>("Min_Eta")),
      max_Eta_(iConfig.getParameter<double>("Max_Eta"))
       {
  m_dtCluster_token = consumes<reco::MuonDTRecHitClusterCollection>(m_dtCluster_tag);
}

HLTDTRechitClusterFilter::~HLTDTRechitClusterFilter() = default;

void HLTDTRechitClusterFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  makeHLTFilterDescription(desc);
  desc.add<edm::InputTag>("dtClusterTag", edm::InputTag("hltdtRechitClusters"));
  desc.add<int>("MinN", 1);
  desc.add<int>("Min_Size", 50);
  desc.add<int>("Max_nMB1", 0);
  desc.add<int>("Max_nMB2", 0);
  desc.add<int>("Min_nStation10", 0);
  desc.add<double>("Min_avgStation10", 0.0);
  desc.add<double>("Min_Eta", -1.0);
  desc.add<double>("Max_Eta", -1.0);
  descriptions.add("hltDTRechitClusterFilter", desc);
}

//
// member functions
//

// ------------ method called on each new Event  ------------
bool HLTDTRechitClusterFilter::hltFilter(edm::Event& iEvent,
                                     const edm::EventSetup& iSetup,
                                     trigger::TriggerFilterObjectWithRefs& filterproduct) const {
  using namespace edm;
  using namespace std;
  using namespace trigger;

  int nClusterPassed = 0;

  edm::Handle<reco::MuonDTRecHitClusterCollection> dtRechitClusters;
  iEvent.getByToken(m_dtCluster_token, dtRechitClusters);

  for (auto cluster : *dtRechitClusters) {
    if ( (cluster.size() >= min_Size_) &&
        (cluster.nMB1() <= max_nMB1_) &&
        (cluster.nMB2() <= max_nMB2_) &&
        (cluster.nStation10() >= min_nStation10_) &&
        (cluster.avgStation10() >= min_avgStation10_) &&
        ((min_Eta_ < 0.0) || (std::abs(cluster.eta()) >= min_Eta_)) &&
        ((max_Eta_ < 0.0) || (std::abs(cluster.eta()) <= max_Eta_)) 
        ){
        nClusterPassed++;
    }
  }

  return (nClusterPassed >= min_N_);
}

// define as a framework module
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTDTRechitClusterFilter);
