#ifndef HLTCSCRechitClusterFilter_h
#define HLTCSCRechitClusterFilter_h
// -*- C++ -*-
//
// Package:    HLTCSCRechitClusterFilter
// Class:      HLTCSCRechitClusterFilter
//
/**\class HLTCSCRechitClusterFilter HLTCSCRechitClusterFilter.cc filter/HLTCSCRechitClusterFilter/src/HLTCSCRechitClusterFilter.cc

Description: Filter to select HCAL abort gap events

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Carlo Battilana
//         Created:  Tue Jan 22 13:55:00 CET 2008
//
//

// include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"

#include "DataFormats/MuonReco/interface/MuonCSCRecHitCluster.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

//
// class declaration
//

class HLTCSCRechitClusterFilter : public HLTFilter {
public:
  explicit HLTCSCRechitClusterFilter(const edm::ParameterSet&);
  ~HLTCSCRechitClusterFilter() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  bool hltFilter(edm::Event&,
                 const edm::EventSetup&,
                 trigger::TriggerFilterObjectWithRefs& filterproduct) const override;

  edm::EDGetTokenT<reco::MuonCSCRecHitClusterCollection> m_cscCluster_token;
  edm::InputTag m_cscCluster_tag;
  int min_N_;
  int min_Size_;
  int max_nME11_;
  int max_nME12_;
  int min_nStation10_;
  double min_avgStation10_;
  double min_Time_;
  double max_Time_;
  double min_Eta_;
  double max_Eta_;
  double max_TimeSpread_;
};

#endif
