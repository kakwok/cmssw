#ifndef HLTDTRechitClusterFilter_h
#define HLTDTRechitClusterFilter_h
// -*- C++ -*-
//
// Package:    HLTDTRechitClusterFilter
// Class:      HLTDTRechitClusterFilter
//
/**\class HLTDTRechitClusterFilter HLTDTRechitClusterFilter.cc filter/HLTDTRechitClusterFilter/src/HLTDTRechitClusterFilter.cc

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

#include "DataFormats/MuonReco/interface/MuonDTRecHitCluster.h"

//
// class declaration
//

class HLTDTRechitClusterFilter : public HLTFilter {
public:
  explicit HLTDTRechitClusterFilter(const edm::ParameterSet&);
  ~HLTDTRechitClusterFilter() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  bool hltFilter(edm::Event&,
                 const edm::EventSetup&,
                 trigger::TriggerFilterObjectWithRefs& filterproduct) const override;

  edm::EDGetTokenT<reco::MuonDTRecHitClusterCollection> m_dtCluster_token;
  edm::InputTag m_dtCluster_tag;
  int min_N_;
  int min_Size_;
  int max_nMB1_;
  int max_nMB2_;
  int min_nStation10_;
  double min_avgStation10_;
  double min_Eta_;
  double max_Eta_;
};

#endif
