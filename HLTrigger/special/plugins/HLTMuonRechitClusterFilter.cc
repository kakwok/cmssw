#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"

#include "DataFormats/MuonReco/interface/MuonRecHitCluster.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

// system include files
#include <vector>
#include <map>
#include <iostream>
#include <memory>

class HLTMuonRechitClusterFilter : public HLTFilter {
public:
  explicit HLTMuonRechitClusterFilter(const edm::ParameterSet&);
  ~HLTMuonRechitClusterFilter() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  bool hltFilter(edm::Event&,
                 const edm::EventSetup&,
                 trigger::TriggerFilterObjectWithRefs& filterproduct) const override;

  edm::EDGetTokenT<reco::MuonRecHitClusterCollection> cluster_token_;
  edm::InputTag cluster_tag_;
  int min_N_;
  int min_Size_;
  int min_SizeMinusMB1_;
  int max_nMB1_;
  int max_nMB2_;
  int max_nME11_;
  int max_nME12_;
  int max_nME41_;
  int max_nME42_;
  int min_nStation_;
  double min_avgStation_;
  double min_Time_;
  double max_Time_;
  double min_Eta_;
  double max_Eta_;
  double max_TimeSpread_;
};
//
// constructors and destructor
//
HLTMuonRechitClusterFilter::HLTMuonRechitClusterFilter(const edm::ParameterSet& iConfig)
    : HLTFilter(iConfig),
      cluster_tag_(iConfig.getParameter<edm::InputTag>("ClusterTag")),
      min_N_(iConfig.getParameter<int>("MinN")),
      min_Size_(iConfig.getParameter<int>("Min_Size")),
      min_SizeMinusMB1_(iConfig.getParameter<int>("Min_SizeMinusMB1")),
      max_nMB1_(iConfig.getParameter<int>("Max_nMB1")),
      max_nMB2_(iConfig.getParameter<int>("Max_nMB2")),
      max_nME11_(iConfig.getParameter<int>("Max_nME11")),
      max_nME12_(iConfig.getParameter<int>("Max_nME12")),
      max_nME41_(iConfig.getParameter<int>("Max_nME41")),
      max_nME42_(iConfig.getParameter<int>("Max_nME42")),
      min_nStation_(iConfig.getParameter<int>("Min_nStation")),
      min_avgStation_(iConfig.getParameter<double>("Min_avgStation")),
      min_Time_(iConfig.getParameter<double>("Min_Time")),
      max_Time_(iConfig.getParameter<double>("Max_Time")),
      min_Eta_(iConfig.getParameter<double>("Min_Eta")),
      max_Eta_(iConfig.getParameter<double>("Max_Eta")),
      max_TimeSpread_(iConfig.getParameter<double>("Max_TimeSpread")) {
  cluster_token_ = consumes<reco::MuonRecHitClusterCollection>(cluster_tag_);
}

HLTMuonRechitClusterFilter::~HLTMuonRechitClusterFilter() = default;

void HLTMuonRechitClusterFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  makeHLTFilterDescription(desc);
  desc.add<edm::InputTag>("ClusterTag", edm::InputTag("hltCSCrechitClusters"));
  desc.add<int>("MinN", 1);
  desc.add<int>("Min_Size", 50);
  desc.add<int>("Min_SizeMinusMB1", 0);
  desc.add<int>("Max_nMB1", 0);
  desc.add<int>("Max_nMB2", 0);
  desc.add<int>("Max_nME11", 0);
  desc.add<int>("Max_nME12", 0);
  desc.add<int>("Max_nME41", 0);
  desc.add<int>("Max_nME42", 0);
  desc.add<int>("Min_nStation", 0);
  desc.add<double>("Min_avgStation", 0.0);
  desc.add<double>("Min_Time", -999);
  desc.add<double>("Max_Time", 999);
  desc.add<double>("Min_Eta", -1.0);
  desc.add<double>("Max_Eta", -1.0);
  desc.add<double>("Max_TimeSpread", 999);
  descriptions.add("hltMuonRechitClusterFilter", desc);
}

//
// member functions
//

// ------------ method called on each new Event  ------------
bool HLTMuonRechitClusterFilter::hltFilter(edm::Event& iEvent,
                                           const edm::EventSetup& iSetup,
                                           trigger::TriggerFilterObjectWithRefs& filterproduct) const {
  using namespace edm;
  using namespace std;
  using namespace trigger;

  int nClusterPassed = 0;

  auto const& rechitClusters = iEvent.get(cluster_token_);

  for (auto const& cluster : rechitClusters) {
    if ((cluster.size() >= min_Size_) && ((cluster.size() - cluster.nMB1()) >= min_SizeMinusMB1_) &&
        (cluster.nMB1() <= max_nMB1_) && (cluster.nMB2() <= max_nMB2_) && (cluster.nME11() <= max_nME11_) &&
        (cluster.nME12() <= max_nME12_) && (cluster.nME41() <= max_nME41_) && (cluster.nME42() <= max_nME42_) &&
        (cluster.nStation() >= min_nStation_) && (cluster.avgStation() >= min_avgStation_) &&
        ((min_Eta_ < 0.0) || (std::abs(cluster.eta()) >= min_Eta_)) &&
        ((max_Eta_ < 0.0) || (std::abs(cluster.eta()) <= max_Eta_)) && (cluster.time() > min_Time_) &&
        (cluster.time() <= max_Time_) && (cluster.timeSpread() <= max_TimeSpread_)) {
      nClusterPassed++;
    }
  }

  return (nClusterPassed >= min_N_);
}

// define as a framework module
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTMuonRechitClusterFilter);
