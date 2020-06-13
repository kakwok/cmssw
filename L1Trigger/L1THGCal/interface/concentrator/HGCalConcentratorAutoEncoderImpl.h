#ifndef __L1Trigger_L1THGCal_HGCalConcentratorAutoEncoderImpl_h__
#define __L1Trigger_L1THGCal_HGCalConcentratorAutoEncoderImpl_h__

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"
#include <vector>

#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"

class HGCalConcentratorAutoEncoderImpl {
public:
  HGCalConcentratorAutoEncoderImpl(const edm::ParameterSet& conf);

  void select(const std::vector<l1t::HGCalTriggerCell>& trigCellVecInput,
              std::vector<l1t::HGCalTriggerCell>& trigCellVecOutput);

  void eventSetup(const edm::EventSetup& es) { triggerTools_.eventSetup(es); }

private:
  std::vector<int> cellRemap_;

  int ae_outputCellU_[48];
  int ae_outputCellV_[48];

  HGCalTriggerTools triggerTools_;
};

#endif
