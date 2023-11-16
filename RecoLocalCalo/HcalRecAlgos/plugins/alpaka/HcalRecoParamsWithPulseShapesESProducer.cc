#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/HcalObjects/interface/alpaka/HcalRecoParamWithPulseShapePortable.h"
#include "CondFormats/HcalObjects/interface/HcalRecoParamWithPulseShapeSoA.h"
#include "CondFormats/DataRecord/interface/HcalRecoParamsRcd.h"

#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseShapes.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/PulseShapeFunctor.h"

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ModuleFactory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"



namespace ALPAKA_ACCELERATOR_NAMESPACE {
  class HcalRecoParamWithPulseShapeESProducer : public ESProducer {
  public:
    HcalRecoParamWithPulseShapeESProducer(edm::ParameterSet const& iConfig) : ESProducer(iConfig) {
      auto cc = setWhatProduced(this);
      recoParamsToken_ = cc.consumes();
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      descriptions.addWithDefaultLabel(desc);
    }

    std::unique_ptr<HcalRecoParamWithPulseShapeHost> produce(HcalRecoParamsRcd const& iRecord) {

      auto const& recoParams     = iRecord.get(recoParamsToken_);

      auto const      containers = recoParams.getAllContainers();
      size_t const totalChannels = recoParams.getAllContainers()[0].second.size() + recoParams.getAllContainers()[1].second.size();

      //Get unique ids
      HcalPulseShapes pulseShapes;
      std::unordered_map<unsigned int, uint32_t> idCache; //<pulseShapeId,arrIdx>
      uint32_t unique_ids = 0;

      auto const& barrelValues = containers[0].second;
      for (uint64_t i = 0; i < barrelValues.size(); ++i) {
          auto const pulseShapeId = barrelValues[i].pulseShapeID();
          if (pulseShapeId == 0)   continue;
          if (auto const iter = idCache.find(pulseShapeId); iter == idCache.end()) {
            unique_ids++;
            // new guy
            idCache[pulseShapeId] = unique_ids;
          }
      }
      auto const& endcapValues = containers[1].second;
      for (uint64_t i = 0; i < endcapValues.size(); ++i) {
        auto const pulseShapeId = endcapValues[i].pulseShapeID();
        if (auto const iter = idCache.find(pulseShapeId); iter == idCache.end()) {
          if (pulseShapeId == 0)   continue;
            unique_ids++;
            idCache[pulseShapeId] = unique_ids;
        }
      }

      auto product = std::make_unique<HcalRecoParamWithPulseShapeHost>(totalChannels,unique_ids);
      auto recoView       = product->recoParam().view();
      auto pulseShapeView = product->pulseShape().view();
      for (uint64_t i = 0; i < barrelValues.size(); ++i) {
        recoView[i].param1()  = barrelValues[i].param1;
        recoView[i].param2()  = barrelValues[i].param2;
        recoView[i].ids()     = barrelValues[i].pulseShapeID();
      }
      //fill pulseShape views
      for (auto& it: idCache) {
          auto const pulseShapeId = it.first;
          auto const arrId = it.second;
          auto const& pulseShape = pulseShapes.getShape(pulseShapeId);                                    
          FitterFuncs::PulseShapeFunctor functor{pulseShape, false, false, false, 1, 0, 0, hcal::constants::maxSamples};  

          for (int i = 0; i < hcal::constants::maxPSshapeBin; i++) {
             acc25nsVec_[offset256 * numShapes + i] = functor.acc25nsVec()[i];
             pulseShapeView[arrId].acc25nsVec[i]  = functor.acc25nsVec()[i]
             pulseShapeView[arrId].diff25nsItvlVec[i]  = functor.diff25nsItvlVec()[i]
           }
      }
 

      return product;
    }

  private:
    edm::ESGetToken<HcalRecoParams , HcalRecoParamsRcd > recoParamsToken_;

  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_EVENTSETUP_ALPAKA_MODULE(HcalRecoParamWithPulseShapeESProducer);
