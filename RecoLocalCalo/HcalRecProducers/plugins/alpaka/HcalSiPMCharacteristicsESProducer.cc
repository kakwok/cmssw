#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/HcalObjects/interface/alpaka/HcalSiPMCharacteristicsPortable.h"
#include "CondFormats/HcalObjects/interface/HcalSiPMCharacteristicsSoA.h"
#include "CondFormats/DataRecord/interface/HcalSiPMCharacteristicsRcd.h"

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ModuleFactory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"



namespace ALPAKA_ACCELERATOR_NAMESPACE {
  class HcalSiPMCharacteristicsESProducer : public ESProducer {
  public:
    HcalSiPMCharacteristicsESProducer(edm::ParameterSet const& iConfig) : ESProducer(iConfig) {
      auto cc = setWhatProduced(this);
      sipmCharacteristicsToken_ = cc.consumes();
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      descriptions.addWithDefaultLabel(desc);
    }

    std::unique_ptr<HcalSiPMCharacteristicsPortableHost> produce(HcalSiPMCharacteristicsRcd const& iRecord) {

      auto const& sipmCharacteristics = iRecord.get(sipmCharacteristicsToken_);

      size_t const totalItems = sipmCharacteristics.getTypes();

      auto product = std::make_unique<HcalSiPMCharacteristicsPortableHost>(totalItems,cms::alpakatools::host());

      auto view = product->view();

      for (uint32_t i = 0; i < sipmCharacteristics.getTypes(); i++) {

        auto vi = view[i];
        auto const type      = sipmCharacteristics.getType(i);
    
        vi.precisionItem() = HcalSiPMCharacteristics::PrecisionItem(type, 
                                                    sipmCharacteristics.getPixels(type),
                                                    sipmCharacteristics.getNonLinearities(type)[0], 
                                                    sipmCharacteristics.getNonLinearities(type)[1], 
                                                    sipmCharacteristics.getNonLinearities(type)[2], 
                                                    sipmCharacteristics.getCrossTalk(i),  
                                                    sipmCharacteristics.getAuxi1(i), 
                                                    sipmCharacteristics.getAuxi2(i)); 

      }
      return product;
    }

  private:
    edm::ESGetToken<HcalSiPMCharacteristics , HcalSiPMCharacteristicsRcd > sipmCharacteristicsToken_;

  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_EVENTSETUP_ALPAKA_MODULE(HcalSiPMCharacteristicsESProducer);
