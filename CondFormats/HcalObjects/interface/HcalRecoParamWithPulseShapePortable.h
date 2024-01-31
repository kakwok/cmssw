#ifndef CondFormats_HcalObjects_interface_HcalRecoParamWithPulseShapePortable_h
#define CondFormats_HcalObjects_interface_HcalRecoParamWithPulseShapePortable_h

#include "CondFormats/HcalObjects/interface/HcalRecoParamWithPulseShapeSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "DataFormats/Portable/interface/PortableCollection.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"

//using HcalRecoParamPortableHost = PortableHostCollection<HcalRecoParamSoA>;
//using HcalPulseShapePortableHost = PortableHostCollection<HcalPulseShapeSoA>;
//class HcalRecoParamWithPulseShapeHost{
//
//    HcalRecoParamPortableHost recoParam_;
//    HcalPulseShapePortableHost pulseShape_;
//
//   // ConstView is passed to the kernel by value
//    class ConstView {
//    
//      float pulseShape(int id) const {
//        // code to go through the indirection from view1_ to view2_/view3_
//      }
//    
//      SOA1ConstView view1_;
//      SOA2ConstView view2_;
//      SOA3ConstView view3_;
//    }; 
//}

template <typename TDev>
class HcalRecoParamWithPulseShapeT {
public:
  using RecoParamCollection = PortableCollection<HcalRecoParamSoA, TDev>;
  using PulseShapeCollection = PortableCollection<HcalPulseShapeSoA, TDev>;

  class ConstView {
    constexpr float pulseShape(int id) const {return 0; }

  private:
    typename RecoParamCollection::View recoParamView_;
    typename PulseShapeCollection::View pulseShapeView_;
  };

  HcalRecoParamWithPulseShapeT(size_t recoSize, size_t pulseSize, TDev const& dev) : recoParam_(recoSize, dev), pulseShape_(pulseSize, dev) {}
  template <typename TQueue, typename = std::enable_if_t<alpaka::isQueue<TQueue>>>
    HcalRecoParamWithPulseShapeT(size_t recoSize, size_t pulseSize, TQueue const& queue) : recoParam_(recoSize, queue), pulseShape_(pulseSize, queue) {}
  HcalRecoParamWithPulseShapeT(RecoParamCollection reco, PulseShapeCollection pulse) : recoParam_(std::move(reco)), pulseShape_(std::move(pulse)) {}

  const RecoParamCollection& recoParam() const { return recoParam_; }
  const PulseShapeCollection& pulseShape() const { return pulseShape_; }
  
  typename RecoParamCollection::View recoParamView() { return recoParam_.view(); }
  typename PulseShapeCollection::View pulseShapeView() { return pulseShape_.view(); }

private:
  RecoParamCollection recoParam_;
  PulseShapeCollection pulseShape_;
};

using HcalRecoParamWithPulseShapeHost = HcalRecoParamWithPulseShapeT<alpaka::DevCpu>;


#endif
