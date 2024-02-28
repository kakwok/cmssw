
#include<vector>

#include "DataFormats/HcalDigi/interface/alpaka/HcalDigiDeviceCollection.h"
#include "DataFormats/HcalRecHit/interface/alpaka/HcalRecHitDeviceCollection.h"
#include "CondFormats/HcalObjects/interface/alpaka/HcalMahiConditionsDevice.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "DeclsForKernels.h"


namespace ALPAKA_ACCELERATOR_NAMESPACE::hcal::reconstruction {


    using IProductTypef01 = hcal::Phase1DigiDeviceCollection;   
    using IProductTypef5 = hcal::Phase0DigiDeviceCollection; 
    using IProductTypef3 = hcal::Phase1DigiDeviceCollection; 
    using OProductType = hcal::RecHitDeviceCollection; 

    void entryPoint(Queue& queue,
                    IProductTypef01 const& f01HEDigis,
                    IProductTypef5 const& f5HBDigis,
                    IProductTypef3 const& f3HBDigis,
                    OProductType& outputGPU,
                    ConditionsProducts const& conditions,
                    ConfigParameters const& configParameters
                    );

} //ALPAKA_ACCELERATOR_NAMESPACE::hcal::reconstruction
