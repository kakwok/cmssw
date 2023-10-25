#ifndef CondFormats_HcalObjects_HcalMahiConditionsSoA_h
#define CondFormats_HcalObjects_HcalMahiConditionsSoA_h

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalConstants.h"

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"


using  HcalPedestalArray   = std::array<float, 4>; // 4 capIds 

GENERATE_SOA_LAYOUT(HcalMahiConditionsSoALayout,
                    SOA_COLUMN(uint32_t, param1),
                    SOA_COLUMN(uint32_t, param2),
                    SOA_COLUMN(HcalPedestalArray, pedestals_value),
                    SOA_COLUMN(HcalPedestalArray, pedestals_width),
                    SOA_COLUMN(HcalPedestalArray, gains_value),
                    SOA_COLUMN(HcalPedestalArray, convertedPedestals),
                    SOA_COLUMN(HcalPedestalArray, convertedPedestalWidths),
                    SOA_COLUMN(float, lutCorrs_values),
                    SOA_COLUMN(float, respCorrs_values),
                    SOA_COLUMN(float, timeCorrs_values),
                    // Use EIGEN_COLUMN for matrix?
                    SOA_COLUMN(float, pedestalWidths_sigma00),
                    SOA_COLUMN(float, pedestalWidths_sigma01),
                    SOA_COLUMN(float, pedestalWidths_sigma02),
                    SOA_COLUMN(float, pedestalWidths_sigma03),
                    SOA_COLUMN(float, pedestalWidths_sigma10),
                    SOA_COLUMN(float, pedestalWidths_sigma11),
                    SOA_COLUMN(float, pedestalWidths_sigma12),
                    SOA_COLUMN(float, pedestalWidths_sigma13),
                    SOA_COLUMN(float, pedestalWidths_sigma20),
                    SOA_COLUMN(float, pedestalWidths_sigma21),
                    SOA_COLUMN(float, pedestalWidths_sigma22),
                    SOA_COLUMN(float, pedestalWidths_sigma23),
                    SOA_COLUMN(float, pedestalWidths_sigma30),
                    SOA_COLUMN(float, pedestalWidths_sigma31),
                    SOA_COLUMN(float, pedestalWidths_sigma32),
                    SOA_COLUMN(float, pedestalWidths_sigma33),
                    SOA_COLUMN(float, gainWidths_value0),
                    SOA_COLUMN(float, gainWidths_value1),
                    SOA_COLUMN(float, gainWidths_value2),
                    SOA_COLUMN(float, gainWidths_value3),
                    SOA_COLUMN(uint32_t, channelQuality_status),
                    SOA_COLUMN(int, qieTypes_values),
                    SOA_COLUMN(int, sipmPar_type),
                    SOA_COLUMN(int, sipmPar_auxi1),
                    SOA_COLUMN(float, sipmPar_fcByPE),
                    SOA_COLUMN(float, sipmPar_darkCurrent),
                    SOA_COLUMN(float, sipmPar_auxi2)
                    )
using HcalMahiConditionsSoA = HcalMahiConditionsSoALayout<>;

//using  HcalPSfunctorArray   = std::array<float, hcal::constants::maxPSshapeBin>; // 256
//using  HcalPSfunctorBXarray = std::array<float, hcal::constants::nsPerBX>;       // 25
//
//GENERATE_SOA_LAYOUT(HcalRecoParamSoALayout,
//                    SOA_COLUMN(uint32_t, param1),
//                    SOA_COLUMN(uint32_t, param2),
//                    SOA_COLUMN(uint32_t, ids)
//                    )
//GENERATE_SOA_LAYOUT(HcalPulseShapeSoALayout,
//                    SOA_COLUMN(HcalPSfunctorArray, acc25nsVec),
//                    SOA_COLUMN(HcalPSfunctorArray, diff25nsItvlVec),
//                    SOA_COLUMN(HcalPSfunctorBXarray, accVarLenIdxMinusOneVec),
//                    SOA_COLUMN(HcalPSfunctorBXarray, diffVarItvlIdxMinusOneVec),
//                    SOA_COLUMN(HcalPSfunctorBXarray, accVarLenIdxZEROVec),
//                    SOA_COLUMN(HcalPSfunctorBXarray, diffVarItvlIdxZEROVec)
//                    )
//
//using HcalRecoParamSoA      = HcalRecoParamSoALayout<>;
//using HcalPulseShapeSoA     = HcalPulseShapeSoALayout<>;
//
//struct HcalRecoParamsWithPulseShapes{
//    HcalRecoParamSoA     const * recoParam;
//    HcalPulseShapeSoA    const * shape;
//
//}
//


#endif
