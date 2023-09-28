#ifndef CondFormats_HcalObjects_HcalMahiConditionsSoA_h
#define CondFormats_HcalObjects_HcalMahiConditionsSoA_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"


GENERATE_SOA_LAYOUT(HcalMahiConditionsSoALayout,
                    SOA_COLUMN(uint32_t, param1),
                    SOA_COLUMN(uint32_t, param2),
                    SOA_COLUMN(float, pedestals_values),
                    SOA_COLUMN(float, pedestals_widths),
                    SOA_COLUMN(float, gains_values),
                    SOA_COLUMN(float, lutcCorrs_values),
                    SOA_COLUMN(float, respCorrs_values),
                    SOA_COLUMN(float, timeCorrs_values),
                    // Use EIGEN_COLUMN for matrix?
                    SOA_COLUMN(float, pedestalsWidths_sigma00),
                    SOA_COLUMN(float, pedestalsWidths_sigma01),
                    SOA_COLUMN(float, pedestalsWidths_sigma02),
                    SOA_COLUMN(float, pedestalsWidths_sigma03),
                    SOA_COLUMN(float, pedestalsWidths_sigma10),
                    SOA_COLUMN(float, pedestalsWidths_sigma11),
                    SOA_COLUMN(float, pedestalsWidths_sigma12),
                    SOA_COLUMN(float, pedestalsWidths_sigma13),
                    SOA_COLUMN(float, pedestalsWidths_sigma20),
                    SOA_COLUMN(float, pedestalsWidths_sigma21),
                    SOA_COLUMN(float, pedestalsWidths_sigma22),
                    SOA_COLUMN(float, pedestalsWidths_sigma23),
                    SOA_COLUMN(float, pedestalsWidths_sigma30),
                    SOA_COLUMN(float, pedestalsWidths_sigma31),
                    SOA_COLUMN(float, pedestalsWidths_sigma32),
                    SOA_COLUMN(float, pedestalsWidths_sigma33),
                    SOA_COLUMN(float, gainWidths_value0),
                    SOA_COLUMN(float, gainWidths_value1),
                    SOA_COLUMN(float, gainWidths_value2),
                    SOA_COLUMN(float, gainWidths_value3),
                    SOA_COLUMN(uint32_t, channelQuality_status),
                    SOA_COLUMN(int, qieTypes_value),
                    SOA_COLUMN(int, sipmPar_type),
                    SOA_COLUMN(int, sipmPar_auxi1),
                    SOA_COLUMN(float, sipmPar_fcByPE),
                    SOA_COLUMN(float, sipmPar_darkCurrent),
                    SOA_COLUMN(float, sipmPar_auxi2)
                    )

using HcalMahiConditionsSoA = HcalMahiConditionsSoALayout<>;

#endif
