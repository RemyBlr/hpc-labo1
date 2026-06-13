//
// Created by Rémy bleuer on 13.06.2026.
//

#ifndef ECG_DEALINATION_ECG_STREAMING_H
#define ECG_DEALINATION_ECG_STREAMING_H

#pragma once
#include "ecg_processing.h"
#include "output_structs.h"

int ecg_analyze_streaming(const ECG_Params *params,
                          const double *signal,
                          int n_samples,
                          ECG_Peaks *peaks,
                          ECG_Intervals *intervals);

#endif //ECG_DEALINATION_ECG_STREAMING_H
