#ifndef __FRANCESCO_ANALYSIS__
#define __FRANCESCO_ANALYSIS__ 1
#include "TH1.h"

float GetMean( int * channel, float inf, float sup);
float GetSigma(int * channel, float * noise, float mean, float inf, float sup);
void CommonNoise( int * channel, float *mean, float *sigma, float *max_sigma, float *noise, int * GEcounter);
double GetCleanedSigma(TH1F *h);
void N_seed(int * ADC, float * mean_2, float * comnoise_real2, float * sigma_3, float significance, int *count);
void Silicon_analysis_cls5(string root_data_file, string typeofparticle, string file_calib, string file_fit_par);
void Series();


#endif
