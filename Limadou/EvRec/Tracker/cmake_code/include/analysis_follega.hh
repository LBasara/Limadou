#ifndef __ANALYSIS_FOLLEGA__
#define __ANALYSIS_FOLLEGA__ 1
#include "TH1.h"

using namespace std;

const int n_chann = 4608; //number of channels 
const int n_steps = 1; // number of steps of 1000 events
const int n_ev = 1000; //number of events per step
const int VA_chan = 64; // number of channel in a VA
const int n_VA = 72; // total number of VA
//const int ADC_CHAN = 3*VA_chan;
//const int LADDER_CHAN = 4*ADC_CHAN;
const double sigmarange = 1.5;
//const double MIN_NCHANSPERBIN_4CN=5;
//const double MIN_SIGMA_CHANNEL_OFF=6;
//const double MAX_SIGMA_CHANNEL_OFF=30;
//const double MAX_SIGMA_NOISE_LEVEL=4; 
const double sigma_cut=6.5;
const double gaussian_threshold = 3.;
const double sigma_cut_gauss = 3.;

float GetMean( int * channel, float inf, float sup);
float GetSigma(int * channel, float * noise, float mean, float inf, float sup);
void CommonNoise( int * channel, float *mean, float *sigma, float *max_sigma, float *noise, int * GEcounter);
double GetCleanedSigma(TH1F *h);
//void N_seed(int * ADC, float * mean_2, float * comnoise_real2, double * sigma_3, float significance, int *count);
void Silicon_analysis_cls5(string root_data_file, string typeofparticle, string file_calib, string file_fit_par);
//void Series();


#endif
