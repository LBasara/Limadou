#ifndef LTRACKERQUICKLOOK 
#define LTRACKERQUICKLOOK 1

#include "detector_const.hh"

#include "TCanvas.h"
#include "TH2.h"

const double HOTCHANNELTHRESHOLD=20.;
const double COLDCHANNELTHRESHOLD=3.;
const double GAUSCHANNELTHRESHOLD=3.5;

const int NCALIBEVENTS_QL=1000;

TCanvas *drawing6_chan_2D(TH2D *ladder[N_LADDER]);
void QuickLook(std::string namefile);

#endif
