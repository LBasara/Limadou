#ifndef __ESTER_CAL2__
#define __ESTER_CAL2__ 1

#include "TCanvas.h"
#include "TH2.h"

#include "detector_const.hh"


TCanvas *drawing6_2D(TH2D *ladder[N_LADDER]);
void analysis(std::string namefile,std::string calib_file,std::string outputname);
void analysis_ondata(std::string inputfile,std::string output_name); //run on a single data set


void run1();//test: analysis of all the data from test beam



#endif
