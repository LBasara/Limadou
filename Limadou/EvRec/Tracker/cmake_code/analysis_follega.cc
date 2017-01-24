#include "LEvRec0.hh"
#include "analysis_const.hh"
#include "detector_const.hh"
#include "LTrackerTools.hh"
#include "LTrackerCluster.hh"
#include "analysis_follega.hh"
#include "LPlotTools.hh"
#include "LTrackerFittingTools.hh"
#include "LTrackerCalibrationManager.hh"
#include "LTrackerMask.hh"
#include <iostream>

struct event {
  int entry;
  std::vector<LTrackerCluster> cls;
};


				//==================================================================================================
				//========================================== MAIN =================================================
				//==================================================================================================


int main(int argc, char *run_calib[]){
	//calibration procedure 
	if(argc!=3) {
    	std::cerr << "Error! Usage:    ./analysis_follega <calRunFile> <calOutFile>" << std::endl;
    	std::cerr << "Aborted." << std::endl;
    return -999;
  	}
	LTrackerCalibrationManager::GetInstance().LoadRun(run_calib[1]);
	LTrackerCalibration *cal =   LTrackerCalibrationManager::GetInstance().Calibrate();
  	cal->Write(run_calib[2]);

  	double MEAN_CALIB[NCHAN];
  	double SIGMA_CALIB[NCHAN];
  	double NOGAUSS_CALIB[NCHAN];


}






		




