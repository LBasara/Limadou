#ifndef __LTRACKERCALIBRATIONMANAGER__
#define __LTRACKERCALIBRATIONMANAGER__ 1

#include "LTrackerCalibration.hh"

const double CHANCLEANINGTHRESHOLD=3.;
const double MINSIGMA1=0.;
const double MAXSIGMA1=50.;
const int NSIGMA1BIN=100;
const double HALFSIGMA1WIDTH=1.5; // ADCs

const double GAUSSIANITYSIGMATHRESHOLD=3.0;
const double GAUSSIANITYEVRACTHRESHOLD=0.003;


class LTrackerCalibrationManager {

public:
  static LTrackerCalibrationManager& GetInstance();
 
  int LoadRun(char *fileInp);
  LTrackerCalibration* Calibrate(int nEvents=-1, int skipEvents=-1);
  int SaveCalibration(char *fileOut);

  LTrackerCalibration* LoadCalibration(char *fileInp); // for future: Load tracker calibration from a total calibration
  LTrackerCalibration* LoadTrackerCalibration(char *fileInp);
  

private:
  char *calRunFileName;
  LTrackerCalibrationManager();
  // CalibrationSlots
  int CalculateCalibrationSlots(int nEvents, int skipEvents, int nEntries,
			    int *pivot);
  LTrackerCalibrationSlot* CalibrateSlot(int StartEntry, int StopEntry);
  void ComputeCNMask(double *sigma1, double &CN_mask) {
  void ComputeCN(double *counts, double *pedestal, bool *CN_mask, double &CN);

  



  // C++ 11
  // =======
  // We can use the better technique of deleting the methods
  // we don't want.
public:
  LTrackerCalibrationManager(LTrackerCalibrationManager const&)               = delete;
  void operator=(LTrackerCalibrationManager const&)  = delete;

};

#endif
