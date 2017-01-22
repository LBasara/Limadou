#ifndef __LTRACKERCALIBRATIONMANAGER__
#define __LTRACKERCALIBRATIONMANAGER__ 1

#include "LTrackerCalibration.hh"
#include "LEvRec0File.hh"

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
  
  int LoadRun(const char *fileInp);
  LTrackerCalibration* Calibrate(const int nEvents=-1, const int skipEvents=-1);
  int SaveCalibration(const char *fileOut);
  
  LTrackerCalibration* LoadCalibration(const char *fileInp) {return 0;}; // for future: Load tracker calibration from a total calibration
  inline LTrackerCalibration* LoadTrackerCalibration(const char *fileInp){return 0;};
  
  
private:
  LEvRec0File *calRunFile;
  LTrackerCalibrationManager();
  
  // CalibrationSlots
  int CalculateCalibrationSlots(const int nEvents, const int skipEvents, const int nEntries,
			    int *pivot);
  LTrackerCalibrationSlot* CalibrateSlot(int StartEntry, const int StopEntry);
  void RawMeanSigma(const int StartEntry, const int StopEntry, double *mean0, double *sigma0);
  void CleanedMeanSigma(const int StartEntry, const int StopEntry, const double *mean0, const double *sigma0, double *mean1, double *sigma1);
  void ComputeCNMask(const double *sigma1, bool *CN_mask);
  void ComputeCN(const short *counts, const double *pedestal, const bool *CN_mask, double *CN);
  void CNCorrectedSigma(const int StartEntry, const int StopEntry, const double *mean1, const double *sigma1, const bool *CN_mask, double *mean2, double *sigma2);
  void GaussianityIndex(const int StartEntry, const int StopEntry, const double *mean2, const double *sigma2, const bool *CN_mask, double *ngindex);

  LTrackerCalibration* CreateTrackerCalibration();
  
  ~LTrackerCalibrationManager();
  



  // C++ 11
  // =======
  // We can use the better technique of deleting the methods
  // we don't want.
public:
  LTrackerCalibrationManager(LTrackerCalibrationManager const&) = delete;
  void operator=(LTrackerCalibrationManager const&) = delete;

};

#endif
