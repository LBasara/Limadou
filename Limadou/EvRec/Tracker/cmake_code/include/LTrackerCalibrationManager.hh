#ifndef __LTRACKERCALIBRATIONMANAGER__
#define __LTRACKERCALIBRATIONMANAGER__ 1

#include "LTrackerCalibration.hh"

class LTrackerCalibrationManager {

public:
  static LTrackerCalibrationManager& GetInstance();
 
  int LoadRun(char *fileInp);
  LTrackerCalibration* Calibrate(int nEvents=-1);
  int SaveCalibration(char *fileOut);

  LTrackerCalibration* LoadCalibration(char *fileInp); // for future: Load tracker calibration from a total calibration
  LTrackerCalibration* LoadTrackerCalibration(char *fileInp);
  

private:
  char *calRunFileName;
  
  LTrackerCalibrationManager();
  
  // C++ 11
  // =======
  // We can use the better technique of deleting the methods
  // we don't want.
public:
  LTrackerCalibrationManager(LTrackerCalibrationManager const&)               = delete;
  void operator=(LTrackerCalibrationManager const&)  = delete;

};

#endif
