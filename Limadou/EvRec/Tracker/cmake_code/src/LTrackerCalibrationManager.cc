#include "LTrackerCalibrationManager.hh"
#include "LEvRec0File.hh"
#include "LEvRec0.hh"


LTrackerCalibrationManager::LTrackerCalibrationManager() {
  calRunFileName=0;
}


LTrackerCalibrationManager& LTrackerCalibrationManager::GetInstance() {
  static LTrackerCalibrationManager instance; // Guaranteed to be destroyed.
                                              // Instantiated on first use.
  return instance;
}


int LTrackerCalibrationManager::LoadRun(char *fileInp) {
  calRunFileName = fileInp;
  return 0;
}

LTrackerCalibration* LTrackerCalibrationManager::Calibrate(int nEvents) {
  
  // Open the run file
  LEvRec0File calRunFile(calRunFileName);
  LEvRec0 calEv;
  calRunFile.SetTheEventPointer(calEv);



  LTrackerCalibration *result = new LTrackerCalibration();
  return result;
}
