#include <iostream>
#include "LTrackerCalibrationManager.hh"

int main(int argc, char *argv[]) {
  if(argc!=2) {
    std::cerr << "Error! Usage:    ./CalibrateTracker <calRunFile> <calOutFile>" << std::endl;
    std::cerr << "Aborted." << std::endl;
  }
  
  
  LTrackerCalibrationManager::GetInstance().LoadRun(argv[1]);
  LTrackerCalibration *cal =   LTrackerCalibrationManager::GetInstance().Calibrate();
  cal->Write(argv[2]);
  
  return 0;
}
