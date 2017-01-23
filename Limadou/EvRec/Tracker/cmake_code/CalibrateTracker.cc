#include <iostream>
#include "LTrackerCalibrationManager.hh"

int main(int argc, char *argv[]) {
  if(argc!=3) {
    std::cerr << "Error! Usage:    ./CalibrateTracker <calRunFile> <calOutFile>" << std::endl;
    std::cerr << "Aborted." << std::endl;
    return -999;
  }
  
  
  LTrackerCalibrationManager::GetInstance().LoadRun(argv[1]);
  LTrackerCalibration *cal =   LTrackerCalibrationManager::GetInstance().Calibrate();
  cal->Write(argv[2]);
  //  cal->Read(argv[2]);
  //cal->Write("ciccio.cal");
  
  return 0;
}
