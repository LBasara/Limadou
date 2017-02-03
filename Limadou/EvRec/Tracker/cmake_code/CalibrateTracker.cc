#include <iostream>
#include "LTrackerCalibrationManager.hh"
#include "LTrackerMask.hh"

int main(int argc, char *argv[]) {
  if(argc!=3) {
    std::cerr << "Error! Usage:    ./CalibrateTracker <calRunFile> <calOutFile>" << std::endl;
    std::cerr << "Aborted." << std::endl;
    return -999;
  }
  
  
  LTrackerCalibrationManager::GetInstance().LoadRun(argv[1]);
  LTrackerCalibration *cal =   LTrackerCalibrationManager::GetInstance().Calibrate();
  cal->Write(argv[2]);

  /*
  // Test read-write
  LTrackerCalibration *cal=LTrackerCalibration::Read(argv[2]);
  cal->Write("ciccio.cal");
  */

  
  // Test LTrackerMask
  auto ms = cal->GetMaskOnSigma(0, 4.,6.);
  ms.Dump();
  auto mngi = cal->GetMaskOnNGI(0, -100.,5.);
  mngi.Dump();
  auto mtot = ms&&mngi;
  mtot.Dump();
  
  return 0;
}
