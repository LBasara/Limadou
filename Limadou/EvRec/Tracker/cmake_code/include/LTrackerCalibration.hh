#ifndef __LTRACKERCALIBRATION__
#define __LTRACKERCALIBRATION__ 1

#include "LTrackerCalibrationSlot.hh"
#include <vector>


class LTrackerCalibration {
public:
  inline LTrackerCalibration(){;};
  void Add(LTrackerCalibrationSlot *lcal);
  
private:
  // Calib infos
  int RunNumber;
  int InitialTargetRun;
  int FinalTargetRun;
  int InitialTargetEvent;
  int FinalTargetEvent;
  std::vector<LTrackerCalibrationSlot> calarray;
  
};


#endif

