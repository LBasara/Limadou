#include "LTrackerCalibration.hh"

void LTrackerCalibration::Add(LTrackerCalibrationSlot *lcal) {
  calarray.push_back(*lcal);
  return;
}


