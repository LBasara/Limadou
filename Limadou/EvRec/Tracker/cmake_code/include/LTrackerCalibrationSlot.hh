#ifndef __LTRACKERCALIBRATIONSLOT__
#define __LTRACKERCALIBRATIONSLOT__ 1

#include "detector_const.hh"

class LTrackerCalibrationSlot {
public:
  LTrackerCalibrationSlot(int StartE, int StopE, double *ped, double *sig, double *ngi);;
  
private:
  // Calib infos
  int StartEvent;
  int StopEvent;
  double pedestal[NCHAN];
  double sigma[NCHAN];
  double ngindex[NCHAN];
};




#endif
