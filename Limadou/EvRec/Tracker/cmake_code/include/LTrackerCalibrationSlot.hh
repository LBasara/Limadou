#ifndef __LTRACKERCALIBRATIONSLOT__
#define __LTRACKERCALIBRATIONSLOT__ 1

#include "detector_const.hh"
#include <fstream>

class LTrackerCalibrationSlot {
public:
  LTrackerCalibrationSlot(int StartE, int StopE, double *ped, double *sig, double *ngi);;
  void Write(std::ofstream *output);
  static LTrackerCalibrationSlot* Read(std::ifstream *input);
  inline int GetStartEvent(){return StartEvent;};
  inline int GetStopEvent(){return StopEvent;};
  inline double* GetPedestal(){return pedestal;};
  inline double* GetSigma(){return sigma;};
  inline double* GetNGIndex(){return ngindex;};
  
private:
  // Calib infos
  int StartEvent;
  int StopEvent;
  double pedestal[NCHAN];
  double sigma[NCHAN];
  double ngindex[NCHAN];
};




#endif
