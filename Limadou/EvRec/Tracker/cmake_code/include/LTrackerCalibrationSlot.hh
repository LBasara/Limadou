#ifndef __LTRACKERCALIBRATIONSLOT__
#define __LTRACKERCALIBRATIONSLOT__ 1

#include "LTrackerMask.hh"
#include "detector_const.hh"
#include <fstream>

class LTrackerCalibrationSlot {
public:
  LTrackerCalibrationSlot(int StartE, int StopE, double *ped, double *sig, double *ngi, bool *cnm);;
  void Write(std::ofstream *output);
  static LTrackerCalibrationSlot* Read(std::ifstream *input);
  inline int GetStartEvent(){return StartEvent;};
  inline int GetStopEvent(){return StopEvent;};
  inline double* GetPedestal(){return pedestal;};
  inline double* GetSigma(){return sigma;};
  inline double* GetNGIndex(){return ngindex;};
  inline bool* GetCNMask(){return CN_mask;};
  LTrackerMask GetMaskOnSigma(const double sigmaMin, const double sigmaMax);
  LTrackerMask GetMaskOnNGI(const double NGIMin, const double NGIMax);
  
private:
  // Calib infos
  int StartEvent;
  int StopEvent;
  double pedestal[NCHAN];
  double sigma[NCHAN];
  double ngindex[NCHAN];
  bool CN_mask[NCHAN];
};




#endif
