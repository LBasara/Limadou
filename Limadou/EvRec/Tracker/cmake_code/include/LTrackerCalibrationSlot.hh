#ifndef __LTRACKERCALIBRATIONSLOT__
#define __LTRACKERCALIBRATIONSLOT__ 1

#include "LTrackerMask.hh"
#include <fstream>


class LTrackerCalibrationSlot {
public:
  LTrackerCalibrationSlot(int StartE, int StopE, double *ped, double *sig, double *ngi, bool *cnm);
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
  int StartEvent=0;
  int StopEvent=0;
  double pedestal[NCHAN]={0};
  double sigma[NCHAN]={0};
  double ngindex[NCHAN]={0};
  bool CN_mask[NCHAN]={0};
};




#endif
