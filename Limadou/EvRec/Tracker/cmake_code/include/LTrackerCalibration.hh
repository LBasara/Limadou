#ifndef __LTRACKERCALIBRATION__
#define __LTRACKERCALIBRATION__ 1

#include "LTrackerCalibrationSlot.hh"
#include <vector>


class LTrackerCalibration {
public:
  LTrackerCalibration(const int RunIdINP, const int InitialTargetRunINP, const int FinalTargetRun);
  void Add(const LTrackerCalibrationSlot *lcal);
  void Write(const char *fileOut);
  static LTrackerCalibration* Read(const char *fileIn);
  inline int GetNSlots(){return nSlots;};
  inline int GetRunId(){return RunId;};
  inline int GetInitialTargetRun(){return InitialTargetRun;};
  inline int GetFinalTargetRun(){return FinalTargetRun;};
  inline double* GetPedestal(const int nSlot){return calarray.at(nSlot).GetPedestal();};
  inline double* GetSigma(const int nSlot){return calarray.at(nSlot).GetSigma();};
  inline double* GetNGIndex(const int nSlot){return calarray.at(nSlot).GetNGIndex();};
  inline bool* GetCNMask(const int nSlot){return calarray.at(nSlot).GetCNMask();};
  inline LTrackerMask GetMaskOnSigma(const int nSlot, const double sigmaMin, const double sigmaMax){return calarray.at(nSlot).GetMaskOnSigma(sigmaMin, sigmaMax);};
  inline LTrackerMask GetMaskOnNGI(const int nSlot, const double ngiMin, const double ngiMax){return calarray.at(nSlot).GetMaskOnNGI(ngiMin, ngiMax);};
  
private:
  // Calib infos
  int RunId;
  int nSlots;
  int InitialTargetRun;
  int FinalTargetRun;
  std::vector<LTrackerCalibrationSlot> calarray;

  LTrackerCalibration();
  
};


#endif

