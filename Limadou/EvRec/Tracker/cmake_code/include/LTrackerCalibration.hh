#ifndef __LTRACKERCALIBRATION__
#define __LTRACKERCALIBRATION__ 1

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

