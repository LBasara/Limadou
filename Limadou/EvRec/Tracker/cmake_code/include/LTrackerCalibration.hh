#ifndef __LTRACKERCALIBRATION__
#define __LTRACKERCALIBRATION__ 1

class LTrackerCalibration {
public:
  inline LTrackerCalibration(){;};
  
private:
  // Calib infos
  int RunNumber;
  int InitialTargetRun;
  int FinalTargetRun;
  int InitialTargetEvent;
  int FinalTargetEvent;

};

#endif

