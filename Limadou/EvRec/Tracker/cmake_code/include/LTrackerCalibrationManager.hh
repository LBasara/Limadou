#ifndef __LTRACKERCALIBRATIONMANAGER__
#define __LTRACKERCALIBRATIONMANAGER__ 1

#include "LTrackerCalibration.hh"
#include "LEvRec0File.hh"
#include "statq.hh"

const double CHANCLEANINGTHRESHOLD=3.;
const double MINSIGMA1=0.;
const double MAXSIGMA1=50.;
const int NSIGMA1BIN=100;
const double HALFSIGMA1WIDTH=1.5; // ADCs

const double GAUSSIANITYSIGMATHRESHOLD=3.0;
const double GAUSSIANITYEVRACTHRESHOLD=0.003;


class LTrackerCalibrationManager {

public:
  static LTrackerCalibrationManager& GetInstance();

  int LoadRun(const char *fileInp);
  void SetTargetRuns(const int InitialRun, const int FinalRun=-1);
  LTrackerCalibration* Calibrate(const int nEvents=-1, const int skipEvents=-1);

private:
  LEvRec0File *calRunFile;  // pointer to the run used for calibration
  int InitialTargetRun;     // Run id of first target run
  int FinalTargetRun;       // Run id of last target run
  LTrackerCalibrationManager();

  // CalibrationSlots
  std::vector<int>  CalculateCalibrationSlots(const int nEvents, const int skipEvents, const int nEntries);
  LTrackerCalibrationSlot* CalibrateSlot(int StartEntry, const int StopEntry);
  std::vector<statq> RawMeanSigma(const int StartEntry, const int StopEntry);
  std::vector<statq> CleanedMeanSigma(const int StartEntry, const int StopEntry, std::vector<statq> rawstat);
  bool*  ComputeCNMask(std::vector<statq> cleanstat);
  std::vector<statq> CNCorrectedSigma(const int StartEntry, const int StopEntry,  std::vector<statq> statclean, const bool* CN_mask);
  std::vector<double> GaussianityIndex (const int StartEntry, const int StopEntry, std::vector<statq> statCNcorr, const bool* CN_mask);

  LTrackerCalibration* CreateTrackerCalibration();

  ~LTrackerCalibrationManager();

  bool verboseFLAG;

  /*
  // C++ 03
  // ========
  // Dont forget to declare these two. You want to make sure they
  // are unacceptable otherwise you may accidentally get copies of
  // your singleton appearing.
  LTrackerCalibrationManager(LTrackerCalibrationManager const&);              // Don't Implement
  void operator=(LTrackerCalibrationManager const&); // Don't implement
  */

  //  /* Following implementation to bepreferred, Not yet fully compatible

  // C++ 11
  // =======
  // We can use the better technique of deleting the methods
  // we don't want.
public:
  LTrackerCalibrationManager(LTrackerCalibrationManager const&) = delete;
  void operator=(LTrackerCalibrationManager const&) = delete;
  // */
};

#endif
