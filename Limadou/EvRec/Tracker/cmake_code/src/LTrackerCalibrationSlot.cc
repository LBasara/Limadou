#include "LTrackerCalibrationSlot.hh"

LTrackerCalibrationSlot::LTrackerCalibrationSlot(int StartE, int StopE, double *ped, double *sig, double *ngi) {
  StartEvent=StartE;
  StopEvent=StopE;
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    pedestal[iChan]=ped[iChan];
    sigma[iChan]=sig[iChan];
    ngindex[iChan]=ngi[iChan];
  }
}

