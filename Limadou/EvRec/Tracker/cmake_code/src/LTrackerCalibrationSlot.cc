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

void LTrackerCalibrationSlot::Write(std::ofstream *output) {
  *output << "EVENTS " << StartEvent << " " << StopEvent << std::endl;
  for(int iChan=0; iChan<NCHAN; ++iChan)
    *output << pedestal[iChan] << " " << sigma[iChan] << " " << ngindex[iChan] << std::endl;
  return;
}

LTrackerCalibrationSlot* LTrackerCalibrationSlot::Read(std::ifstream *input) {
  char word[100];
  int StartEventST, StopEventST;
  double pedestalST[NCHAN], sigmaST[NCHAN], ngindexST[NCHAN];
  *input >> word >> StartEventST >>  StopEventST;
  for(int iChan=0; iChan<NCHAN; ++iChan)
    *input >> pedestalST[iChan] >>  sigmaST[iChan] >> ngindexST[iChan];
  
  LTrackerCalibrationSlot *result = new LTrackerCalibrationSlot(StartEventST, StopEventST, pedestalST, sigmaST, ngindexST);
  return result;
}
