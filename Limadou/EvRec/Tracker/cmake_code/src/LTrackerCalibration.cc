#include "LTrackerCalibration.hh"
#include <fstream>
#include <iostream>

void LTrackerCalibration::Add(const LTrackerCalibrationSlot *lcal) {
  calarray.push_back(*lcal);
  ++nSlots;
  std::cout << "LTrackerCalibration: " << "current slot number " << nSlots << std::endl;;
  return;
}

LTrackerCalibration::LTrackerCalibration() {
  RunId=-99999;
  nSlots=0;
  InitialTargetRun=-99999;
  FinalTargetRun=-99999;
}

LTrackerCalibration::LTrackerCalibration(const int RunIdINP,  const int InitialTargetRunINP, const int FinalTargetRunINP) {
  RunId=RunIdINP;
  nSlots=0;
  InitialTargetRun=InitialTargetRunINP;
  FinalTargetRun=FinalTargetRunINP;
}

void LTrackerCalibration::Write(const char *fileOut) {
  std::ofstream output(fileOut, std::ofstream::out); 
  
  output << RunId << std::endl;
  output << InitialTargetRun << " " << FinalTargetRun << std::endl;
  output << nSlots << std::endl;

  for(auto cslotit : calarray) cslotit.Write(&output);
  
  output << std::endl;
  output.close();
  
  return;
}


LTrackerCalibration* LTrackerCalibration::Read(const char *fileIn) {
  
  std::ifstream input(fileIn, std::ifstream::in); 
  
  int RunIdST, InitialTargetRunST, FinalTargetRunST, nSlotsST;

  input >> RunIdST;
  input >> InitialTargetRunST >> FinalTargetRunST;

  LTrackerCalibration *result = new LTrackerCalibration(RunIdST, InitialTargetRunST, FinalTargetRunST);

  input >> nSlotsST;
  for(int is=0; is<nSlotsST; ++is) result->Add(LTrackerCalibrationSlot::Read(&input));
  
  input.close();
  

  return result;
}
 
