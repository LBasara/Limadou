#include "LTrackerMask.hh"
#include <iostream>

LTrackerMask::LTrackerMask() {
  for(auto el : m) el=true;
}

LTrackerMask::LTrackerMask(const bool *mIN) {
  for(int iChan=0; iChan<NCHAN; ++iChan) m[iChan]=mIN[iChan];
}

bool* LTrackerMask::GetBool(void) {
  bool *result = new bool[NCHAN];
  for(int iChan=0; iChan<NCHAN; ++iChan) result[iChan]=m[iChan];
  return result;
}

LTrackerMask& LTrackerMask::operator=(const LTrackerMask& other) {
  for(int iChan=0; iChan<NCHAN; ++iChan) m[iChan]=other.Get(iChan);
  return *this;
}

LTrackerMask operator&&(const LTrackerMask& m1, const LTrackerMask& m2) {
  bool mresult[NCHAN];
  for(int iChan=0; iChan<NCHAN; ++iChan) mresult[iChan] = m1.Get(iChan)&&m2.Get(iChan);
  LTrackerMask result(mresult);
  return result;
}


LTrackerMask operator||(const LTrackerMask& m1, const LTrackerMask& m2) { 
  bool mresult[NCHAN];
  for(int iChan=0; iChan<NCHAN; ++iChan) mresult[iChan] = m1.Get(iChan)||m2.Get(iChan);
  LTrackerMask result(mresult);
  return result;
}

LTrackerMask operator!(const LTrackerMask& m1) {
  bool mresult[NCHAN];
  for(int iChan=0; iChan<NCHAN; ++iChan) mresult[iChan] = !m1.Get(iChan);
  LTrackerMask result(mresult);
  return result;
} 


void LTrackerMask::Dump() {
  for(auto ib : m) std::cout << ib << " ";
  std::cout << std::endl;
  return;
}
