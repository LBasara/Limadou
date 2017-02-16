
#include "LEvRec0.hh"
#include <iostream>



void LEvRec0::DumpStrip(void) {
  std::cout << "strip" << std::endl;
  for(int i=0; i<NCHAN;++i) std::cout << strip[i] << " ";
  std::cout << std::endl;
  return;
}


void LEvRec0::DumpEventIndex() {
  std::cout << "event_index " << event_index << std::endl;
  return;
}
