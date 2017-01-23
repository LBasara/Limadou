#include "detector_const.hh"
#include "analysis_const.hh"

#ifndef __LEVREC0__
#define __LEVREC0__ 1

class LEvRec0 {

public:
  LEvRec0();
  
  short strip[NCHAN];
  unsigned int trigger_index;
  unsigned int hepd_time;
  unsigned int event_index;
  unsigned short event_length;
  unsigned short pmt_high[NPMT];
  unsigned short pmt_low[NPMT];
  unsigned short rate_meter[NRATEMETER];
  bool trigger_flag[NPMT];
  unsigned int alive_time;
  unsigned int dead_time;

  void DumpStrip(void);
  void DumpEventIndex();
};


#endif
