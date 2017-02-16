#include "detector_const.hh"
#include "analysis_const.hh"

#ifndef __LEVREC0__
#define __LEVREC0__ 1

class LEvRec0 {

public:
  LEvRec0() {};

  short strip[NCHAN]={0};
  unsigned int trigger_index=0;
  unsigned int hepd_time=0;
  unsigned int event_index=0;
  unsigned short event_length=0;
  unsigned short pmt_high[NPMT]={0};
  unsigned short pmt_low[NPMT]={0};
  unsigned short rate_meter[NRATEMETER]={0};
  bool trigger_flag[NPMT]={0};
  unsigned int alive_time=0;
  unsigned int dead_time=0;

  void DumpStrip(void);
  void DumpEventIndex();
};


#endif
