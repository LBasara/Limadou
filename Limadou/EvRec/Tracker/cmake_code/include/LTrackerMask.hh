#ifndef __LTRACKERMASK__
#define __LTRACKERMASK__ 1

#include "detector_const.hh"
#include <cstddef>
#include <vector>

class LTrackerMask {
public:
  LTrackerMask();
  LTrackerMask(const bool *mIN);
  LTrackerMask(std::vector<bool> mIN);
  bool* GetBool(void);
  inline bool& operator[](std::size_t idx) {return m[idx];};
  inline bool Get(std::size_t idx) const {return m[idx];};
  LTrackerMask& operator=(const LTrackerMask& other);
  friend LTrackerMask operator&&(const LTrackerMask& m1, const LTrackerMask& m2);
  friend LTrackerMask operator||(const LTrackerMask& m1, const LTrackerMask& m2);
  friend LTrackerMask operator!(const LTrackerMask& m1);
  void Dump();
private:
  bool m[NCHAN];
};
#endif
