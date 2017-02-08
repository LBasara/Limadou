#ifndef __LTRACKERTOOLS__
#define __LTRACKERTOOLS__ 1

#include <vector>
#include "LTrackerCluster.hh"

int ChanToLadder(const int nStrip);
int ChanToADC(const int nStrip);
int ChanToLadderADC(const int nStrip);
int ChanToADCVA(const int nStrip);
int ChanToVA(const int nStrip);
int ChanToVAChan(const int Chan);
int ChanToSide(const int nStrip);  // 0 p - 1 n
int ChanToPlane(const int nStrip); // 0 external - 1 internal
int ChanToLadderPlane (const int nChan); // return 0,1,....11
bool SameLadderPlane(const int Chan1, const int Chan2);
int ChanToLadderChan(const int Chan); // channel inside the ladder
int ChanToSideChan(const int Chan); // channel inside the side
std::vector<LTrackerCluster>* GetClusters(const double* cont, const double *sigma, const bool *mask);
inline std::vector<LTrackerCluster>* GetClusters(const double* cont, const double *sigma){return GetClusters(cont, sigma, 0);};
void ComputeCN(const short *counts, const double *pedestal, const bool *CN_mask, double *CN);


template<class T> T* Symmetrise(const int inpsize, const T* inp) {
  if(inpsize<0) return 0;
  int CENTER=inpsize/2;
  T *result = new T[inpsize];
  for(int i=0; i<CENTER; ++i) {
    result[i]=(inp[i]+inp[inpsize-1-i])*0.5;
    result[inpsize-1-i]=result[i];
  }
  
  if(inpsize%2) result[CENTER] = inp[CENTER]; // odd case
  
  return result;
} 



#endif
