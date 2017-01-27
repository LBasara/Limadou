#ifndef __LTRACKERTOOLS__
#define __LTRACKERTOOLS__ 1

#include <vector>
#include "LTrackerCluster.hh"

int ChanToLadder(const int nStrip);
int ChanToADC(const int nStrip);
int ChanToLadderADC(const int nStrip);
int ChanToADCVA(const int nStrip);
int ChanToVA(const int nStrip);
int ChanToSide(const int nStrip);  // 0 p - 1 n
int ChanToPlane(const int nStrip); // 0 external - 1 internal
int ChanToLadderPlane (const int nChan); // return 0,1,....11
bool SameLadderPlane(const int Chan1, const int Chan2);
int ChanToLadderChan(const int Chan); // channel inside the ladder
int ChanToSideChan(const int Chan); // channel inside the side
std::vector<LTrackerCluster>* GetClusters(const double* cont, const double *sigma);
void ComputeCN(const short *counts, const double *pedestal, const bool *CN_mask, double *CN);
double* Symmetrise(const int inpsize, const double* inp);


#endif
