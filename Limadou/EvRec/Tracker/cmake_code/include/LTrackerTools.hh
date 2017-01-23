#ifndef __LTRACKERTOOLS__
#define __LTRACKERTOOLS__ 1

#include <vector>
#include "LTrackerCluster.hh"

int ChanToLadder(int nStrip);
int ChanToADC(int nStrip);
int ChanToLadderADC(int nStrip);
int ChanToADCVA(int nStrip);
int ChanToVA(int nStrip);
int ChanToSide(int nStrip);  // 0 p - 1 n
int ChanToPlane(int nStrip); // 0 external - 1 internal
int ChanToLadderPlane (int nChan); // return 0,1,....11
bool SameLadderPlane(int Chan1, int Chan2);
int ChanToLadderChan(int Chan); // channel inside the ladder
std::vector<LTrackerCluster>* GetClusters(double* cont, double *sigma);
void ComputeCN(const short *counts, const double *pedestal, const bool *CN_mask, double *CN);


#endif
