#include "LTrackerTools.hh"
#include "detector_const.hh"
#include "LTrackerCluster.hh"
#include <vector>
#include <algorithm>
#include <math.h>


int ChanToLadder(int nStrip) {
  if(nStrip<0 || nStrip>=NCHAN)
    return -1;
  else 
    return nStrip/LADDER_CHAN;
}

int ChanToADC(int nStrip) {
  int nLadder = ChanToLadder(nStrip);
  int reducednStrip = nStrip - LADDER_CHAN*nLadder;
  return reducednStrip/ADC_CHAN;
}

int ChanToVA(int nStrip) {
  int nADC = ChanToADC(nStrip);
  int reducednStrip = nStrip - ADC_CHAN*nADC;
  return reducednStrip/VA_CHAN;
}

int ChanToSide(int nStrip) { // 0 p - 1 n
  int nADC = ChanToADC(nStrip);
  return (nADC/2)%2;
}

int ChanToPlane(int nStrip) { // 0 external - 1 internal
  int scale=ChanToLadder(nStrip);
  if(scale%2) return 1;
  else return 0;
}

int ChanToLadderPlane (int nChan) { // return 0,1,....11
  return ChanToLadder(nChan)*2+ChanToSide(nChan);
}

bool SameLadderPlane(int Chan1, int Chan2) {
  return (ChanToLadderPlane(Chan1) == ChanToLadderPlane(Chan2));
}

int ChanToLadderChan(int Chan) {
  int result= Chan-LADDER_CHAN*ChanToLadder(Chan);
  return result;
}


std::vector<LTrackerCluster>* GetClusters(double* cont, double *sigma) {
  
  auto result= new std::vector<LTrackerCluster>;
  double sn[NCHAN];
  for(int ich=0; ich<NCHAN; ++ich) sn[ich]=cont[ich]/sigma[ich];
  
  // Main loop
  for(int ich=0; ich<NCHAN; ++ich) {
    if(sn[ich]<CLFINDTHRESHOLD) continue;
 
    // Check if there's a higher maximum  - up to two chans ahead
    int maxindex1=ich;
    double max1 = sn[ich];
    int extentAhead=std::min(LADDER_CHAN-1-ChanToLadderChan(ich),2);
    for(int index=1; index<=extentAhead; ++index) {
      if(sn[ich+index]>max1) {
	maxindex1=ich+index;
	max1=sn[ich+index];
      }
    } // Endl loop for max1

    // Find the highest SN adjacent to the maximum
    int maxindex2;
    int lchmax = ChanToLadderChan(maxindex1);
    if(lchmax==0) maxindex2=maxindex1+1;
    else if(lchmax==LADDER_CHAN-1) maxindex2=maxindex1-1;
    else maxindex2=(sn[maxindex1-1]>sn[maxindex1+1] ? maxindex1-1 : maxindex1+1);

    // Check if the maximum pair is suitable for cluster finding
    // Compute the estimator
    double numerator=(cont[maxindex1]+cont[maxindex2]);
    double denominator=sqrt(sigma[maxindex1]*sigma[maxindex1]+
			    sigma[maxindex2]*sigma[maxindex2]);
    double meter=numerator/denominator;
    // Compare with the threshold
    if(meter<CLSNTHRESHOLD) {
      ich+=2; // already explored up to ich+2
      continue;
    }
    LTrackerCluster mycl(maxindex1, cont, sigma);
    result->push_back(mycl);
    ich=maxindex1+2; // already explored up there.

  } // end of the main loop:: Warning, possible overlap between clusters' boundaries

  return result;
}
