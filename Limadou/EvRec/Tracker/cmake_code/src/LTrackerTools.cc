#include "LTrackerTools.hh"
#include "detector_const.hh"
#include "LTrackerCluster.hh"
#include "LTrackerMask.hh"
#include <algorithm>
#include <math.h>
#include <iostream>

int ChanToLadder(const int nStrip) {
  if(nStrip<0 || nStrip>=NCHAN)
    return -1;
  else 
    return nStrip/LADDER_CHAN;
}

int ChanToLadderADC(const int nStrip) {
  int nLadder = ChanToLadder(nStrip);
  int reducednStrip = nStrip - LADDER_CHAN*nLadder;
  return reducednStrip/ADC_CHAN;
}

int ChanToADC(const int nStrip) {
  return nStrip/ADC_CHAN;
}

int ChanToADCVA(const int nStrip) {
  int nADC = ChanToADC(nStrip);
  int reducednStrip = nStrip - ADC_CHAN*nADC;
  return reducednStrip/VA_CHAN;
}

int ChanToVA(const int nStrip) {
  return nStrip/VA_CHAN;
}

int ChanToVAChan(const int Chan) {
  return Chan%VA_CHAN;
}

int ChanToSide(const int nStrip) { // 0 p - 1 n
  int nADC = ChanToLadderADC(nStrip);
  return (nADC/2)%2;
}

int ChanToPlane(const int nStrip) { // 0 external - 1 internal
  int scale=ChanToLadder(nStrip);
  if(scale%2) return 1;
  else return 0;
}

int ChanToLadderPlane (const int nChan) { // return 0,1,....11
  return ChanToLadder(nChan)*2+ChanToSide(nChan);
}

bool SameLadderPlane(const int Chan1, const int Chan2) {
  return (ChanToLadderPlane(Chan1) == ChanToLadderPlane(Chan2));
}

int ChanToLadderChan(const int Chan) {
  int result= Chan-LADDER_CHAN*ChanToLadder(Chan);
  return result;
}

int ChanToSideChan(const int Chan) {
  int result= ChanToLadderChan(Chan)%SIDE_CHAN;
  return result;
}


std::vector<LTrackerCluster>* GetClusters(const double* cont, const double *sigma, const bool *maskIn=0) {
  
  auto result= new std::vector<LTrackerCluster>;
  double sn[NCHAN];
  for(int ich=0; ich<NCHAN; ++ich) sn[ich]=cont[ich]/sigma[ich];

  // Prepare mask even for default case
  double mask[NCHAN];
  if(maskIn==0) for(int ich=0; ich<NCHAN; ++ich) mask[ich]=true;
  else for(int ich=0; ich<NCHAN; ++ich) mask[ich]=maskIn[ich];
  // Apply the mask
  for(int ich=0; ich<NCHAN; ++ich) sn[ich]*=(static_cast<double>(mask[ich]));
  
  // Main loop
  for(int ich=0; ich<NCHAN; ++ich) {
    //for(int ich=NCHAN-1; ich>-1; --ich) {
    if(sn[ich]<CLFINDTHRESHOLD) continue;
 
    int maxindex1=ich;
    double max1 = sn[ich];
    // Check if there's a higher maximum  - up to two chans ahead
    int newmaxindex1 = maxindex1;
    double newmax1 = max1;
    while(1) {
      int extentAhead=std::min(SIDE_CHAN-1-ChanToSideChan(maxindex1),2);
      for(int index=1; index<=extentAhead; ++index) {
	if(sn[maxindex1+index]>newmax1) {
	  newmaxindex1=maxindex1+index;
	  newmax1=sn[maxindex1+index];
	}
      } // Endl loop for newmax1
      if(newmaxindex1 == maxindex1) {
	break;
      } else {
	maxindex1 = newmaxindex1;
	max1 = newmax1;
      }
    }

    // Find the highest SN adjacent to the maximum
    int maxindex2;
    int schmax = ChanToSideChan(maxindex1);
    if(schmax==0) maxindex2=maxindex1+1;
    else if(schmax==SIDE_CHAN-1) maxindex2=maxindex1-1;
    else maxindex2=(sn[maxindex1-1]>sn[maxindex1+1] ? maxindex1-1 : maxindex1+1);

    // Check if the maximum pair is suitable for cluster finding
    // Compute the estimator
    double numerator=(cont[maxindex1]*static_cast<double>(mask[maxindex1])+
		      cont[maxindex2]*static_cast<double>(mask[maxindex2]));
    double denominator=sqrt(
			    sigma[maxindex1]*sigma[maxindex1]*static_cast<double>(mask[maxindex1])+
			    sigma[maxindex2]*sigma[maxindex2]*static_cast<double>(mask[maxindex2])
			    );
    double meter=numerator/denominator;
    // Compare with the threshold
    if(meter<CLSNTHRESHOLD) {
      ich+=2; // already explored up to ich+2
      //ich-=2;
      continue;
    }
    LTrackerCluster mycl(maxindex1, cont, sigma);
    result->push_back(mycl);
    ich=maxindex1+2; // already explored up there.
    //ich=maxindex1-2; // already explored up there.
  } // end of the main loop:: Warning, possible overlap between clusters' boundaries

  return result;
}


void ComputeCN(const short *counts, const double *pedestal, const bool *CN_mask, double *CN) {
  double sumVA[N_VA];
  int countVA[N_VA];
  for(int iVA=0; iVA<N_VA; ++iVA) {
    sumVA[iVA]=0.;
    countVA[iVA]=0;
  }
    
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    if(CN_mask[iChan]==false) continue;
    int iVA=ChanToVA(iChan);
    sumVA[iVA]+=(static_cast<double>(counts[iChan])-pedestal[iChan]);
    ++countVA[iVA];
  }

  for(int iVA=0; iVA<N_VA; ++iVA) {
    CN[iVA]=(sumVA[iVA]/countVA[iVA]);
    if(countVA[iVA]==0) std::cout << "Warning! ComputeCN dividing by zero!" << std::endl;
  }

  return;
}


