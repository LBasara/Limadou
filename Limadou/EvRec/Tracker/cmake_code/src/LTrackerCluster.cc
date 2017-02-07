#include "detector_const.hh"
#include "analysis_const.hh"

#include "LTrackerTools.hh"
#include "LTrackerCluster.hh"
#include <iostream>
#include <cmath>

double LTrackerCluster::ComputeEta() {
  int seedIndex=CLUSTERCHANNELS/2;
  int max2Index=seedIndex;
  if(IsAtStartBorder()) ++max2Index;
  else if(IsAtStopBorder()) --max2Index;
  else max2Index += (sn[max2Index-1]>sn[max2Index+1] ? -1 : 1);
  
  double denominator=count[seedIndex]+count[max2Index];
  double numerator=(max2Index>seedIndex ? count[max2Index] : count[seedIndex]);
 
  eta = numerator/denominator;
  etaCounts=denominator;
  
  return eta;
}

double LTrackerCluster::ComputeEta3() {
  int seedIndex=CLUSTERCHANNELS/2;
  
   double numerator=/*0*count[seedIndex-1]+*/count[seedIndex]+2*count[seedIndex+1];
  double denominator=count[seedIndex-1]+count[seedIndex]+count[seedIndex+1];

  eta = numerator/denominator;
  etaCounts=denominator;
  
  return eta;
}


bool LTrackerCluster::IsAtStartBorder() {
  bool result = (ChanToLadderChan(seed)==0);
  return result;
}

bool LTrackerCluster::IsAtStopBorder() {
  bool result = (ChanToLadderChan(seed)==LADDER_CHAN-1);
  return result;
}

bool LTrackerCluster::IsAtBorder() {
  bool result = (IsAtStartBorder() || IsAtStopBorder());
  return result;
}


LTrackerCluster::LTrackerCluster(){
  seed=-999;
  for(int i=0; i<CLUSTERCHANNELS; ++i) {
    count[i]=-999.;
    sigma[i]=-999.;
    sn[i]=-999.;
  }
  eta=-999.;
  etaCounts=-999.;
  
}

LTrackerCluster::LTrackerCluster(const int inpSeed, const double *inpCont, const double *inpSigma){
  if(inpSeed<0 || inpSeed>NCHAN) {
    LTrackerCluster();
    std::cerr << "LTrackerCluster::LTrackerCluster(int inpSeed, double *inpCont, double *inpSigma)   error:"
	     << " inpSeed=" << inpSeed << ". Default constructor used."
	     << std::endl;
    return;
  }

  seed=inpSeed;
  for(int i=0; i<CLUSTERCHANNELS; ++i) {
    int index=ChanToLadderChan(inpSeed)-CLUSTERCHANNELS/2+i;
    if(index<0 || index>LADDER_CHAN-1) {
      count[i]=-999.;
      sigma[i]=-999.;
      sn[i]=-999.;
    } else {
      count[i]=inpCont[inpSeed-CLUSTERCHANNELS/2+i];
      sigma[i]=inpSigma[inpSeed-CLUSTERCHANNELS/2+i];;
      sn[i]=count[i]/sigma[i];
    }
  }
  ComputeEta();
  //ComputeEta3();
}

double LTrackerCluster::GetSides(const double SideThreshold){
    double result=0.;
    if(sn[0]>SideThreshold) result += count[CLUSTERCHANNELS/2-1];
    if(sn[2]>SideThreshold) result += count[CLUSTERCHANNELS/2+1];
    return result;
  }

double LTrackerCluster::GetCounts(const double SideThreshold){
  double result=0.;
  for(int i=0;i<CLUSTERCHANNELS;++i)
    if(sn[i]>SideThreshold) result += count[i];
  
  return result;
}

double LTrackerCluster::ChargeCenter(const double SideThreshold){
  double result=0.;
  double tot_sn=0.;
  for(int i=0;i<CLUSTERCHANNELS;++i){
    if(sn[i]>SideThreshold) {
      result += i*count[i];
      tot_sn +=count[i];
    }
  }
  return result/tot_sn;
}

int LTrackerCluster::ClusterSize(const double SideThreshold){
  int size=0;
  for(int i=0;i<CLUSTERCHANNELS;++i)
    if(sn[i]>SideThreshold) ++size;
  
  return size;
}

void LTrackerCluster::Dump() {
  std::cout << "***************" << std::endl
            <<"\t seed " << seed << " (side " << ChanToSide(seed) << ", ladder " << ChanToLadder(seed) << ")" << std::endl
            <<"\t eta " << eta << std::endl << "\t";
for(int i = 0; i < CLUSTERCHANNELS; ++i)  std::cout << count[i] << " ";
std::cout << std::endl << "\t";
for(int i = 0; i < CLUSTERCHANNELS; ++i)  std::cout << sn[i] << " ";
 
  std::cout << std::endl << "***************" << std::endl;
  return;
}
