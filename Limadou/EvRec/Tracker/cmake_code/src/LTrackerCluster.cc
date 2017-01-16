#include "detector_const.hh"
#include "analysis_const.hh"

#include "LTrackerTools.hh"
#include "LTrackerCluster.hh"
#include <iostream>

bool LTrackerCluster::CheckIfAtBorder() {
  bool result = (ChanToLadderChan(seed)==0 || ChanToLadderChan(seed)==LADDER_CHAN-1);
  return result;
}


LTrackerCluster::LTrackerCluster(){
  seed=-999;
  for(int i=0; i<CLUSTERCHANNELS; ++i) {
    count[i]=-999.;
    sigma[i]=-999.;
    sign[i]=-999.;
  }
  eta=-999.;
  
}

LTrackerCluster::LTrackerCluster(int inpSeed, double *inpCont, double *inpSigma){
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
      sign[i]=-999.;
    } else {
      count[i]=inpCont[inpSeed-CLUSTERCHANNELS/2+i];
      sigma[i]=inpSigma[inpSeed-CLUSTERCHANNELS/2+i];;
      sign[i]=count[i]/sigma[i];
    }
  }
  eta=-999.;
}

double LTrackerCluster::GetSides(const double SideThreshold){
    double result=0.;
    if(sign[0]>SideThreshold) result += count[CLUSTERCHANNELS/2-1];
    if(sign[2]>SideThreshold) result += count[CLUSTERCHANNELS/2+1];
    return result;
  }

double LTrackerCluster::GetCounts(const double SideThreshold){
  double result=0.;
  for(int i=0;i<CLUSTERCHANNELS;++i)
    if(sign[i]>SideThreshold) result += count[i];
  
  return result;
}

double LTrackerCluster::ChargeCenter(const double SideThreshold){
  double result=0.;
  double tot_sign=0.;
  for(int i=0;i<CLUSTERCHANNELS;++i){
    if(sign[i]>SideThreshold) {
      result += i*count[i];
      tot_sign +=count[i];
    }
  }
  return result/tot_sign;
}

int LTrackerCluster::ClusterSize(const double SideThreshold){
  int size=0;
  for(int i=0;i<CLUSTERCHANNELS;++i)
    if(sign[i]>SideThreshold) ++size;
  
  return size;
}
