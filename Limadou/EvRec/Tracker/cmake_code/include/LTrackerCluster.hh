#ifndef __LTRACKERCLUSTER__
#define __LTRACKERCLUSTER__ 1

const int CLUSTERCHANNELS=5;
const int NEV_WINDOW = CLUSTERCHANNELS+2;
const double MAX_SIGMA2=300.;


class LTrackerCluster {
public:
  LTrackerCluster();
  int seed;
  double count[CLUSTERCHANNELS];
  double sigma[CLUSTERCHANNELS];
  double sign[CLUSTERCHANNELS];
  double eta;

  double GetSides(const double SideThreshold);
  double GetCounts(const double SideThreshold);
  double ChargeCenter(const double SideThreshold);
  int ClusterSize(const double SideThreshold);
};


#endif
