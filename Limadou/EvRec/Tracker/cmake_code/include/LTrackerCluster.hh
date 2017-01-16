#ifndef __LTRACKERCLUSTER__
#define __LTRACKERCLUSTER__ 1

const int CLUSTERCHANNELS=5;
const int NEV_WINDOW = CLUSTERCHANNELS+2;
const double MAX_SIGMA2=300.;

const double CLFINDTHRESHOLD=3.;
const double CLSNTHRESHOLD=3.5;


class LTrackerCluster {
public:
  LTrackerCluster();
  LTrackerCluster(int inpSeed, double *inpCont, double *inpSigma);
  int seed;
  double count[CLUSTERCHANNELS];
  double sigma[CLUSTERCHANNELS];
  double sign[CLUSTERCHANNELS];
  double eta;

  double GetSides(const double SideThreshold);
  double GetCounts(const double SideThreshold);
  double ChargeCenter(const double SideThreshold);
  int ClusterSize(const double SideThreshold);
  bool CheckIfAtBorder();
};


#endif
