#ifndef LIMADOU_ANALYSIS_CONST
#define LIMADOU_ANALYSIS_CONST
const int NCALIBEVENTS = 1000;
const double SIGMA_CUT_CN=5.; //at first 3
const double SIGMA_CUT_SIGMA=3.;
const double SIGMA_SEL=1.5;
const double GE_THRESHOLD=6.5; //good events threshold
const int BIN_SIGMA=500;
const double MAX_SIGMA=100.;
const double SIGMA_STEP=MAX_SIGMA/(double)BIN_SIGMA;

const double MAX_SIGMA2=300.;
const int CLUSTERSIZE=5;
const int NEV_WINDOW = CLUSTERSIZE+2;

const double MIN_NCHANSPERBIN_4CN=5;
const double MIN_SIGMA_CHANNEL_OFF=6;
const double MAX_SIGMA_CHANNEL_OFF=30.;
const double MAX_SIGMA_NOISE_LEVEL=4;

const double NSIGMABOOL=3.;
const double GAUSSTHRESHOLD=0.0027;
const double GAUSINDEX_CUT=3.;
const double SIGMAHOT=30.;



//eta manage

const int ETASTEP=200;
const int ETAMAX=1;
const int ETAMIN=0;
const int ETARANGE=ETAMAX-ETAMIN;
#endif
