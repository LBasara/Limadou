#include "LTrackerCalibrationManager.hh"
#include "LEvRec0File.hh"
#include "LEvRec0.hh"
#include "LTrackerTools.hh"


#include <iostream>
#include <math.h>

LTrackerCalibrationManager::LTrackerCalibrationManager() {
  calRunFile=0;
  InitialTargetRun=-1;
  FinalTargetRun=-1;
  verboseFLAG=true;
}


LTrackerCalibrationManager& LTrackerCalibrationManager::GetInstance() {
  static LTrackerCalibrationManager instance; // Guaranteed to be destroyed.
                                              // Instantiated on first use.
  return instance;
}


int LTrackerCalibrationManager::LoadRun(const char *fileInp) {
  // Reset previous loaded runs
  if(calRunFile) {
    if(calRunFile->IsOpen()) calRunFile->Close();
    calRunFile = 0;
    // Sure we want also to reset the target runs? Today we reply yes... check
    InitialTargetRun=-1;
    FinalTargetRun=-1;
  }
  
  calRunFile = new LEvRec0File(fileInp);
  if(!calRunFile || !(calRunFile->IsOpen())) {
    std::cerr << "Error! Attempt to load a calibration run, but the file cannot be opened."
	      << std::endl;
    calRunFile = 0;
    return -999;
  }
  return 0;
}

void LTrackerCalibrationManager::SetTargetRuns(const int InitialRun, const int FinalRun) {
  InitialTargetRun=InitialRun;
  FinalTargetRun=FinalRun;
  return;
}

LTrackerCalibration* LTrackerCalibrationManager::Calibrate(const int nEvents, const int skipEvents) {
  
  if(calRunFile==0 || !(calRunFile->IsOpen())) {
    std::cerr << "Error! Attempt to call the \"Calibrate\" method, but no calibration run loaded."
	      << std::endl;
    return 0;
  }
  
  int nEntries=calRunFile->GetEntries();
  if(nEntries<skipEvents+nEvents) return 0;
  int *pivot=0;
  int nSlots = CalculateCalibrationSlots(nEvents, skipEvents, nEntries, pivot);

  LTrackerCalibration *result = CreateTrackerCalibration();
  for(int is=0; is<nSlots; ++is) result->Add(CalibrateSlot(pivot[is],pivot[is+1]));

  return result;
}

int LTrackerCalibrationManager::CalculateCalibrationSlots(const int nEvents, const int skipEvents, const int nEntries, int* &pivot) {
  int nEv=nEvents;
  int skipEv=skipEvents;
  if(nEvents==-1) nEv=nEntries;
  if(skipEvents==-1) skipEv=0;
  const int nSlots = (nEntries-skipEv)/nEv;
  pivot = new int[nSlots+1];
  for(int is=0; is<nSlots+1; ++is) pivot[is]=skipEv+is*nEv;

  return nSlots;
}

LTrackerCalibrationSlot* LTrackerCalibrationManager::CalibrateSlot(const int StartEntry, const int StopEntry) {

  if(calRunFile==0 || !(calRunFile->IsOpen())) {
    std::cerr << "Error! Attempt to call the \"CalibrateSlot\" method, but no calibration run loaded."
	      << std::endl;
    return 0;
  }

  // RawMeanSigma
  double mean0[NCHAN];
  double sigma0[NCHAN];
  RawMeanSigma(StartEntry, StopEntry, mean0, sigma0);
  
  // First cleaning
  double mean1[NCHAN];
  double sigma1[NCHAN];
  CleanedMeanSigma(StartEntry, StopEntry, mean0, sigma0, mean1, sigma1);

  // Compute CN mask 
  bool CN_mask[NCHAN];
  ComputeCNMask(sigma1,CN_mask);

  // CNCorrectedSigma
  double mean2[NCHAN];
  double sigma2[NCHAN];
  CNCorrectedSigma(StartEntry, StopEntry, mean1, sigma1, CN_mask, mean2, sigma2);
  
  // Gaussianity
  double ngindex[NCHAN];
  GaussianityIndex(StartEntry, StopEntry, mean2, sigma2, CN_mask, ngindex);

  // Start and stop events
  LEvRec0 cev;
  calRunFile->SetTheEventPointer(cev);  
  calRunFile->GetEntry(StartEntry);
  int StartEvent=static_cast<int>(cev.event_index);
  calRunFile->GetEntry(StopEntry-1);
  int StopEvent=static_cast<int>(cev.event_index);

  // Result
  LTrackerCalibrationSlot *result = new LTrackerCalibrationSlot(StartEvent, StopEvent, mean2, sigma2, ngindex, CN_mask);
  return result;
}




void LTrackerCalibrationManager::ComputeCNMask(const double *sigma1, bool *CN_mask) {
  // Create "histograms" of sigmas per VA and set them to zero
  int hSigma1[N_VA][NSIGMA1BIN+1]; // overflow
  for(int iVA=0; iVA<N_VA; ++iVA)
    for(int iBin=0; iBin<NSIGMA1BIN; ++iBin)
      hSigma1[iVA][iBin] = 0;
  
  double xBin[NSIGMA1BIN+1];
  for(int iBin=0; iBin<NSIGMA1BIN+1; ++iBin) xBin[iBin]=MINSIGMA1+(MAXSIGMA1-MINSIGMA1)*iBin/static_cast<double>(NSIGMA1BIN);

  // Fill histograms
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    int iVA=ChanToVA(iChan);
    for(int iBin=NSIGMA1BIN; iBin>-1; --iBin) {
      if(sigma1[iChan]>xBin[iBin]) {
	++hSigma1[iVA][iBin];
	break;
      }
    }
  }

  // Find histos' maxima
  double hMaxima[N_VA];
  for(int iVA=0; iVA<N_VA; ++iVA) {
    int tmpMax=-99999;
    hMaxima[iVA]=-99999.;
    for(int iBin=0; iBin<NSIGMA1BIN; ++iBin) {
      if(hSigma1[iVA][iBin]>tmpMax) {
	hMaxima[iVA]=xBin[iBin]+0.5*(MAXSIGMA1-MINSIGMA1)/NSIGMA1BIN;
	tmpMax=hSigma1[iVA][iBin];
      }
    }
  }
  
  // Compute CN mask 
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    int iVA=ChanToVA(iChan);
    if(std::fabs(sigma1[iChan]-hMaxima[iVA])>HALFSIGMA1WIDTH) CN_mask[iChan]=false;
    else CN_mask[iChan]=true;
  }  
  
  if(verboseFLAG) {
    std::cout << "CNmask computed" << std::endl;
  }
  return;
}

 
void LTrackerCalibrationManager::RawMeanSigma(const int StartEntry, const int StopEntry, double *mean0, double *sigma0) {
  
  LEvRec0 cev;
  calRunFile->SetTheEventPointer(cev);
  
  // Set up arrays - 0 level
  double sumsq0[NCHAN];
  int counter0[NCHAN];
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    mean0[iChan]=0.;
    sumsq0[iChan]=0.;
    counter0[iChan]=0;
  }
  
  // Average counts and squares
  for(int iEntry=StartEntry; iEntry<StopEntry; ++iEntry) {
    calRunFile->GetEntry(iEntry);
    for(int iChan=0; iChan<NCHAN; ++iChan) {
      auto x = static_cast<double>(cev.strip[iChan]);
      mean0[iChan]+=x;
      sumsq0[iChan]+=(x*x);
      ++counter0[iChan];
    }
  }
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    mean0[iChan]/=counter0[iChan];
    sigma0[iChan]=sqrt(sumsq0[iChan]/counter0[iChan]-mean0[iChan]*mean0[iChan]);
  }
  
  if(verboseFLAG) std::cout << "RawMeanSigma computed" << std::endl;
  return;
}



void LTrackerCalibrationManager::CleanedMeanSigma(const int StartEntry, const int StopEntry, const double *mean0, const double *sigma0, double *mean1, double *sigma1) {

  LEvRec0 cev;
  calRunFile->SetTheEventPointer(cev);
  
  // Set up arrays - 1 level
  double sumsq1[NCHAN];
  int counter1[NCHAN];
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    mean1[iChan]=0.;
    sumsq1[iChan]=0.;
    counter1[iChan]=0;
  }
  
  // Average counts and squares
  for(int iEntry=StartEntry; iEntry<StopEntry; ++iEntry) {
    calRunFile->GetEntry(iEntry);
    for(int iChan=0; iChan<NCHAN; ++iChan) {
      double x = static_cast<double>(cev.strip[iChan]);
      double diff = (x-mean0[iChan]);
      if(std::fabs(diff)>CHANCLEANINGTHRESHOLD*sigma0[iChan]) continue;
      mean1[iChan]+=x;
      sumsq1[iChan]+=(x*x);
      ++counter1[iChan];
    }
  }
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    mean1[iChan]/=counter1[iChan];
    sigma1[iChan]=sqrt(sumsq1[iChan]/counter1[iChan]-mean1[iChan]*mean1[iChan]);
  }
  
  if(verboseFLAG) std::cout << "CleanedMeanSigma computed" << std::endl;
  return;
}



void LTrackerCalibrationManager::CNCorrectedSigma(const int StartEntry, const int StopEntry, const double *mean1, const double *sigma1, const bool *CN_mask, double *mean2, double *sigma2) {
  
  LEvRec0 cev;
  calRunFile->SetTheEventPointer(cev);
  
  // Set up arrays - CN corrected sigma
  double sumsq2[NCHAN];
  int counter2[NCHAN];
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    mean2[iChan]=0.;
    sumsq2[iChan]=0.;
    counter2[iChan]=0;
  }
  
  // Average counts and squares
  for(int iEntry=StartEntry; iEntry<StopEntry; ++iEntry) {
    calRunFile->GetEntry(iEntry);
    double CN[N_VA];
    ComputeCN(cev.strip,mean1,CN_mask,CN);
    for(int iChan=0; iChan<NCHAN; ++iChan) {
      double x = static_cast<double>(cev.strip[iChan]);
      double diff = (x-mean1[iChan]);
      if(std::fabs(diff)>CHANCLEANINGTHRESHOLD*sigma1[iChan]) continue;
      double y = (x-CN[ChanToVA(iChan)]);  // CN corrected!
      mean2[iChan]+=y;
      sumsq2[iChan]+=(y*y);
      ++counter2[iChan];
    }
  }
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    mean2[iChan]/=counter2[iChan];
    sigma2[iChan]=sqrt(sumsq2[iChan]/counter2[iChan]-mean2[iChan]*mean2[iChan]);
  }
  
  if(verboseFLAG) std::cout << "CNCorrectedSigma computed" << std::endl;
  return;
}


void LTrackerCalibrationManager::GaussianityIndex(const int StartEntry, const int StopEntry, const double *mean2, const double *sigma2, const bool *CN_mask, double *ngindex) {
  
  LEvRec0 cev;
  calRunFile->SetTheEventPointer(cev);
  
  // Gaussianity index
  int ngcounter[NCHAN];
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    ngindex[iChan]=0.;
    ngcounter[iChan]=0;
  }
  
  for(int iEntry=StartEntry; iEntry<StopEntry; ++iEntry) {
    calRunFile->GetEntry(iEntry);
    double CN[N_VA];
    ComputeCN(cev.strip,mean2,CN_mask,CN);
    for(int iChan=0; iChan<NCHAN; ++iChan) {
      double x = (static_cast<double>(cev.strip[iChan])-mean2[iChan]-CN[ChanToVA(iChan)]);
      if(std::fabs(x)>GAUSSIANITYSIGMATHRESHOLD*sigma2[iChan]) ++ngindex[iChan];
      ++ngcounter[iChan];
    }
  }
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    double outliers_expected=GAUSSIANITYEVRACTHRESHOLD*ngcounter[iChan];
    double delta=ngindex[iChan]-outliers_expected;
    double denominator=sqrt(outliers_expected+ngindex[iChan]-2*GAUSSIANITYEVRACTHRESHOLD*ngindex[iChan]);
    //double denominator=ngindex[iChan]+outliers_expected;
    ngindex[iChan]=(delta/denominator);
  }
  
  if(verboseFLAG) std::cout << "GaussianityIndex computed" << std::endl;
  return;
}


LTrackerCalibrationManager::~LTrackerCalibrationManager() {
  // do not care about singleton destructor
}


LTrackerCalibration* LTrackerCalibrationManager::CreateTrackerCalibration() {
  if(calRunFile==0 || !(calRunFile->IsOpen())) {
    std::cerr << "Error! Attempt to create a tracker calibration but no calibration run loaded."
	      << std::endl;
    return 0;
  }
  int RunId = calRunFile->GetRunId();
  
  if(InitialTargetRun==-1 || FinalTargetRun==-1) {
    std::cerr << "Warning! Target run interval of current calibration not defined."
	      << std::endl;
  }
  LTrackerCalibration *result = new LTrackerCalibration(RunId, InitialTargetRun, FinalTargetRun);

  return result;
}
