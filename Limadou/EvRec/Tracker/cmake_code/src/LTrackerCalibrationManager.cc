#include "LTrackerCalibrationManager.hh"
#include "LEvRec0File.hh"
#include "LEvRec0.hh"

#include<math.h>

LTrackerCalibrationManager::LTrackerCalibrationManager() {
  calRunFileName=0;
}


LTrackerCalibrationManager& LTrackerCalibrationManager::GetInstance() {
  static LTrackerCalibrationManager instance; // Guaranteed to be destroyed.
                                              // Instantiated on first use.
  return instance;
}


int LTrackerCalibrationManager::LoadRun(char *fileInp) {
  calRunFileName = fileInp;
  return 0;
}


LTrackerCalibration* LTrackerCalibrationManager::Calibrate(int nEvents, int skipEvents) {
  // Open the run file and get the number of entries
  LEvRec0File calRunFile(calRunFileName);
  const int nEntries=calRunFile.GetEntries();
  if(nEntries<skipEvents+nEvents) return 0;
  int nSlots = CalculateCalibrationSlots(nEvents, skipEvents, nEntries, pivot);

  LTrackerCalibration *result = CreateTrackerCalibration();
  for(int is=0; is<nSlots+1; ++is) result->Add(CalibrateSlot(pivot[is],pivot[is+1]));
  return result;
}

void LTrackerCalibrationManager::CalculateCalibrationSlots(int nEvents, int skipEvents, int nEntries, int *pivot) {
  if(nEvents==-1) nEvents=nEntries;
  if(skipEvents==-1) skipEvents=0;
  const int nSlots = (nEntries-skipEvents)/nEvents;
  pivot = new int[nSlots+1];
  for(int is=0; is<nSlots+1; ++is) pivot[is]=skipEvents+is*nEvents;
  return nSlots;
}

LTrackerCalibrationSlot* LTrackerCalibrationManager::CalibrateSlot(int StartEntry, int StopEntry) {
  // Open the run file
  LEvRec0File calRunFile(calRunFileName);
  const int nEntries=calRunFile.GetEntries();
  //CalculateCalibrationSlots();

  LEvRec0 cev;
  calRunFile.SetTheEventPointer(cev);

  // Set up arrays - 0 level
  double mean0[NCHAN];
  double sumsq0[NCHAN];
  double sigma0[NCHAN];
  int counter0[NCHAN];
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    mean0[iChan]=0.;
    sumsq0[iChan]=0.;
    count0[iChan]=0;
  }
  
  // Average counts and squares
  for(int iEntry=StartEntry; iEntry<StopEntry; ++iEntry) {
    calRunFile.GetEntry(iEntry);
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

  // Set up arrays - 1 level
  double mean1[NCHAN];
  double sumsq1[NCHAN];
  double sigma1[NCHAN];
  int counter1[NCHAN];
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    mean1[iChan]=0.;
    sumsq1[iChan]=0.;
    counter1[iChan]=0;
  }
  
  // Average counts and squares
  for(int iEntry=StartEntry; iEntry<StopEntry; ++iEntry) {
    calRunFile.GetEntry(iEntry);
    for(int iChan=0; iChan<NCHAN; ++iChan) {
      double x = (static_cast<double>(cev.strip[iChan])-mean0[iChan]);
      if(abs(x)>CHANCLEANINGTHRESHOLD*sigma0[iChan]) continue;
      mean1[iChan]+=x;
      sumsq1[iChan]+=(x*x);
      ++counter1[iChan];
    }
  }
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    mean1[iChan]/=counter1[iChan];
    sigma1[iChan]=sqrt(sumsq1[iChan]/counter1[iChan]-mean1[iChan]*mean1[iChan]);
  }
  
  // Compute CN mask 
  bool CN_mask[NCHAN];
  ComputeCNMask(sigma1,&CN_mask[0]);

  // Set up arrays - 2 level
  double mean2[NCHAN];
  double sumsq2[NCHAN];
  double sigma2[NCHAN];
  int counter2[NCHAN];
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    mean2[iChan]=0.;
    sumsq2[iChan]=0.;
    counter2[iChan]=0;
  }
  
  // Average counts and squares
  for(int iEntry=StartEntry; iEntry<StopEntry; ++iEntry) {
    calRunFile.GetEntry(iEntry);
    double CN[N_VA];
    ComputeCN(cev.strip,mean0,CN_mask,&CN[0]);
    for(int iChan=0; iChan<NCHAN; ++iChan) {
      double x = (static_cast<double>(cev.strip[iChan])-mean0[iChan]-CN[iChan]);
      mean2[iChan]+=x;
      sumsq2[iChan]+=(x*x);
      ++counter2[iChan];
    }
  }
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    mean2[iChan]/=counter2[iChan];
    sigma2[iChan]=sqrt(sumsq2[iChan]/counter2[iChan]-mean2[iChan]*mean2[iChan]);
  }
  
  // Gaussianity
  double ngindex[NCHAN];
  int ngcounter[NCHAN];
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    ngindex[iChan]=0.;
    ngcounter[iChan]=0;
  }
  
  for(int iEntry=StartEntry; iEntry<StopEntry; ++iEntry) {
    calRunFile.GetEntry(iEntry);
    double CN[N_VA];
    ComputeCN(cev.strip,mean0,CN_mask,&CN[0]);
    for(int iChan=0; iChan<NCHAN; ++iChan) {
      double x = (static_cast<double>(cev.strip[iChan])-mean0[iChan]-CN[iChan]);
      if(abs(x)>GAUSSIANITYSIGMATHRESHOLD*sigma2[iChan]) ++ngindex[iChan];
      ++ngcounter[iChan];
    }
  }
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    double outliers_expected=GAUSSIANITYEVRACTHRESHOLD*ngcounter[iChan];
    double delta=ngindex[iChan]-outliers_expected;
    ngindex[iChan]=(delta/sqrt(outliers_expected));
  }

  calRunFile.GetEntry(StartEntry);
  int StartEvent=cev.event_index;
  calRunFile.GetEntry(StopEntry-1);
  int StopEvent=cev.event_index;
  LTrackerCalibrationSlot *result = new LTrackerCalibrationSlot(StartEvent, StopEvent, mean0, sigma2, ngindex);
  return result;
}




void LTrackerCalibrationManager::ComputeCNMask(double *sigma1, double &CN_mask) {
  // Create "histograms" of sigmas per VA and set them to zero
  int hSigma1[N_VA][NSIGMA1BIN+1]; // overflow
  for(int iVA=0; iVA<N_VA; ++iVA)
    for(int iBin=0; iBin<NSIGMA1BIN; ++iBin)
      hSigma1[iVA][iBin] = 0;
  
  double xBin[NSIGMA1BIN+1];
  for(int iBin=0; iBin<NSIGMA1BIN; ++iBin)
    xBin[iBin]=MINSIGMA1+(MAXSIGMA1-MINSIGMA1)*iBin;
  
  // Fill histograms
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    int iVA=ChanToVA(iChan);
    int bin=NSIGMA1BIN;
    for(int iBin=NSIGMA1BIN; iBin>-1; --iBin) {
      if(sigma2[iChan]>xBin[iBin]) {
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
	hMaxima[iVA]=xBin[iBin]+0.5*(MAXSIGMA1-MINSIGMA1);
	tmpMax=hSigma1[iVA][iBin];
      }
    }
  }
  
  // Compute CN mask 
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    int iVA=ChanToVA(iChan);
    if(abs(sigma2[iChan]-hMaxima[iVA])>HALFSIGMA1WIDTH) CN_mask[iChan]=false;
    else CN_mask[iChan]=true;
  }  
  
  return;
}



void LTrackerCalibrationManager::ComputeCN(double *counts, double *pedestal, bool *CN_mask, double &CN) {
  double sumVA[N_VA];
  int countVA[N_VA];
  for(int iVA=0; iVA<N_VA; ++iVA) {
    sumVA[iVA]=0.;
    countVA[iVA]=0;
  }
    
  for(int iChan=0; iChan<NCHAN; ++iChan) {
    int iVA=ChanToVA(iChan);
    if(CN_mask[iChan]==false) continue;
    sumVA[iVA]+=(counts[iChan]-pedestal[iChan]);
    ++countVA[iVA];
  }

  for(int iVA=0; iVA<N_VA; ++iVA) {
    CN[iVA]=(sumVA[iVA]/countVA[iVA]);
  }
  return;
}
