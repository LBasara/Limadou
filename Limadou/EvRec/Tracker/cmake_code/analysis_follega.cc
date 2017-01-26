#include "LEvRec0.hh"
#include "analysis_const.hh"
#include "detector_const.hh"
#include "LTrackerTools.hh"
#include "LTrackerCluster.hh"
#include "analysis_follega.hh"
#include "LPlotTools.hh"
#include "LTrackerFittingTools.hh"
#include "LTrackerCalibrationManager.hh"
#include "LTrackerMask.hh"
#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

struct event {
  int entry;
  std::vector<LTrackerCluster> cls;
};


				//==================================================================================================
				//========================================== MAIN =================================================
				//==================================================================================================


int main(int argc, char *run[]){
//calibration procedure 
  if(argc!=3) {
    std::cerr << "Error! Usage:    ./analysis_follega <calRunFile> <calOutFile> <RunFile>" << std::endl;
    std::cerr << "Aborted." << std::endl;
    return -999;
  }
  LTrackerCalibrationManager::GetInstance().LoadRun(run[1]);
  LTrackerCalibration *cal =   LTrackerCalibrationManager::GetInstance().Calibrate();
  cal->Write(run[2]);

  const double * PED_CALIB = cal->GetPedestal(0);
  const double * SIGMA_CALIB = cal->GetSigma(0);
  const double * NOGAUSS_CALIB =  cal->GetNGIndex(0);
  const bool * CNMASK = cal->GetCNMask(0);


  TFile * data = new TFile(run[3],"READ");
  TTree * datatree = (TTree*)data->Get("T");

  short strip[n_chann];
  UInt_t event_index;
  datatree->SetBranchAddress("strip[4608]",&strip);// choose the branch 
  datatree->SetBranchAddress("event_index",&event_index);
  datatree->SetBranchStatus("*",kFALSE);// close the other branch
  datatree->SetBranchStatus("strip[4608]",kTRUE);

  short * ADC_COUNTS[NCALIBEVENTS];
  double * CNOISE[NCALIBEVENTS];
  for (int iev = 0; iev < NCALIBEVENTS; ++iev) ADC_COUNTS[iev] = new short[NCHAN];
  for (int iev = 0; iev < NCALIBEVENTS; ++iev) CNOISE[iev] = new double[NCHAN];


  int NSLOTS;
  NSLOTS= int(datatree->GetEntries()/NCALIBEVENTS);

  for( int slot = 0; slot < NSLOTS; ++slot){

    for (int iev = 0; iev < NCALIBEVENTS; ++iev){
      datatree->GetEntry(slot*NCALIBEVENTS+iev);

      for (int ichan = 0; ichan < NCHAN; ++ichan){
        ADC_COUNTS[iev][ichan] = strip[ichan];
      }

      ComputeCN(ADC_COUNTS[iev], PED_CALIB, CNMASK, &CNOISE[iev][0]);

    }



  }





}






		




