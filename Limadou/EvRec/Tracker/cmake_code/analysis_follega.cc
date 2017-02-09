#include "LTrackerTools.hh"
#include "analysis_follega.hh"
#include "LPlotTools.hh"
#include "LTrackerFittingTools.hh"
#include "LTrackerCalibrationManager.hh"

#include "TStyle.h"

struct event {
  int entry;
  std::vector<LTrackerCluster> cls;
};

				//==================================================================================================
				//========================================== MAIN =================================================
				//==================================================================================================

int main(int argc, char *run[]){
//calibration procedure 
  if((argc!=4)&&(argc!=5)) {
    std::cerr << "Error! Usage:    ./analysis_follega <calRunFile> <calOutFile> <RunFile>" << std::endl;
    std::cerr << "Aborted." << std::endl;
    return -999;
  }

  int pom = 1;
  if(argc==5) pom = 0;
  
  LTrackerCalibrationManager::GetInstance().LoadRun(run[1]);
  LTrackerCalibration *cal;


  double * PED_CALIB;
  double * SIGMA_CALIB;
  double * NOGAUSS_CALIB;
  bool * CNMASK;

  if (pom==0){
    cal = LTrackerCalibrationManager::GetInstance().Calibrate(NCALIBEVENTS,0);
    cal->Write(run[2]);
  }
  if (pom==1){ 
    cal = LTrackerCalibrationManager::GetInstance().Calibrate();
    cal->Write(run[2]);
    PED_CALIB = cal->GetPedestal(0);
    SIGMA_CALIB = cal->GetSigma(0);
    NOGAUSS_CALIB =  cal->GetNGIndex(0);
    CNMASK = cal->GetCNMask(0);
  }

  std::string NAME = "";
  
  TFile * data = new TFile(run[3],"READ");
  TTree * datatree = (TTree*)data->Get("T");
  TH1F * h_eta_distribution[2][6];
  TH2F * h_counts_distribution[2][6];
  TH2F * h_sn_distribution[2][6];
  TProfile * h_counts_profile[2][6];
  TH1F * h_sumcounts_distribution[2][6];
  TH1F * h_sumcounts_distribution_corr[2][6];
  for (int side = 0; side < 2; side ++){
    std::string NAMESIDE = "p";
    if (side==1) NAMESIDE = "n";
    for (int w=0; w<6; w++){
      NAME = "eta_distribution_"+to_string(w)+"_side_"+NAMESIDE;
      h_eta_distribution[side][w] = new TH1F(NAME.c_str(),NAME.c_str(),100,-0.1,1.1);
      NAME ="sn_distribution_"+to_string(w)+"_side_"+NAMESIDE;
      h_sn_distribution[side][w] = new TH2F(NAME.c_str(),NAME.c_str(),200,-0.1,1.1,300,-20,50);
      NAME ="counts_distribution_"+to_string(w)+"_side_"+NAMESIDE;
      h_counts_distribution[side][w] = new TH2F(NAME.c_str(),NAME.c_str(),200,-0.1,1.1,300,-20,600);
      NAME ="sumcounts_distribution_"+to_string(w)+"_side_"+NAMESIDE;
      h_sumcounts_distribution[side][w] = new TH1F(NAME.c_str(),NAME.c_str(),300,0,600);
      NAME ="sumcounts_distribution_corr"+to_string(w)+"_side_"+NAMESIDE;
      h_sumcounts_distribution_corr[side][w] = new TH1F(NAME.c_str(),NAME.c_str(),300,0,600);
      NAME ="counts_profdistribution_"+to_string(w)+"_side_"+NAMESIDE;
      h_counts_profile[side][w] = new TProfile(NAME.c_str(),NAME.c_str(),10,0.,1.);
    }
  }

  short strip[n_chann];
  UInt_t event_index;
  datatree->SetBranchAddress("strip[4608]",&strip);// choose the branch 
  datatree->SetBranchAddress("event_index",&event_index);
  datatree->SetBranchStatus("*",kFALSE);// close the other branch
  datatree->SetBranchStatus("strip[4608]",kTRUE);

  short * ADC_COUNTS[NCALIBEVENTS];
  double * CNOISE[NCALIBEVENTS];
  double * CLEARCOUNTS[NCALIBEVENTS];
  for (int iev = 0; iev < NCALIBEVENTS; ++iev) ADC_COUNTS[iev] = new short[NCHAN];
  for (int iev = 0; iev < NCALIBEVENTS; ++iev) CNOISE[iev] = new double[NCHAN];
  for (int iev = 0; iev < NCALIBEVENTS; ++iev) CLEARCOUNTS[iev] = new double[NCHAN];

  int NSLOTS;
  NSLOTS= int(datatree->GetEntries()/NCALIBEVENTS);
  std::vector< event > clev;

  for( int slot = 0; slot < NSLOTS; ++slot){
    if(pom==0){
      PED_CALIB = cal->GetPedestal(slot);
      SIGMA_CALIB = cal->GetSigma(slot);
      NOGAUSS_CALIB =  cal->GetNGIndex(slot);
      CNMASK = cal->GetCNMask(slot);
    }

    for (int iev = 0; iev < NCALIBEVENTS; ++iev){
      datatree->GetEntry(slot*NCALIBEVENTS+iev);
      for (int ichan = 0; ichan < NCHAN; ++ichan){
        ADC_COUNTS[iev][ichan] = strip[ichan];
      }
      ComputeCN(ADC_COUNTS[iev], PED_CALIB, CNMASK, &CNOISE[iev][0]);
      for(int ichan=0; ichan<NCHAN; ++ichan){
        CLEARCOUNTS[iev][ichan] = ADC_COUNTS[iev][ichan]-PED_CALIB[ichan]-CNOISE[iev][ChanToVA(ichan)];
      }
    }

//********** Improved cluster finding algorithm ************ 

    for(int iev=0;iev<NCALIBEVENTS;++iev){
      std::vector<LTrackerCluster> *clusters = GetClusters(CLEARCOUNTS[iev],SIGMA_CALIB);
      //std::cout<<"size di clusters "<< clusters->size() <<std::endl;
      event myevent;
      myevent.entry=iev;
      //cit is the cluster
      for(auto cit : *clusters){
        myevent.cls.push_back(cit);
        //eta and sum of the counts for clusters
        int side = ChanToSide(cit.seed);
        int ladder = ChanToLadder(cit.seed);
        double eta = cit.GetEta();
        double sumclscount = 0;
        double sumclssn = 0;
        for(auto clssn : cit.sn){
          sumclssn += clssn;
        }
        sumclscount = cit.GetCounts(0.);

        for (auto chansn : cit.sn) {
          sumclssn += chansn;
        }
        // check the borders
        if(sumclssn<-600 || sumclscount<-600) continue; 

        h_eta_distribution[side][ladder]->Fill(eta);
        h_sumcounts_distribution[side][ladder]->Fill(sumclscount);
        h_counts_distribution[side][ladder]->Fill(eta,sumclscount);
        h_sn_distribution[side][ladder]->Fill(eta,sumclssn);

        if(((eta>-0.1)&&(eta<0.1))||((eta>0.4)&&(eta<0.6))||((eta>0.9)&&(eta<1.1))){
          h_counts_profile[side][ladder]->Fill(eta,sumclscount);
        }

      }
      clev.push_back(myevent);
    }
  }
  //*************** End cluster finding algorithm and of the loop over all the events ****************
  //********** Fitting and Correction, we need to loop over the clusters again...what a waste! ************
  double P0[2][6];
  double P1[2][6];
  double P2[2][6];
  
  Fitting_func(h_counts_profile[0], &P0[0][0], &P1[0][0], &P2[0][0]);
  Fitting_func(h_counts_profile[1], &P0[1][0], &P1[1][0], &P2[1][0]);
  
  for(auto ientry : clev){
    for(auto cit : ientry.cls){
      int side = ChanToSide(cit.seed);
      int ladder = ChanToLadder(cit.seed);
      double sumclscount = 0;
      sumclscount = cit.GetCounts(0.);
      if(sumclscount<-600) continue; 
      // check the borders

      if(P2[side][ladder]>0){
        h_sumcounts_distribution_corr[side][ladder]->Fill(sumclscount*Correction_1(cit.GetEta(),P0[side][ladder],P1[side][ladder],P2[side][ladder])); 
      }
      else{
        h_sumcounts_distribution_corr[side][ladder]->Fill(sumclscount);
      } 
    }    
  } 

  double gaus_Mean[2][6];
  double gaus_Sigma[2][6];
  double chi2red[2][6];

  if(pom == 1){
    Langau(h_sumcounts_distribution_corr[0], &gaus_Mean[0][0], &gaus_Sigma[0][0], &chi2red[0][0],10.,600.);
    Langau(h_sumcounts_distribution_corr[1], &gaus_Mean[1][0], &gaus_Sigma[1][0], &chi2red[1][0],10.,600.);
  }
  else{
    Langau(h_sumcounts_distribution[0], &gaus_Mean[0][0], &gaus_Sigma[0][0], &chi2red[0][0],0.,110.);
    Langau(h_sumcounts_distribution[1], &gaus_Mean[1][0], &gaus_Sigma[1][0], &chi2red[1][0],0.,110.);
  }
  //*************** End Fitting and Correction ****************
  
  std::ofstream langau_paremters;
  langau_paremters.open("langau_param.txt", fstream::app);

  for(int i_par=0; i_par<6; i_par++){
    langau_paremters<<gaus_Mean[0][i_par]<<" "<<gaus_Sigma[0][i_par]<<" "<<chi2red[0][i_par]<<" "<<gaus_Mean[1][i_par]<<" "<<gaus_Sigma[1][i_par]<<" "<<chi2red[0][i_par]<<std::endl;
  }

  TCanvas*c1 = new TCanvas();
  gStyle->SetPalette(1);
  string NAMEPLOTbegin = run[2];
  NAMEPLOTbegin += ".pdf[";
  string NAMEPLOT = run[2];
  NAMEPLOT += ".pdf";
  string NAMEPLOTend = run[2];
  NAMEPLOTend += ".pdf]";

  c1->Print(NAMEPLOTbegin.c_str());
  for (int side=0; side<2; ++side){
    Plot_6histo_1d(h_eta_distribution[side],"eta","entries","")->Print(NAMEPLOT.c_str());
    Plot_6histo_2d(h_counts_distribution[side],"eta","counts")->Print(NAMEPLOT.c_str());
    Plot_6histo_2d(h_sn_distribution[side],"eta","S/N")->Print(NAMEPLOT.c_str());
    Plot_6histo_profile(h_counts_profile[side],"eta","counts","")->Print(NAMEPLOT.c_str());
    Plot_6histo_1d(h_sumcounts_distribution[side],"sumclscount","entries","")->Print(NAMEPLOT.c_str());
    Plot_6histo_1d(h_sumcounts_distribution_corr[side],"sumclscount","entries","")->Print(NAMEPLOT.c_str());
  }
  c1->Print(NAMEPLOTend.c_str());
}






		




