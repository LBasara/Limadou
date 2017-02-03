#include "detector_const.hh"
#include "analysis_const.hh"

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TROOT.h"
#include "TRint.h"

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

#include "LEvRec0File.hh"
#include "LEvRec0.hh"
#include "LTrackerCluster.hh"
#include "LTrackerTools.hh"
#include "LTrackerCalibrationManager.hh"
#include "LTrackerCalibration.hh"


struct event {
  int entry;
  std::vector<LTrackerCluster> cls;
};

void f_eta(int size,double *x_axis,int* dn_deta,double *result){
  double partial_sum=0.;
  double step_int=0.;
  double step=(double)ETARANGE/(double)ETASTEP;
  for(int i=0;i<size;++i){
    step_int=dn_deta[i]*step;
    partial_sum+=step_int;
    result[i]=partial_sum;
  }
  for(int i=0;i<size;++i){
    if(partial_sum>0.)    result[i]/=partial_sum;
  }
  return;
}
//graphic tools

TCanvas *display_ladders1D(TH1D *full_histo,std::string name,std::string title,double xrange=LADDER_BIN){
  TCanvas *display= new TCanvas();
  
  display->Divide(3,2);
  TH1D *ladder[N_LADDER];
  for(int ld=0;ld<N_LADDER;++ld){
     std::stringstream Stream;
     Stream<<name<<"ladder_" << ld;
    ladder[ld]=new TH1D(Stream.str().c_str(),(name+title).c_str(),xrange,0,xrange);
     for(int bn=0;bn<LADDER_BIN;++bn)
      ladder[ld]->SetBinContent(bn+1,full_histo->GetBinContent(ld*xrange+bn));
     ladder[ld]->SetStats(0);
  }
  display->cd(1);
  ladder[0]->Draw("same");
  display->cd(4);
  ladder[1]->Draw("same");
  display->cd(2);
  ladder[2]->Draw("same");
  display->cd(5);
  ladder[3]->Draw("same");
  display->cd(3);
  ladder[4]->Draw("same");
  display->cd(6);
  ladder[5]->Draw("same");
  return display;
}

TCanvas *display_ladders2D(TH2D *full_histo,std::string name,std::string title){
  TCanvas *display= new TCanvas();
  display->Divide(3,2);
  int xbin=(full_histo->GetNbinsX())/N_LADDER;
  int ybin=full_histo->GetNbinsY();

  double ymin=full_histo->GetYaxis()->GetXmin();
  double ymax=full_histo->GetYaxis()->GetXmax();  
  
  TH2D *ladder[N_LADDER];
  for(int ld=0;ld<N_LADDER;++ld){
    std::stringstream Stream;
    Stream<<name<<"ladder_" << ld;
    ladder[ld]=new TH2D(Stream.str().c_str(),(name+title).c_str(),xbin,0,xbin,ybin,ymin,ymax);
    for(int xch=0;xch<=xbin;++xch){
      for(int ych=0;ych<ybin;++ych)
	ladder[ld]->SetBinContent(xch,ych,full_histo->GetBinContent(ld*xbin+xch,ych));  
    }
    ladder[ld]->SetStats(0);
  }
  display->cd(1);
  gPad->cd(1);
  gPad->SetLogz();

  ladder[0]->Draw("zcol");
  display->cd(4);
  gPad->cd(4);
  gPad->SetLogz();
  ladder[1]->Draw("zcol");
  display->cd(2);
  gPad->cd(2);
  gPad->SetLogz();
  ladder[2]->Draw("zcol");
  display->cd(5);
  gPad->cd(5);
  gPad->SetLogz();
  ladder[3]->Draw("zcol");
  display->cd(3);
  gPad->cd(3);
  gPad->SetLogz();
  ladder[4]->Draw("zcol");
  display->cd(6);
  gPad->cd(6);
  gPad->SetLogz();
  ladder[5]->Draw("zcol");
  return display;
}


TCanvas* drawing1D(TH1D *histo,bool log){
  TCanvas *display=new TCanvas();
  histo->SetStats(0);
  if (log) gPad->SetLogy();
  histo->Draw("");
  return display;
}

TCanvas* drawing2D(TH2D *histo,std::string title){
  TCanvas *display=new TCanvas();
  histo->SetTitle(title.c_str());
  histo->SetStats(0);
  gPad->SetLogz();
  histo->Draw("zcol");
  return display;
}
TCanvas *drawing6_1D(TH1D *ladder[N_LADDER]){
  TCanvas *display=new TCanvas();
 display->Divide(3,2);
  for(int i=0;i<N_LADDER;++i)
    ladder[i]->SetStats(0);
  display->cd(1);
  ladder[0]->Draw("");
  display->cd(4);
  ladder[1]->Draw("");
  display->cd(2);
  ladder[2]->Draw("");
  display->cd(5);
  ladder[3]->Draw("");
  display->cd(3);
  ladder[4]->Draw("");
  display->cd(6);
  ladder[5]->Draw("");
  return display;
}
TCanvas *drawing6_1D_same(TH1D *ladder[N_LADDER],TH1D *ladder2[N_LADDER]){
  TCanvas *display=new TCanvas();
 display->Divide(3,2);
 for(int i=0;i<N_LADDER;++i){
    ladder[i]->SetStats(0);
    ladder2[i]->SetStats(0);
    ladder2[i]->SetLineColor(2);
 }
  display->cd(1);
  ladder[0]->Draw("");
  ladder2[0]->Draw("same");
  display->cd(4);
  ladder[1]->Draw("");
  ladder2[1]->Draw("same");
  display->cd(2);
  ladder[2]->Draw("");
  ladder2[2]->Draw("same");
  display->cd(5);
  ladder[3]->Draw("");
  ladder2[3]->Draw("same");
  display->cd(3);
  ladder[4]->Draw("");
  ladder2[4]->Draw("same");
  display->cd(6);
  ladder[5]->Draw("");
  ladder2[5]->Draw("same");
  return display;
}
TCanvas *drawing6_2D(TH2D *ladder[N_LADDER]){
  TCanvas *display=new TCanvas();
 display->Divide(3,2);
 for(int i=0;i<N_LADDER;++i){
    ladder[i]->SetStats(0);
    gPad->cd(i+1);
    gPad->SetLogz();
 }
  display->cd(1);
  gPad->cd(1);
  gPad->SetLogz();
  ladder[0]->Draw("zcol");
  display->cd(4);
  gPad->cd(4);
  gPad->SetLogz();
  ladder[1]->Draw("zcol");
  display->cd(2);
  gPad->cd(2);
  gPad->SetLogz();
  ladder[2]->Draw("zcol");
  display->cd(5);
  gPad->cd(5);
  gPad->SetLogz();
  ladder[3]->Draw("zcol");
  display->cd(3);
  gPad->cd(3);
  gPad->SetLogz();
  ladder[4]->Draw("zcol");
  display->cd(6);
  gPad->cd(6);
  gPad->SetLogz();
  ladder[5]->Draw("zcol");
  return display;
}

void PrintPDF(std::vector <std::vector <event>> storage, int MAXEVENTS, int eta_dist_p[N_LADDER][ETASTEP], int eta_dist_n[N_LADDER][ETASTEP], TFile *outputroot,std::string outputname){
  
  std::stringstream Stream_out;
  Stream_out<<outputname<<"_events.txt";
  std:: ofstream outputfile(Stream_out.str().c_str());

  outputroot->cd();

  TH2D *significativit_histo=new TH2D("significativit","",NCHAN,0,NCHAN,200,-50,50.);
  TH2D *clusterseed_histo=new TH2D("clusterseed","",1000,0,1000,700,-200,500);
  TH2D *seed_m1_histo=new TH2D("seed_m1","",1000,0,1000,700,-200,500);
  TH2D *seed_p1_histo=new TH2D("seed_p1","",1000,0,1000,700,-200,500);

  TH2D *clustersigma_histo=new TH2D("clustersigma","",MAX_SIGMA2,0,MAX_SIGMA2,200,-20,100);
  TH2D *clustersigma_m1_histo=new TH2D("clustersigma_m1","",MAX_SIGMA2,0,MAX_SIGMA2,200,-20,100);
  TH2D *clustersigma_p1_histo=new TH2D("clustersigma_p1","",MAX_SIGMA2,0,MAX_SIGMA2,200,-20,100);

  TH2D *clusterseed_corr_histo=new TH2D("clusterseed_corr","",1000,0,1000,700,-200,500);
  TH2D *clusterseed_m1_corr_histo=new TH2D("seed_m1_corr","",1000,0,1000,700,-200,500);
  TH2D *clusterseed_p1_corr_histo=new TH2D("seed_p1_corr","",1000,0,1000,700,-200,500);

  TH1D *deltax_central_ladder_p_hist=new TH1D("deltax_central_ladder_p_hist","deltax_central_ladder_p;deltax;counts",100,-50,50);
  TH1D *deltax_central_ladder_n_hist=new TH1D("deltax_central_ladder_n_hist","deltax_central_ladder_n;deltax;counts",100,-50,50);

TH1D *deltax_central_ladder_noeta_p_hist=new TH1D("deltax_central_ladder_noeta_p_hist","deltax_central_ladder_noeta_p;deltax;counts",200,-5,5);
TH1D *deltax_central_ladder_noeta_n_hist=new TH1D("deltax_central_ladder_noeta_n_hist","deltax_central_ladder_noeta_n;deltax;counts",200,-5,5);

 TH1D *theta_central_ladder_p_hist=new TH1D("theta_central_ladder_p","theta_central_ladder_p;theta;counts",100,-90,90);
 TH1D *theta_central_ladder_n_hist=new TH1D("theta_central_ladder_n","theta_central_ladder_n;theta;counts",100,-90,90);

 TH1D *n3_cluster_place=new TH1D("n3_cluster_position","n3_cluster_position",SIDE_CHAN,0,SIDE_CHAN);

 TH2D *check_timing_right_histo=new TH2D("check_timing_right","check_timing_righ;ADC_seed;(ADC_right-ADC_seed)/ADC_seed",100,0,500,100,-1,2);
 TH2D *check_timing_left_histo=new TH2D("check_timing_left","check_timing_left;ADC_seed;(ADC_left-ADC_seed)/ADC_seed",100,0,500,100,-1,2);

  TH1D *sum_adc_p_hist[N_LADDER];
  TH1D *sum_adc_n_hist[N_LADDER];
  TH1D *sum_adc_p_corr_hist[N_LADDER];
  TH1D *sum_adc_n_corr_hist[N_LADDER];
  TH1D *sum_sign_p_hist[N_LADDER];
  TH1D *sum_sign_n_hist[N_LADDER];

  TH1D *sum_m1_p_hist[N_LADDER];
  TH1D *sum_m1_n_hist[N_LADDER];
  TH1D *sum_p1_p_hist[N_LADDER];
  TH1D *sum_p1_n_hist[N_LADDER];
  TH1D *sum_seed_p_hist[N_LADDER];
  TH1D *sum_seed_n_hist[N_LADDER];

  TH2D *seed_frac_sum_p_histo[N_LADDER];
  TH2D *seed_frac_sum_n_histo[N_LADDER];
  TH2D *seed_frac_m1_p_histo[N_LADDER];
  TH2D *seed_frac_m1_n_histo[N_LADDER];
  TH2D *seed_frac_p1_p_histo[N_LADDER];
  TH2D *seed_frac_p1_n_histo[N_LADDER];

  TH2D *clusterseed_p_histo[N_LADDER];
  TH2D *clusterseed_n_histo[N_LADDER];
  TH2D *clusterseed_m1_p_histo[N_LADDER];
  TH2D *clusterseed_m1_n_histo[N_LADDER];
  TH2D *clusterseed_p1_p_histo[N_LADDER];
  TH2D *clusterseed_p1_n_histo[N_LADDER];

  TH2D *sign_frac_sum_p_histo[N_LADDER];
  TH2D *sign_frac_sum_n_histo[N_LADDER];
  TH2D *sign_frac_m1_p_histo[N_LADDER];
  TH2D *sign_frac_m1_n_histo[N_LADDER];
  TH2D *sign_frac_p1_p_histo[N_LADDER];
  TH2D *sign_frac_p1_n_histo[N_LADDER];

  TH2D *adc_p1_m1_p_histo[N_LADDER];
  TH2D *adc_p1_m1_n_histo[N_LADDER];

  TH1D *charge_center_p_hist[N_LADDER];
  TH1D *charge_center_n_hist[N_LADDER];
  TH1D *charge_center_corr_p_hist[N_LADDER];
  TH1D *charge_center_corr_n_hist[N_LADDER];

  TH1D *clustersize_p_hist[N_LADDER];
  TH1D *clustersize_n_hist[N_LADDER];

  TH1D *eta_first_p_hist[N_LADDER];
  TH1D *eta_first_n_hist[N_LADDER];
  TH2D *eta_ADC_p_histo[N_LADDER];
  TH2D *eta_ADC_n_histo[N_LADDER];
  TH1D *f_eta_p_hist[N_LADDER];
  TH1D *f_eta_n_hist[N_LADDER];

  TH1D *real_cluster_pos_p_hist[N_LADDER];
  TH1D *real_cluster_pos_n_hist[N_LADDER];

  TH2D *cluster_shape_p_histo[N_LADDER];
  TH2D *cluster_shape_n_histo[N_LADDER];
  
  TH1D *real_charge_center_p_hist[N_LADDER];
  TH1D *real_charge_center_n_hist[N_LADDER];
  
  for(int ld=0;ld<N_LADDER;++ld){
    sum_adc_p_hist[ld]=new TH1D(Form("ADC_entries_p_%d",ld),Form("ADC_entries_p_%d;ADC;",ld),500,0,500);
    sum_adc_n_hist[ld]=new TH1D(Form("ADC_entries_n_%d",ld),Form("ADC_entries_n_%d;ADC;",ld),500,0,500);
    sum_adc_p_corr_hist[ld]=new TH1D(Form("ADC_entries_p_corr_%d",ld),Form("ADC_entries_p_corr_%d;ADC;",ld),500,0,500);
    sum_adc_n_corr_hist[ld]=new TH1D(Form("ADC_entries_n_corr_%d",ld),Form("ADC_entries_n_corr_%d;ADC;",ld),500,0,500);
    sum_sign_p_hist[ld]=new TH1D(Form("sign_entries_p_%d",ld),Form("sign_entries_p_%d;ADC;",ld),100,0,100);
    sum_sign_n_hist[ld]=new TH1D(Form("sign_entries_n_%d",ld),Form("sign_entries_n_%d;ADC;",ld),100,0,100);

    sum_m1_p_hist[ld]=new TH1D(Form("sum_m1_p_%d",ld),Form("sum_m1_p_%d;ADC;",ld),500,0,500);
    sum_m1_n_hist[ld]=new TH1D(Form("sum_m1_n_%d",ld),Form("sum_m1_n_%d;ADC;",ld),500,0,500);
    sum_p1_p_hist[ld]=new TH1D(Form("sum_p1_p_%d",ld),Form("sum_p1_p_%d;ADC;",ld),500,0,500);
    sum_p1_n_hist[ld]=new TH1D(Form("sum_p1_n_%d",ld),Form("sum_p1_n_%d;ADC;",ld),500,0,500);
    sum_seed_p_hist[ld]=new TH1D(Form("sum_seed_p_%d",ld),Form("sum_seed_p_%d;ADC;",ld),500,0,500);
    sum_seed_n_hist[ld]=new TH1D(Form("sum_seed_n_%d",ld),Form("sum_seed_n_%d;ADC;",ld),500,0,500);

    seed_frac_sum_p_histo[ld]=new TH2D(Form("seed_frac_sum_p_histo_%d",ld),Form("seed_frac_sum_p_histo_%d;ADC;ADC",ld),1000,0,1000,100,-20,20);
    seed_frac_sum_n_histo[ld]=new TH2D(Form("seed_frac_sum_n_histo_%d",ld),Form("seed_frac_sum_n_histo_%d;ADC;ADC",ld),1000,0,1000,100,-20,20);
    seed_frac_m1_p_histo[ld]=new TH2D(Form("seed_frac_m1__p_histo_%d",ld),Form("seed_frac_m1__p_histo_%d;ADC;ADC",ld),1000,0,1000,100,-20,20);
    seed_frac_m1_n_histo[ld]=new TH2D(Form("seed_frac_m1_n_histo_%d",ld),Form("seed_frac_m1__n_histo_%d;ADC;ADC",ld),1000,0,1000,100,-20,20);
    seed_frac_p1_p_histo[ld]=new TH2D(Form("seed_frac_p1_p_histo_%d",ld),Form("seed_frac_p1_p_histo_%d;ADC;ADC",ld),1000,0,1000,100,-20,20);
    seed_frac_p1_n_histo[ld]=new TH2D(Form("seed_frac_p1_n_histo_%d",ld),Form("seed_frac_p1_n_histo_%d",ld),1000,0,1000,100,-20,20);

    clusterseed_p_histo[ld]=new TH2D(Form("clusterseed_p_%d",ld),Form("clusterseed_p_%d;seed [ADC]; (seed+1)+(seed-1) [ADC]",ld),500,0,500,400,-200,200);
    clusterseed_n_histo[ld]=new TH2D(Form("clusterseed_n_%d",ld),Form("clusterseed_n_%d;seed [ADC]; (seed+1)+(seed-1) [ADC]",ld),500,0,500,400,-200,200);
    clusterseed_m1_p_histo[ld]=new TH2D(Form("clusterseed_m1_p_%d",ld),Form("clusterseed_m1_p_%d;seed [ADC]; seed-1 [ADC]",ld),500,0,500,400,-200,200);
    clusterseed_m1_n_histo[ld]=new TH2D(Form("clusterseed_m1_n_%d",ld),Form("clusterseed_m1_n_%d;seed [ADC]; seed-1 [ADC]",ld),500,0,500,400,-200,200);
    clusterseed_p1_p_histo[ld]=new TH2D(Form("clusterseed_p1_p_%d",ld),Form("clusterseed_p1_p_%d;seed [ADC]; seed+1 [ADC]",ld),500,0,500,400,-200,200);
clusterseed_p1_n_histo[ld]=new TH2D(Form("clusterseed_p1_n_%d",ld),Form("clusterseed_p1_n_%d;seed [ADC]; seed+1 [ADC]",ld),500,0,500,400,-200,200);


    sign_frac_sum_p_histo[ld]=new TH2D(Form("sign_frac_sum_p_histo_%d",ld),Form("sign_frac_sum_p_histo_%d",ld),100,0,100,100,-20,20);
    sign_frac_sum_n_histo[ld]=new TH2D(Form("sign_frac_sum_n_histo_%d",ld),Form("sign_frac_sum_n_histo_%d",ld),100,0,100,100,-20,20);
    sign_frac_m1_p_histo[ld]=new TH2D(Form("sign_frac_m1_p_histo_%d",ld),Form("sign_frac_m1_p_histo_%d",ld),100,0,100,100,-20,20);
    sign_frac_m1_n_histo[ld]=new TH2D(Form("sign_frac_m1_n_histo_%d",ld),Form("sign_frac_m1_n_histo_%d",ld),100,0,100,100,-20,20);
    sign_frac_p1_p_histo[ld]=new TH2D(Form("sign_frac_p1_p_histo_%d",ld),Form("sign_frac_p1_p_histo_%d",ld),100,0,100,100,-20,20);
    sign_frac_p1_n_histo[ld]=new TH2D(Form("sign_frac_p1_n_histo_%d",ld),Form("sign_frac_p1_n_histo_%d",ld),100,0,100,100,-20,20);

    adc_p1_m1_p_histo[ld]=new TH2D(Form("adc_p1_m1_p_%d",ld),Form("adc_p1_m1_p_%d;ADC_p1[ADC];ADC_m1[ADC]",ld),160,-10,150,160,-10,150);
    adc_p1_m1_n_histo[ld]=new TH2D(Form("adc_p1_m1_n_%d",ld),Form("adc_p1_m1_n_%d;ADC_p1[ADC];ADC_m1[ADC]",ld),160,-10,150,160,-10,150);

    charge_center_p_hist[ld]=new TH1D(Form("charge_center_p_%d",ld),Form("charge_center_p_%d;cluster_chan;counts",ld),CLUSTERCHANNELS*8,0,CLUSTERCHANNELS);
    charge_center_n_hist[ld]=new TH1D(Form("charge_center_n_%d",ld),Form("charge_center_n_%d;cluster_chan;counts",ld),CLUSTERCHANNELS*8,0,CLUSTERCHANNELS);
    charge_center_corr_p_hist[ld]=new TH1D(Form("charge_center_corr_p_%d",ld),Form("charge_center_corr_p_%d;cluster_chan;counts",ld),CLUSTERCHANNELS*8,0,CLUSTERCHANNELS);
    charge_center_corr_n_hist[ld]=new TH1D(Form("charge_center_corr_n_%d",ld),Form("charge_center_corr_n_%d;cluster_chan;counts",ld),CLUSTERCHANNELS*8,0,CLUSTERCHANNELS);

    clustersize_p_hist[ld]=new TH1D(Form("clustersize_p_%d",ld),Form("clusersize_p_%d;clustersize;counts",ld),7,0,7);
    clustersize_n_hist[ld]=new TH1D(Form("clustersize_n_%d",ld),Form("clusersize_n_%d;clustersize;counts",ld),7,0,7);

    eta_first_p_hist[ld]=new TH1D(Form("eta_p_%d",ld),Form("eta_p_%d;eta;dN/deta",ld),ETASTEP,ETAMIN,ETAMAX);
    eta_first_n_hist[ld]=new TH1D(Form("eta_n_%d",ld),Form("eta_n_%d;eta;dN/deta",ld),ETASTEP,ETAMIN,ETAMAX);
    eta_ADC_p_histo[ld]=new TH2D(Form("eta_ADC_p_%d",ld),Form("eta_ADC_p_%d;eta;ADC_cluster",ld),ETASTEP,ETAMIN-0.2,ETAMAX+0.2,500,0,500);
    eta_ADC_n_histo[ld]=new TH2D(Form("eta_ADC_n_%d",ld),Form("eta_ADC_n_%d;eta;ADC_cluster",ld),ETASTEP,ETAMIN-0.2,ETAMAX+0.2,500,0,500);
    f_eta_p_hist[ld]=new TH1D(Form("f_eta_p_%d",ld),Form("f_eta_p_%d;eta;f(eta)",ld),ETASTEP,ETAMIN,ETAMAX);
    f_eta_n_hist[ld]=new TH1D(Form("f_eta_n_%d",ld),Form("f_eta_n_%d;eta;f(eta)",ld),ETASTEP,ETAMIN,ETAMAX);

  real_cluster_pos_p_hist[ld]=new TH1D(Form("real_cluster_pos_p_%d",ld),Form("real_cluster_pos_p_%d;chan;counts",ld),SIDE_CHAN*4,0,SIDE_CHAN);
  real_cluster_pos_n_hist[ld]=new TH1D(Form("real_cluster_pos_n_%d",ld),Form("real_cluster_pos_n_%d;chan;counts",ld),SIDE_CHAN*4,0,SIDE_CHAN);

  cluster_shape_p_histo[ld]=new TH2D(Form("cluster_shape_p_%d",ld),Form("cluster_shape_p_%d;position;SN",ld),5,0,5,100,-10,50);
  cluster_shape_n_histo[ld]=new TH2D(Form("cluster_shape_n_%d",ld),Form("cluster_shape_n_%d;position;SN",ld),5,0,5,100,-10,50);
  }

 //eta vectors definition
  double step_eta=(double)ETARANGE/(double)ETASTEP;
  double basement[ETASTEP];
  for(int ch=0;ch<ETASTEP;++ch){
    basement[ch]=ch*step_eta;
  }

  //f_eta calculation
  double f_eta_p[N_LADDER][ETASTEP];
  double f_eta_n[N_LADDER][ETASTEP];
  for(int ladder=0;ladder<N_LADDER;++ladder){
    //diff_p=1.e10;
    f_eta(ETASTEP,basement,eta_dist_p[ladder],f_eta_p[ladder]);
    f_eta(ETASTEP,basement,eta_dist_n[ladder],f_eta_n[ladder]);
    for (int ch=0;ch<ETASTEP;++ch){      
      f_eta_p_hist[ladder]->SetBinContent(ch+1,f_eta_p[ladder][ch]);
      f_eta_n_hist[ladder]->SetBinContent(ch+1,f_eta_n[ladder][ch]);
      //eta_first_p_hist[ladder]->SetBinContent(ch+1,eta_symm_p[ladder][ch]);
      //eta_first_n_hist[ladder]->SetBinContent(ch+1,eta_symm_n[ladder][ch]);
    }
  }

  int pcounts[N_LADDER]={0};
  int ncounts[N_LADDER]={0};
  for(int vec=0;vec<(int)storage.size();++vec){
      for(int nev=0;nev<(int)(storage.at(vec).size());++nev){// I want to correct values for f_eta and fill all the histograms here.
	double x_p[N_LADDER]={0.};
	double x_n[N_LADDER]={0.};
	double sn_p[N_LADDER]={0.};
	double sn_n[N_LADDER]={0.};
	int seed_p[N_LADDER]={0 };
	int seed_n[N_LADDER]={0 };
	 bool pcheck[N_LADDER];
	 bool ncheck[N_LADDER];

	  for(int ild=0;ild<N_LADDER;++ild){
	    pcheck[ild]=0;
	    ncheck[ild]=0;
	  }
	for(int ncl=0;ncl<(int)(storage.at(vec).at(nev).cls.size());++ncl){
	  
	  
	  int ev_seed=storage.at(vec).at(nev).cls.at(ncl).seed;
	  double ev_adc_seed=storage.at(vec).at(nev).cls.at(ncl).count[CLUSTERCHANNELS/2];
	  double ev_adc_p1=storage.at(vec).at(nev).cls.at(ncl).count[CLUSTERCHANNELS/2+1];
	  double ev_adc_m1=storage.at(vec).at(nev).cls.at(ncl).count[CLUSTERCHANNELS/2-1];
	  double ev_sign_seed=storage.at(vec).at(nev).cls.at(ncl).sn[CLUSTERCHANNELS/2];
	  double ev_sign_p1=storage.at(vec).at(nev).cls.at(ncl).sn[CLUSTERCHANNELS/2+1];
	  double ev_sign_m1=storage.at(vec).at(nev).cls.at(ncl).sn[CLUSTERCHANNELS/2-1];

	  
	  double chargecenter=0.;
	  double totcharge=1.;
	  double eta=storage.at(vec).at(nev).cls.at(ncl).GetEta();
	  
	  double sign_cluster=0.;
	  double sign_cluster_num=0.;
	  double sign_cluster_den=0.;
	 
	  for(int cl=0;cl<CLUSTERCHANNELS;++cl){
	    sign_cluster_num+=storage.at(vec).at(nev).cls.at(ncl).count[cl];
	    sign_cluster_den+=storage.at(vec).at(nev).cls.at(ncl).sigma[cl]*storage.at(vec).at(nev).cls.at(ncl).sigma[cl];
	  }
	  sign_cluster=sign_cluster_num/sqrt(sign_cluster_den);


	  
	  int ladder=ChanToLadder(ev_seed);
	  int side=ChanToSide(ev_seed);
	 
	 
	  if(side){//n side
	    if(!ncheck[ladder]){
	      ++ncounts[ladder];
	      ncheck[ladder]=1;
	    }
	    
	    if(eta>=ETAMIN && eta<=ETAMAX ){
	      if(ETAMAX==2)
		chargecenter=CLUSTERCHANNELS/2-0.5+f_eta_n[ladder][(int)((eta/*+shift_f_eta_n[ladder]*/)*(ETASTEP/ETARANGE))];
	      else
		chargecenter=(ev_sign_p1>ev_sign_m1 ? CLUSTERCHANNELS/2+f_eta_n[ladder][(int)(eta*ETASTEP)] : CLUSTERCHANNELS/2-1+f_eta_n[ladder][(int)(eta*ETASTEP)]);
	       //chargecenter=(ev_sign_p1>ev_sign_m1 ? CLUSTERCHANNELS/2+eta : CLUSTERCHANNELS/2-1+eta);
	    }
	    else chargecenter=9999.;
	    
	    //chargecenter=9999.;
	    if(!x_n[ladder] || sign_cluster>sn_n[ladder]){
	    x_n[ladder]=(double)(ev_seed%SIDE_CHAN)+(chargecenter/totcharge-CLUSTERCHANNELS/2);
	    sn_n[ladder]=sign_cluster; //Warning!!!! Should use the squared sum??
	    }
	    
	    seed_n[ladder]=ev_seed%SIDE_CHAN;
	    //outputfile << nev <<"\t"<< ladder <<"\t"<<side<<"\t"<<x_n[ladder] <<std::endl;


	    //Filling histograms
	    charge_center_corr_n_hist[ladder]->Fill(chargecenter/totcharge);


	    clusterseed_n_histo[ladder]->Fill(ev_adc_seed,ev_adc_p1+ev_adc_m1);
	    clusterseed_m1_n_histo[ladder]->Fill(ev_adc_seed,ev_adc_m1);
	    clusterseed_p1_n_histo[ladder]->Fill(ev_adc_seed,ev_adc_p1);
	  
	    //sum_adc_n_hist[ladder]->Fill((ev_adc_seed+ev_adc_p1+ev_adc_m1)*correction(eta));
	    //sum_adc_n_corr_hist[ladder]->Fill(clev.at(nev).cls.at(ncl).GetCounts(3.)*correction(eta));

	    sum_sign_n_hist[ladder]->Fill(sqrt(ev_sign_seed*ev_sign_seed+ev_sign_p1*ev_sign_p1+ev_sign_m1*ev_sign_m1));
	    seed_frac_sum_n_histo[ladder]->Fill(ev_adc_seed,(ev_adc_p1+ev_adc_m1)/ev_adc_seed);
	    seed_frac_m1_n_histo[ladder]->Fill(ev_adc_seed,ev_adc_m1/ev_adc_seed);
	    seed_frac_p1_n_histo[ladder]->Fill(ev_adc_seed,ev_adc_p1/ev_adc_seed);
	    sign_frac_sum_n_histo[ladder]->Fill(ev_sign_seed,sqrt(ev_sign_p1*ev_sign_p1+ev_sign_m1*ev_sign_m1)/ev_sign_seed);
	    sign_frac_m1_n_histo[ladder]->Fill(ev_sign_seed,ev_sign_m1/ev_sign_seed);
	    sign_frac_p1_n_histo[ladder]->Fill(ev_sign_seed,ev_sign_p1/ev_sign_seed);
	    adc_p1_m1_n_histo[ladder]->Fill(ev_adc_p1,ev_adc_m1);
	
	    sum_m1_n_hist[ladder]->Fill(ev_adc_m1);
	    sum_p1_n_hist[ladder]->Fill(ev_adc_p1);
	    sum_seed_n_hist[ladder]->Fill(ev_adc_seed);
	    charge_center_n_hist[ladder]->Fill(storage.at(vec).at(nev).cls.at(ncl).ChargeCenter(3.));
	    clustersize_n_hist[ladder]->Fill(storage.at(vec).at(nev).cls.at(ncl).ClusterSize(3.));
	    eta_first_n_hist[ladder]->Fill(eta);
	    eta_ADC_n_histo[ladder]->Fill(eta,storage.at(vec).at(nev).cls.at(ncl).GetEtaCounts());
	real_cluster_pos_n_hist[ladder]->Fill(ev_seed%SIDE_CHAN+eta*CLUSTERCHANNELS/2);
	for(int i=0;i<CLUSTERCHANNELS;++i)
	  cluster_shape_n_histo[ladder]->Fill(i,storage.at(vec).at(nev).cls.at(ncl).sn[i]);
	

	  }
	  else{//p side
	    if(!pcheck[ladder]){
	      ++pcounts[ladder];
	      pcheck[ladder]=1;
	    }
	    if(eta>=ETAMIN && eta<=ETAMAX ){
	      if(ETAMAX==2)
		chargecenter=CLUSTERCHANNELS/2-0.5+f_eta_p[ladder][(int)((eta/*+shift_f_eta_p[ladder]*/)*(ETASTEP/ETARANGE))];
	      else
		chargecenter=(ev_sign_p1>ev_sign_m1 ? CLUSTERCHANNELS/2+f_eta_p[ladder][(int)(eta*ETASTEP)] : CLUSTERCHANNELS/2-1+f_eta_p[ladder][(int)(eta*ETASTEP)]); //2 chan eta!!!
	      //chargecenter=(ev_sign_p1>ev_sign_m1 ? CLUSTERCHANNELS/2+eta : CLUSTERCHANNELS/2-1+eta);
	      }
	      else chargecenter=9999.;
	    
	    //chargecenter=9999.;
	    if(!x_p[ladder] || sign_cluster>sn_p[ladder]){ // ?????????
	    x_p[ladder]=(double)(ev_seed%SIDE_CHAN)+(chargecenter/totcharge-CLUSTERCHANNELS/2);
	    sn_p[ladder]=sign_cluster; //warning! see n case for more information
	    }
   	    
	    seed_p[ladder]=ev_seed%SIDE_CHAN;
	    //outputfile << nev <<"\t"<< ladder <<"\t"<<side<<"\t"<<x_p[ladder] <<std::endl;

	    //Filling histograms
	    charge_center_corr_p_hist[ladder]->Fill(chargecenter/totcharge);

	    clusterseed_p_histo[ladder]->Fill(ev_adc_seed,ev_adc_p1+ev_adc_m1);
	    clusterseed_m1_p_histo[ladder]->Fill(ev_adc_seed,ev_adc_m1);
	    clusterseed_p1_p_histo[ladder]->Fill(ev_adc_seed,ev_adc_p1);
	  
	    //sum_adc_p_hist[ladder]->Fill((ev_adc_seed+ev_adc_p1+ev_adc_m1)*correction(eta));
	    //sum_adc_p_corr_hist[ladder]->Fill(clev.at(nev).cls.at(ncl).GetCounts(3.)*correction(eta));

	    sum_sign_p_hist[ladder]->Fill(sqrt(ev_sign_seed*ev_sign_seed+ev_sign_p1*ev_sign_p1+ev_sign_m1*ev_sign_m1));
	    seed_frac_sum_p_histo[ladder]->Fill(ev_adc_seed,(ev_adc_p1+ev_adc_m1)/ev_adc_seed);
	    seed_frac_m1_p_histo[ladder]->Fill(ev_adc_seed,ev_adc_m1/ev_adc_seed);
	    seed_frac_p1_p_histo[ladder]->Fill(ev_adc_seed,ev_adc_p1/ev_adc_seed);
	    sign_frac_sum_p_histo[ladder]->Fill(ev_sign_seed,sqrt(ev_sign_p1*ev_sign_p1+ev_sign_m1*ev_sign_m1)/ev_sign_seed);
	    sign_frac_m1_p_histo[ladder]->Fill(ev_sign_seed,ev_sign_m1/ev_sign_seed);
	    sign_frac_p1_p_histo[ladder]->Fill(ev_sign_seed,ev_sign_p1/ev_sign_seed);
	    adc_p1_m1_p_histo[ladder]->Fill(ev_adc_p1,ev_adc_m1);
	
	    sum_m1_p_hist[ladder]->Fill(ev_adc_m1);
	    sum_p1_p_hist[ladder]->Fill(ev_adc_p1);
	    sum_seed_p_hist[ladder]->Fill(ev_adc_seed);
	    charge_center_p_hist[ladder]->Fill(storage.at(vec).at(nev).cls.at(ncl).ChargeCenter(3.));
	clustersize_p_hist[ladder]->Fill(storage.at(vec).at(nev).cls.at(ncl).ClusterSize(3.));
	eta_first_p_hist[ladder]->Fill(eta);
	
	eta_ADC_p_histo[ladder]->Fill(eta,storage.at(vec).at(nev).cls.at(ncl).GetEtaCounts());
	real_cluster_pos_p_hist[ladder]->Fill(ev_seed%SIDE_CHAN+eta*CLUSTERCHANNELS/2);
	for(int i=0;i<CLUSTERCHANNELS;++i)
	  cluster_shape_p_histo[ladder]->Fill(i,storage.at(vec).at(nev).cls.at(ncl).sn[i]);
	    
	  }
	  if(x_p[2]!=0. && x_p[3]!=0.){
	    deltax_central_ladder_p_hist->Fill(x_p[3]-x_p[2]);
	    deltax_central_ladder_noeta_p_hist->Fill(seed_p[3]-seed_p[2]);
	    theta_central_ladder_p_hist->Fill(atan2((x_p[3]-x_p[2])*0.0182,2.)*57.2958);
	  }
	    	    
	  if(x_n[2]!=0. && x_n[3]!=0.){
	    deltax_central_ladder_n_hist->Fill(x_n[3]-x_n[2]);
	    deltax_central_ladder_noeta_n_hist->Fill(seed_n[3]-seed_n[2]);
	    theta_central_ladder_n_hist->Fill(atan2((x_n[3]-x_n[2])*0.0182,2.)*57.2958);
	  }

	  if(ev_sign_seed<MAX_SIGMA2 && ev_sign_m1<MAX_SIGMA2 && ev_sign_p1<MAX_SIGMA2){
	  clustersigma_histo->Fill(ev_sign_seed,sqrt(ev_sign_m1*ev_sign_m1+ev_sign_p1*ev_sign_p1));
	  clustersigma_m1_histo->Fill(ev_sign_seed,ev_sign_m1);
	  clustersigma_p1_histo->Fill(ev_sign_seed,ev_sign_p1);
	}
	


	  significativit_histo->Fill(ev_seed, ev_sign_seed);
	  clusterseed_histo->Fill(ev_adc_seed,ev_adc_p1+ev_adc_m1);
	  seed_m1_histo->Fill(ev_adc_seed,ev_adc_m1);
	  seed_p1_histo->Fill(ev_adc_seed,ev_adc_p1);
	  clusterseed_corr_histo->Fill(ev_adc_seed,storage.at(vec).at(nev).cls.at(ncl).GetSides(3.));
	  
	  check_timing_right_histo->Fill(ev_adc_seed,(ev_adc_seed-ev_adc_p1)/ev_adc_seed);
	  check_timing_left_histo->Fill(ev_adc_seed,(ev_adc_seed-ev_adc_m1)/ev_adc_seed);
	  
	}
	if(x_p[2]!=0. && x_p[3]!=0. && x_n[2]!=0. && x_n[3]!=0.)
	    outputfile << nev <<"\t"<< x_p[2] <<"\t"<<x_p[3]<<"\t"<<x_n[2]<<"\t"<<x_n[3] <<std::endl;
      }
    }
 significativit_histo->Write();
  clusterseed_histo->Write();
  seed_m1_histo->Write();
  seed_p1_histo->Write();
  clustersigma_histo->Write();
  clustersigma_m1_histo->Write();
  clustersigma_p1_histo->Write();
  clusterseed_corr_histo->Write();
  deltax_central_ladder_p_hist->Write();
  deltax_central_ladder_n_hist->Write();
  deltax_central_ladder_noeta_p_hist->Write();
  deltax_central_ladder_noeta_n_hist->Write();
  theta_central_ladder_p_hist->Write();
  theta_central_ladder_n_hist->Write();
  check_timing_left_histo->Write();
  check_timing_right_histo->Write();
  n3_cluster_place->Write();
  for(int ld=0;ld<N_LADDER;++ld){

    sum_adc_p_hist[ld]->Write();
    sum_adc_n_hist[ld]->Write();
    sum_adc_p_corr_hist[ld]->Write();
    sum_adc_n_corr_hist[ld]->Write();
    sum_sign_p_hist[ld]->Write();
    sum_sign_n_hist[ld]->Write();
    clusterseed_p_histo[ld]->Write();
    clusterseed_n_histo[ld]->Write();
    clusterseed_m1_p_histo[ld]->Write();
    clusterseed_m1_n_histo[ld]->Write();
    clusterseed_p1_p_histo[ld]->Write();
    clusterseed_p1_n_histo[ld]->Write();
    adc_p1_m1_p_histo[ld]->Write();
    adc_p1_m1_n_histo[ld]->Write();
    sum_m1_p_hist[ld]->Write();
    sum_p1_p_hist[ld]->Write();
    sum_m1_n_hist[ld]->Write();
    sum_p1_n_hist[ld]->Write();
    sum_seed_p_hist[ld]->Write();
    sum_seed_n_hist[ld]->Write();
    charge_center_p_hist[ld]->Write();
    charge_center_n_hist[ld]->Write();
    clustersize_p_hist[ld]->Write();
    clustersize_n_hist[ld]->Write();
    eta_first_p_hist[ld]->Write();
    eta_first_n_hist[ld]->Write();
    eta_ADC_p_histo[ld]->Write();
    eta_ADC_n_histo[ld]->Write();
    f_eta_p_hist[ld]->Write();
    f_eta_n_hist[ld]->Write();
    real_cluster_pos_p_hist[ld]->Write();
    real_cluster_pos_n_hist[ld]->Write();
    cluster_shape_p_histo[ld]->Write();
    cluster_shape_n_histo[ld]->Write();
    charge_center_corr_p_hist[ld]->Write();
    charge_center_corr_n_hist[ld]->Write(); 
  }

    TCanvas *out=new TCanvas();
  /*
 std::stringstream Stream2;
  std::stringstream Stream2_open;
  std::stringstream Stream2_close;
  Stream2<<outputname <<"_results.pdf";
  Stream2_open<<outputname<<"_results.pdf[";
  Stream2_close<<outputname<<"_results.pdf]";

  out->Print(Stream2_open.str().c_str());
  
  drawing2D(clusterseed_histo,"correlation plot;seed [ADC];(seed+1)+(seed-1) [ADC]")->Print(Stream2.str().c_str());
  drawing2D(clusterseed_corr_histo,"correlation plot corr 3*sigma;seed [ADC];(seed+1)+(seed-1) [ADC]")->Print(Stream2.str().c_str());
  drawing6_2D(clusterseed_p_histo)->Print(Stream2.str().c_str());
  drawing6_2D(clusterseed_n_histo)->Print(Stream2.str().c_str());
  drawing2D(seed_p1_histo,"correlation plot seed+1;seed [ADC];(seed+1) [ADC]")->Print(Stream2.str().c_str());
  drawing6_2D(clusterseed_p1_p_histo)->Print(Stream2.str().c_str());
  drawing6_2D(clusterseed_p1_n_histo)->Print(Stream2.str().c_str());
  drawing2D(seed_m1_histo,"correlation plot seed-1;seed [ADC];(seed-1) [ADC]")->Print(Stream2.str().c_str());
  drawing6_2D(clusterseed_m1_p_histo)->Print(Stream2.str().c_str());
  drawing6_2D(clusterseed_m1_n_histo)->Print(Stream2.str().c_str());
  drawing2D(clustersigma_histo,"correlation plot sign;sign_seed;sqrt(sign_(seed+1)^2+sign_(seed-1)^2)")->Print(Stream2.str().c_str());
  drawing2D(clustersigma_p1_histo,"correlation plot sign seed+1;sign_seed ;sign_(seed+1)")->Print(Stream2.str().c_str());
  drawing2D(clustersigma_m1_histo,"correlation plot sign seed-1;sign_seed ;sign_(seed-1)")->Print(Stream2.str().c_str());
  //drawing6_1D(meanmean_hist,0)->Print(Stream2.str().c_str());
  //display_ladders1D(calibration_comp_mean,"mean2_cal-mean2_data",";;ADC",LADDER_BIN/2)->Print(Stream2.str().c_str());
  //display_ladders1D(calibration_comp_sigma,"sigma3_cal-sigma3_data",";;",LADDER_BIN/2)->Print(Stream2.str().c_str());


  out->Print(Stream2_close.str().c_str());
  
  std::stringstream Stream3;
  std::stringstream Stream3_open;
  std::stringstream Stream3_close;
 
  Stream3 << outputname << "_cluster_analysis.pdf";
  Stream3_open << outputname << "_cluster_analysis.pdf[";
  Stream3_close << outputname << "_cluster_analysis.pdf]";

  out->Print(Stream3_open.str().c_str());
  //  drawing6_1D(sum_adc_p_hist)->Print(Stream3.str().c_str());
  //drawing6_1D(sum_adc_n_hist)->Print(Stream3.str().c_str());
  //drawing6_1D(sum_adc_p_corr_hist)->Print(Stream3.str().c_str());
  //drawing6_1D(sum_adc_n_corr_hist)->Print(Stream3.str().c_str());
  drawing6_1D(sum_seed_p_hist)->Print(Stream3.str().c_str());
  drawing6_1D(sum_seed_n_hist)->Print(Stream3.str().c_str());
  drawing6_1D(sum_m1_p_hist)->Print(Stream3.str().c_str());
  drawing6_1D(sum_m1_n_hist)->Print(Stream3.str().c_str());
  drawing6_1D(sum_p1_p_hist)->Print(Stream3.str().c_str());
  drawing6_1D(sum_p1_n_hist)->Print(Stream3.str().c_str());
  drawing6_1D(sum_sign_p_hist)->Print(Stream3.str().c_str());
  drawing6_1D(sum_sign_n_hist)->Print(Stream3.str().c_str());
  drawing6_2D(seed_frac_sum_p_histo)->Print(Stream3.str().c_str());
  drawing6_2D(seed_frac_sum_n_histo)->Print(Stream3.str().c_str());
  drawing6_2D(seed_frac_m1_p_histo)->Print(Stream3.str().c_str());
  drawing6_2D(seed_frac_m1_n_histo)->Print(Stream3.str().c_str());
  drawing6_2D(seed_frac_p1_p_histo)->Print(Stream3.str().c_str());
  drawing6_2D(seed_frac_p1_n_histo)->Print(Stream3.str().c_str());
  drawing6_2D(sign_frac_sum_p_histo)->Print(Stream3.str().c_str());
  drawing6_2D(sign_frac_sum_n_histo)->Print(Stream3.str().c_str());
  drawing6_2D(sign_frac_m1_p_histo)->Print(Stream3.str().c_str());
  drawing6_2D(sign_frac_m1_n_histo)->Print(Stream3.str().c_str());
  drawing6_2D(sign_frac_p1_p_histo)->Print(Stream3.str().c_str());
  drawing6_2D(sign_frac_p1_n_histo)->Print(Stream3.str().c_str());

  drawing6_2D(adc_p1_m1_p_histo)->Print(Stream3.str().c_str());
  drawing6_2D(adc_p1_m1_n_histo)->Print(Stream3.str().c_str());
  out->Print(Stream3_close.str().c_str());
*/
  std::stringstream Stream4;
  std::stringstream Stream4_open;
  std::stringstream Stream4_close;
  Stream4<<outputname <<"_eta.pdf";
  Stream4_open<<outputname<<"_eta.pdf[";
  Stream4_close<<outputname<<"_eta.pdf]";
  out->Print(Stream4_open.str().c_str());

  drawing6_1D(charge_center_p_hist)->Print(Stream4.str().c_str());
  drawing6_1D(charge_center_n_hist)->Print(Stream4.str().c_str());
  drawing6_1D(clustersize_p_hist)->Print(Stream4.str().c_str());
  drawing6_1D(clustersize_n_hist)->Print(Stream4.str().c_str());
  drawing6_2D(cluster_shape_p_histo)->Print(Stream4.str().c_str());
  drawing6_2D(cluster_shape_n_histo)->Print(Stream4.str().c_str());

  drawing6_2D(eta_ADC_p_histo)->Print(Stream4.str().c_str());
  drawing6_2D(eta_ADC_n_histo)->Print(Stream4.str().c_str());
  drawing6_1D(eta_first_p_hist)->Print(Stream4.str().c_str());
  drawing6_1D(eta_first_n_hist)->Print(Stream4.str().c_str());
  drawing6_1D(f_eta_p_hist)->Print(Stream4.str().c_str());
  drawing6_1D(f_eta_n_hist)->Print(Stream4.str().c_str());

  //drawing6_1D(real_cluster_pos_p_hist)->Print(Stream4.str().c_str());
  //drawing6_1D(real_cluster_pos_n_hist)->Print(Stream4.str().c_str());
  drawing6_1D(charge_center_corr_p_hist)->Print(Stream4.str().c_str());
  drawing6_1D(charge_center_corr_n_hist)->Print(Stream4.str().c_str());
  drawing6_1D_same(charge_center_p_hist,charge_center_corr_p_hist)->Print(Stream4.str().c_str());
  drawing6_1D_same(charge_center_n_hist,charge_center_corr_n_hist)->Print(Stream4.str().c_str());
  //drawing1D(deltax_central_ladder_noeta_p_hist,0)->Print(Stream4.str().c_str());
  //drawing1D(deltax_central_ladder_noeta_n_hist,0)->Print(Stream4.str().c_str());
  drawing1D(deltax_central_ladder_p_hist,0)->Print(Stream4.str().c_str());
  drawing1D(deltax_central_ladder_n_hist,0)->Print(Stream4.str().c_str());
  drawing1D( theta_central_ladder_p_hist,0)->Print(Stream4.str().c_str());
  drawing1D( theta_central_ladder_n_hist,0)->Print(Stream4.str().c_str());
  //drawing1D(n3_cluster_place,0)->Print(Stream4.str().c_str());
  drawing2D( check_timing_right_histo,"check_timing_righ;ADC_seed;(ADC_right-ADC_seed)/ADC_seed")->Print(Stream4.str().c_str());
  drawing2D( check_timing_left_histo,"check_timing_left;ADC_seed;(ADC_left-ADC_seed)/ADC_seed")->Print(Stream4.str().c_str());
  out->Print(Stream4_close.str().c_str());


  for(int ld=0;ld<N_LADDER;++ld){
    std::cout<<"Ladder "<<ld <<" side p: "<<pcounts[ld]<<" Overall efficiency: "<<pcounts[ld]/(double)MAXEVENTS<< " Error: "<<sqrt(pcounts[ld]*(1-pcounts[ld]/(double)MAXEVENTS)) <<std::endl;
    std::cout<<"Ladder "<<ld <<" side n: "<<ncounts[ld]<<" Overall efficiency: "<<ncounts[ld]/(double)MAXEVENTS<< " Error: "<<sqrt(pcounts[ld]*(1-pcounts[ld]/(double)MAXEVENTS)) <<std::endl;

    }
}

void analysis(std::string namefile,std::string calib_file,std::string outputname){
  gROOT->Reset();
  gDirectory->GetList()->Delete();
  /*
   std::ifstream infile (calib_file.c_str(),std::ifstream::in);
  std::string buffer;   //avoiding header
  for(int j=0;j<4;++j)
    std::getline(infile,buffer);
  */
  LTrackerCalibration *cal=LTrackerCalibration::Read(calib_file.c_str());
  
  
  //int nslot=cal->GetNSlots();
  double *mean2_calib= cal->GetPedestal(0);
  double *sigma3_calib= cal->GetSigma(0);
  double *isgaus= cal->GetNGIndex(0);
  bool *CNmask= cal->GetCNMask(0);
  bool IsAlive[NCHAN];
  
  LTrackerMask hotmask=cal->GetMaskOnSigma(0,-999,SIGMAHOT);
  bool *evmask=hotmask.GetBool();
  /*
  double mean2_calib[NCHAN];
  double sigma3_calib[NCHAN];
  double gausindex[NCHAN];
  bool CNmask[NCHAN];
   bool IsAlive[NCHAN];

  
  for(int ichan=0;ichan<NCHAN;++ichan){
    mean2_calib[ichan]=0.;
    sigma3_calib[ichan]=0.;
    gausindex[ichan]=0.;
    CNmask[ichan]=0;
  }
     */
  
     for(int ichan=0;ichan<NCHAN;++ichan){
    // infile >> mean2_calib[ichan] >> sigma3_calib[ichan] >>gausindex[ichan]>> CNmask[ichan];
     if(isgaus[ichan]<=GAUSINDEX_CUT&& sigma3_calib[ichan]<SIGMAHOT) IsAlive[ichan]=1;
     else IsAlive[ichan]=0;
     }
  
  // Opening data file
  LEvRec0File input(namefile.c_str());
  
  LEvRec0 ev;
  input.SetTheEventPointer(ev);
  gStyle->SetPalette(1);
  //gROOT->SetBatch(kTRUE);
  const int MAXEVENTS = input.GetEntries();
  const int N_PKG=MAXEVENTS/NCALIBEVENTS;

  //Histogram defining
  std::stringstream Stream_root;
  Stream_root<<outputname<<".root";
  
  TFile *output2=new TFile(Stream_root.str().c_str(),"RECREATE");
  //definig only histos that are not cluster-related
  TH2D *common_noise_total=new TH2D("common_noise_total","",N_VA,0,N_VA,500,-100,100);
 
  //matrices definition
  short *data[NCHAN];
  for(int ichan=0;ichan<NCHAN;++ichan)
    data[ichan]=new short[NCALIBEVENTS];
  short *data2[NCALIBEVENTS];
  double *CN_matrix_clean[NCALIBEVENTS];
  double *counts_clean[NCALIBEVENTS];
  for(int iev=0;iev<NCALIBEVENTS;++iev){
    data2[iev]=new short[NCHAN];
    CN_matrix_clean[iev]= new double[N_VA];
    counts_clean[iev]=new double[NCHAN];
  }

  //eta vectors definition
  double step_eta=(double)ETARANGE/(double)ETASTEP;
  //double basement[ETASTEP];
  int eta_dist_p[N_LADDER][ETASTEP];
  int eta_dist_n[N_LADDER][ETASTEP];
  for(int ch=0;ch<ETASTEP;++ch){
    //basement[ch]=ch*step_eta;
    for(int ladder=0;ladder<N_LADDER;++ladder){
      eta_dist_p[ladder][ch]=0;
      eta_dist_n[ladder][ch]=0;
    }
  }

  std::stringstream Stream2;
  std::stringstream Stream2_open;
  std::stringstream Stream2_close;
  Stream2<<outputname <<"_results.pdf";
  Stream2_open<<outputname<<"_results.pdf[";
  Stream2_close<<outputname<<"_results.pdf]";


  
  std::vector <std::vector <event> > storage; //good events storage
  for(int ipk=0;ipk<N_PKG;++ipk){//loop on group of events
    //reading the data from file
    for(int iev=0; iev<NCALIBEVENTS; ++iev){
      input.GetEntry(NCALIBEVENTS*ipk+iev);
      for(int ichan=0;ichan<NCHAN;ichan++){
	data[ichan][iev]=ev.strip[ichan];
      }
    }
    for(int iev=0;iev<NCALIBEVENTS;++iev){
      for(int ichan=0;ichan<NCHAN;ichan++)
	data2[iev][ichan]=data[ichan][iev];
    }
    //CN calculation
    for(int evt=0;evt<NCALIBEVENTS;evt++)
      ComputeCN(data2[evt],mean2_calib,CNmask,&CN_matrix_clean[evt][0]);
    
    //clean events
    for(int iev=0;iev<NCALIBEVENTS;++iev){
      for(int ichan=0;ichan<NCHAN;++ichan){
	counts_clean[iev][ichan]=(data[ichan][iev]-mean2_calib[ichan]-CN_matrix_clean[iev][ichan/VA_CHAN]);
      }
    }

    
    std::vector< event > clev; // clever OR clean event 
    for(int iev=0;iev<NCALIBEVENTS;++iev){
      std::vector<LTrackerCluster> *clusters=GetClusters(counts_clean[iev],sigma3_calib,evmask);
      event myevent;
      for(int ev=0;ev<clusters->size();++ev){
	int ladder=ChanToLadder(clusters->at(ev).seed);
	double eta=clusters->at(ev).GetEta();
	if(eta>=ETAMIN && eta<=ETAMAX){
	  if(ChanToSide(clusters->at(ev).seed)) //n cases
	    ++eta_dist_n[ladder][(int)(ETASTEP*(eta-ETAMIN)/ETARANGE)];
	  else 
	    ++eta_dist_p[ladder][(int)(ETASTEP*(eta-ETAMIN)/ETARANGE)];
	}
	if(ChanToLadder(clusters->at(ev).seed)==2 && ChanToSide(clusters->at(ev).seed)==1){ //clusters->at(ev).Dump();
	
	  //n3_cluster_place->Fill(ChanToSideChan(clusters->at(ev).seed));
	}
	myevent.cls.push_back(clusters->at(ev));
      }
      if(myevent.cls.size()>0) clev.push_back(myevent);
    }
    storage.push_back(clev);
    for(int iev=0;iev<NCALIBEVENTS;iev++){
      for(int iva=0;iva<N_VA;iva++){
	common_noise_total->Fill(iva,CN_matrix_clean[iev][iva]);
      }
    }
  }//ipkg ends here!!

  PrintPDF(storage,MAXEVENTS,eta_dist_p,eta_dist_n,output2,outputname);
  
  //common_noise_total->Write();
  
  //TCanvas *out=new TCanvas();
  /*
 std::stringstream Stream2;
  std::stringstream Stream2_open;
  std::stringstream Stream2_close;
  Stream2<<outputname <<"_general_results.pdf";
  Stream2_open<<outputname<<"_general_results.pdf[";
  Stream2_close<<outputname<<"_general_results.pdf]";

  out->Print(Stream2_open.str().c_str());
  display_ladders2D(common_noise_total,"common_noise",";VA;ADC")->Print(Stream2.str().c_str());
  //display_ladders2D(sigma1_histo,"sigma_1",";chan;")->Print(Stream2.str().c_str());
  //display_ladders2D(sigma2_histo,"sigma_2",";chan;")->Print(Stream2.str().c_str());
  //display_ladders2D(sigma3_histo,"sigma_3",";chan;")->Print(Stream2.str().c_str());
  // display_ladders2D(total_counts_adc,"tot_counts",";chan;ADC")->Print(Stream2.str().c_str());
  //display_ladders2D(mean2_histo,"mean_2",";chan;ADC")->Print(Stream2.str().c_str());
  //display_ladders2D(total_counts_mean,"total_counts_mean",";chan;")->Print(Stream.str().c_str());
  //  display_ladders2D(total_counts_clean,"total_counts_clean",";chan;ADC clean")->Print(Stream2.str().c_str());
  display_ladders2D(significativit_histo,"SN",";chan;SN")->Print(Stream2.str().c_str());
  //display_ladders1D(good_chan_hist,"good_channel",";chan;")->Print(Stream2.str().c_str());
  */
  
}

void analysis_ondata(std::string namefile,std::string outputname){
  gROOT->Reset();
  gDirectory->GetList()->Delete();
  gStyle->SetPalette(1);
  
// Opening data file
  LEvRec0File input(namefile.c_str());
  LEvRec0 ev;
  input.SetTheEventPointer(ev);
  const int MAXEVENTS = input.GetEntries();
  const int N_PKG=MAXEVENTS/NCALIBEVENTS;
  //Histogram defining
  std::stringstream Stream_root;
  Stream_root<<outputname<<".root";
  
  TFile *output2=new TFile(Stream_root.str().c_str(),"RECREATE");
  
  
  /*
  TH1D *good_chan_hist=new TH1D("good_chan_hist","",NCHAN,0,NCHAN);
  TH1D *meanmean_hist[N_LADDER];
  for(int ld=0;ld<N_LADDER;++ld)
    meanmean_hist[ld]=new TH1D(Form("Mean_value_ladder_%d",ld),Form("Mean_value_ladder_%d",ld),N_PKG,0,N_PKG);
  
  
  TH2D *mean1_histo=new TH2D("mean1_2d","",NCHAN,0.,NCHAN,3000,0,3000);
  TH2D *sigma1_histo=new TH2D("sigma1_2d","",NCHAN,0,NCHAN,500,0,40.);
  TH2D *mean2_histo=new TH2D("mean2_2d","",NCHAN,0.,NCHAN,3000,0,3000);
  TH2D *sigma2_histo=new TH2D("sigma2_2d","",NCHAN,0,NCHAN,500,0,40.);
  TH2D *sigma3_histo=new TH2D("sigma3_2d","",NCHAN,0.,NCHAN,500,0.,40.);
  */
 
  TH2D *common_noise_total=new TH2D("common_noise_total","",N_VA,0,N_VA,500,-100,100);

  //matrices definition
   short *data[NCHAN];
  for(int ichan=0;ichan<NCHAN;++ichan)
    data[ichan]=new short[NCALIBEVENTS];
  short *data2[NCALIBEVENTS];
  for(int iev=0;iev<NCALIBEVENTS;++iev)
    data2[iev]=new short[NCHAN];

  double *counts_clean[NCALIBEVENTS];
  double *CN_matrix_clean[NCALIBEVENTS];
  for(int iev=0;iev<NCALIBEVENTS;++iev){
    CN_matrix_clean[iev]= new double[N_VA];
    counts_clean[iev]=new double[NCHAN];
  }
 
  double step_eta=(double)ETARANGE/(double)ETASTEP;
    double basement[ETASTEP];
    int eta_dist_p[N_LADDER][ETASTEP];
    int eta_dist_n[N_LADDER][ETASTEP];
    for(int ch=0;ch<ETASTEP;++ch){
      basement[ch]=ch*step_eta;
      for(int ladder=0;ladder<N_LADDER;++ladder){
	eta_dist_p[ladder][ch]=0;
	eta_dist_n[ladder][ch]=0;
      }
    }
    std::stringstream Stream_out;
    Stream_out<<outputname<<"_events.txt";
    std:: ofstream outputfile(Stream_out.str().c_str());

    std::vector <std::vector <event> > storage;  //good events storage

    LTrackerCalibrationManager::GetInstance().LoadRun(namefile.c_str());
     LTrackerCalibration *cal=LTrackerCalibrationManager::GetInstance().Calibrate(NCALIBEVENTS,0);
   
    //processing events

  for(int ipk=0;ipk<N_PKG;++ipk){
    //std::cout<<"Processing events. Step "<<ipk<<std::endl;
    for(int iev=0; iev<NCALIBEVENTS; ++iev){
      input.GetEntry(NCALIBEVENTS*ipk+iev);
      for(int ichan=0;ichan<NCHAN;ichan++){
	data[ichan][iev]=ev.strip[ichan];

      }
    }
    for(int iev=0;iev<NCALIBEVENTS;++iev){
      for(int ichan=0;ichan<NCHAN;ichan++)
	data2[iev][ichan]=data[ichan][iev];
    }


     double *mean2_calib= cal->GetPedestal(ipk);
     double *sigma3_calib= cal->GetSigma(ipk);
     double *isgaus=cal->GetNGIndex(ipk);
     bool *CNmask=cal->GetCNMask(ipk);
     bool IsAlive[NCHAN];

    for(int evt=0;evt<NCALIBEVENTS;evt++){
      ComputeCN(data2[evt],mean2_calib,CNmask,&CN_matrix_clean[evt][0]);
    }

    for(int iev=0;iev<NCALIBEVENTS;++iev){
      for(int ichan=0;ichan<NCHAN;++ichan)
	counts_clean[iev][ichan]=0.;
    }
    
    
    for(int ichan=0;ichan<NCHAN;++ichan){
      if(isgaus[ichan]<=GAUSINDEX_CUT&& sigma3_calib[ichan]<SIGMAHOT) IsAlive[ichan]=1;
      else IsAlive[ichan]=0;	
    }
    for(int iev=0;iev<NCALIBEVENTS;++iev){
      for(int ichan=0;ichan<NCHAN;++ichan){
	counts_clean[iev][ichan]=(data[ichan][iev]-mean2_calib[ichan]-CN_matrix_clean[iev][ichan/VA_CHAN]);
       }
     }
    LTrackerMask hotmask=cal->GetMaskOnSigma(ipk,-999,SIGMAHOT);
    bool *evmask=hotmask.GetBool();
    std::vector< event > clev; // clever OR clean event 
    //std::cout<<"Check"<<std::endl;
    for(int iev=0;iev<NCALIBEVENTS;++iev){
      std::vector<LTrackerCluster> *clusters=GetClusters(counts_clean[iev],sigma3_calib,evmask);
      //std::cout<<"Check. Event N "<<iev<<" Package "<<ipk<<std::endl;
      event myevent;
      for(int ev=0;ev<clusters->size();++ev){
	int ladder=ChanToLadder(clusters->at(ev).seed);
	double eta=clusters->at(ev).GetEta();
	if(eta>=ETAMIN && eta<=ETAMAX){
	  if(ChanToSide(clusters->at(ev).seed)) //n cases
	    ++eta_dist_n[ladder][(int)(ETASTEP*(eta)/ETARANGE)];
	  else 
	    ++eta_dist_p[ladder][(int)(ETASTEP*(eta)/ETARANGE)];
	}
	
	myevent.cls.push_back(clusters->at(ev));
      }
      clev.push_back(myevent);
    }
    storage.push_back(clev);
    //}
    
    for(int iev=0;iev<NCALIBEVENTS;iev++){
      for(int iva=0;iva<N_VA;iva++){
	common_noise_total->Fill(iva,CN_matrix_clean[iev][iva]);
      }
    }
  }//ipkg ends here
   PrintPDF(storage,MAXEVENTS,eta_dist_p,eta_dist_n,output2,outputname);
  output2->cd();
  //mean1_histo->Write();
  //sigma1_histo->Write();
  //mean2_histo->Write();
  //sigma2_histo->Write();
  //sigma3_histo->Write();
  //common_noise_total->Write();
  //total_counts_adc->Write();
  //total_counts_mean->Write();
  //total_counts_clean->Write();
  //good_chan_hist->Write();
  

 
  /*
  std::stringstream Stream2;
  std::stringstream Stream2_open;
  std::stringstream Stream2_close;
  Stream2<<outputname <<"_general_results.pdf";
  Stream2_open<<outputname<<"_general_results.pdf[";
  Stream2_close<<outputname<<"_general_results.pdf]";
  */

  //TCanvas *out=new TCanvas();
  /*
  out->Print(Stream2_open.str().c_str());
  display_ladders2D(common_noise_total,"common_noise",";VA;ADC")->Print(Stream2.str().c_str());
  //display_ladders2D(sigma1_histo,"sigma_1",";chan;")->Print(Stream2.str().c_str());
  //display_ladders2D(sigma2_histo,"sigma_2",";chan;")->Print(Stream2.str().c_str());
  //display_ladders2D(sigma3_histo,"sigma_3",";chan;")->Print(Stream2.str().c_str());
  // display_ladders2D(total_counts_adc,"tot_counts",";chan;ADC")->Print(Stream2.str().c_str());
  //display_ladders2D(mean2_histo,"mean_2",";chan;ADC")->Print(Stream2.str().c_str());
  //display_ladders2D(total_counts_mean,"total_counts_mean",";chan;")->Print(Stream2.str().c_str());
  //display_ladders2D(total_counts_clean,"total_counts_clean",";chan;ADC clean")->Print(Stream2.str().c_str());
  
  //display_ladders1D(good_chan_hist,"good_channel",";chan;")->Print(Stream2.str().c_str());
  
  out->Print(Stream2_close.str().c_str());
  
  */
  }
