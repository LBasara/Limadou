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

#include "LEvent.hh"
#include "LTrackerCluster.hh"
#include "LTrackerTools.hh"

//structures

struct event {
  int entry;
  std::vector<LTrackerCluster> cls;
};

//Getting data from the Limadou TTree

void SetBranchAddresses(TTree *tree, LEvent &ev) {
  tree->SetBranchAddress("strip[4608]",&ev.strip);
  tree->SetBranchAddress("trigger_index", &ev.trigger_index);
  tree->SetBranchAddress("hepd_time", &ev.hepd_time);
  tree->SetBranchAddress("event_index", &ev.event_index);
  tree->SetBranchAddress("event_length", &ev.event_length);
  tree->SetBranchAddress("pmt_high[64]", &ev.pmt_high);
  tree->SetBranchAddress("pmt_low[64]", &ev.pmt_low);
  tree->SetBranchAddress("rate_meter[9]", &ev.rate_meter);
  tree->SetBranchAddress("trigger_flag[64]", &ev.trigger_flag);
  tree->SetBranchAddress("alive_time", &ev.alive_time);
  tree->SetBranchAddress("dead_time", &ev.dead_time);
  return;
}

//mean, sigma, common noise functions

double mean_chan(short *vector, double mean=0,double threshold=0,int length=NCALIBEVENTS){
    double result=0;
  int min,max;
  int count=0;
  if((mean-threshold)<0){
    min=0;
  }else
    min=mean-threshold;
  if((mean+threshold)>NADC || (threshold)==0){
    max=NADC;
  }
  else
    max=mean+threshold;
  for(int i=0;i<length;i++){
    if(vector[i]>min && vector[i]<max){
      result+=vector[i];
      count++;
      }
  }
	
  result/=(count);
 

  return result;
}

double sigma_chan(short *vector, double mean,double threshold=0.,int length=NCALIBEVENTS){
  double sigma=0;
  int min,max;
  int count=0;  
   if((mean-threshold)<0){
    min=0;
  }else
    min=mean-threshold;
  if((mean+threshold)>NADC || threshold==0){
    max=NADC;
  }
  else
    max=mean+threshold;

  for(int i=0;i<length;i++){
    if(vector[i]<max && vector[i]>min){
      sigma+=(vector[i]-mean)*(vector[i]-mean);
      count++;
    }
  }
  sigma/=count;
  sigma=sqrt(sigma);
  return sigma;
}

double sigma_clean(short *vector,double *cmn,double mean,double threshold,double ev){
  double result=0.;
  int min=999.;
  int max=999.;
  int count=0;
   if((mean-threshold)<0){
    min=0;
  }else
    min=mean-threshold;
  if((mean+threshold)>NADC){
    max=NADC;
  }
  else
    max=mean+threshold;
  for(int i=0;i<ev;++i){
    if(vector[i]<max && vector[i]>min){
      result+=(vector[i]-mean-cmn[i])*(vector[i]-mean-cmn[i]);
      ++count;
    }
      
  }
  result/=(double)count;
  result=sqrt(result);
  //std::cout<<count<<"\t"<<result<<std::endl;
  return result;

}


//R.I.
double GetCleanedSigma(int *h_va /*sigma_VA*/,int length=BIN_SIGMA) {
    // Put the histo into something workable on
    int nb=length;

    double max=MIN_NCHANSPERBIN_4CN-1.; // minimum number of chans we want for CN calculation: careful! If split into more bins you cannot demand too much...
    double xmax=99999.;
    for(int ib=0; ib<nb; ++ib) {
        double x=ib*SIGMA_STEP;
        double y=h_va[ib];
           
        if(x<=MIN_SIGMA_CHANNEL_OFF || x>=MAX_SIGMA_CHANNEL_OFF) continue; // dead channels. Go ahead.
        if(x>xmax+MAX_SIGMA_NOISE_LEVEL ) break; // you have gone too far. The first 'structure' is the one you are interested in
        if(y>max) { // update max infos
            max=y;
            xmax=x;
        }
    }
    
    return xmax;
} // it returns 99999. by default 


void common_noise_clean(short *adc,double *mean,double *sigma_chan,double *sigma_mean,double *cmn,bool *IsAlive){
  for(int i=0;i<N_VA;++i){
    double result = 0.;
    
    int count = 0;
    double min=-999.;
    double max=-999.;
    double min_s=-999.;
    double max_s=-999.;
    for(int j=0;j<VA_CHAN;++j){
      if(!IsAlive[i*VA_CHAN+j]) continue;
      min=mean[i*VA_CHAN+j]-SIGMA_CUT_CN*sigma_chan[i*VA_CHAN+j];
      max=mean[i*VA_CHAN+j]+SIGMA_CUT_CN*sigma_chan[i*VA_CHAN+j];
      min_s=sigma_mean[i]-SIGMA_SEL;
      max_s=sigma_mean[i]+SIGMA_SEL;				    
      if(adc[i*VA_CHAN+j]<=max && adc[i*VA_CHAN+j]>=min){
	if(sigma_chan[i*VA_CHAN+j]<=max_s && sigma_chan[i*VA_CHAN+j]>=min_s && sigma_chan[i*VA_CHAN+j]>=MIN_SIGMA_CHANNEL_OFF && sigma_chan[i*VA_CHAN+j]<=MAX_SIGMA_CHANNEL_OFF){
	  result += (adc[i*VA_CHAN+j]-mean[i*VA_CHAN+j]);
	  ++count;
	}
      }
    }
    if(count==0)
      result=0.;
    else
      result/=(double)count;
    cmn[i]=result;
  }
  
  return;
}

//Check for channel "gaussianity"
double IsGauss(short *data,double mean,double sigma,double *cn,int nev){
  int count=0;
  for(int iev=0;iev<nev;++iev)
    if(abs(data[iev]-mean-cn[iev])/sigma>NSIGMABOOL) ++count;
  return (count-GAUSSTHRESHOLD*nev)/sqrt(GAUSSTHRESHOLD*nev);
}

//Graphics
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


//locate chan
double correction(double eta){
  double result=-9999.;
  double sum=71.54;
  result=7.909/sum*eta*eta+
    0.774/sum*eta+
    63.64/sum;
  
  return 1/result;
}

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

//Production of calibration file
void calibration(std::string namefile_calib,std::string outputname){
  TFile *input_calib=TFile::Open(namefile_calib.c_str());
  //calibration run processing...
  std::cout<<"Calibration run processing..."<<std::endl;
 
  TTree *tree_calib=(TTree*)input_calib->Get("T");
  LEvent ev_calib;
  SetBranchAddresses(tree_calib,ev_calib);
  tree_calib->SetBranchStatus("*",kFALSE);
  tree_calib->SetBranchStatus("strip[4608]",kTRUE);

  const int NCALIBRUN=tree_calib->GetEntries();  
  short *data_calib[NCHAN];

  for(int ichan=0;ichan<NCHAN;++ichan)
    data_calib[ichan]=new short[NCALIBRUN];
  

  short *data2_calib[NCALIBRUN];
  for(int iev=0;iev<NCALIBRUN;++iev)
    data2_calib[iev]=new short[NCHAN];
  double *CN_matrix_calib[NCALIBRUN];
  for(int iev=0;iev<NCALIBRUN;++iev)
    CN_matrix_calib[iev]= new double[N_VA];

  for(int iev=0; iev<NCALIBRUN; ++iev){
    tree_calib->GetEntry(iev);
    for(int ichan=0;ichan<NCHAN;ichan++){
      data_calib[ichan][iev]=ev_calib.strip[ichan];
    }
  }
  for(int iev=0;iev<NCALIBRUN;++iev){
    for(int ichan=0;ichan<NCHAN;ichan++){
      data2_calib[iev][ichan]=data_calib[ichan][iev];
      
    }
  }


  double mean1_calib[NCHAN];
  double sigma1_calib[NCHAN];
  double mean2_calib[NCHAN];
  double sigma2_calib[NCHAN];
  double sigma3_calib[NCHAN];
  bool IsAlive[NCHAN];

  double co_no_calib[NCALIBRUN];


for(int iev=0;iev<NCALIBRUN;++iev)
    co_no_calib[iev]=0.;
  
  for(int ichan=0;ichan<NCHAN;++ichan){
    mean1_calib[ichan]=9999.;
    sigma1_calib[ichan]=0.;
    mean2_calib[ichan]=0.;
    sigma2_calib[ichan]=0.;
    sigma3_calib[ichan]=0;
    IsAlive[ichan]=1;
  }
  for(int ichan=0;ichan<NCHAN;ichan++){
    mean1_calib[ichan]=mean_chan(data_calib[ichan],0.,0.,NCALIBRUN);
    sigma1_calib[ichan]=sigma_chan(data_calib[ichan],mean1_calib[ichan],0.,NCALIBRUN);
    mean2_calib[ichan]=mean_chan(data_calib[ichan],mean1_calib[ichan],(SIGMA_CUT_SIGMA*sigma1_calib[ichan]),NCALIBRUN);
    sigma2_calib[ichan]=sigma_clean(data_calib[ichan],co_no_calib,mean2_calib[ichan],(SIGMA_CUT_SIGMA*sigma1_calib[ichan]),NCALIBRUN);
  }
    
    //Calculation of the mean value of the distribution of the sigma2 for each channel
    int sigma_dist[BIN_SIGMA];
    double sigma_mean[N_VA];
    for(int iva=0;iva<N_VA;++iva)
      sigma_mean[iva]=0.;
    
    for(int iva=0;iva<N_VA;++iva){
      for(int ibin=0;ibin<BIN_SIGMA;++ibin)
	sigma_dist[ibin]=0;

      for(int ich=0;ich<VA_CHAN;++ich){
	if(sigma2_calib[(iva*VA_CHAN+ich)]<100.)
	  sigma_dist[(int)(sigma2_calib[(iva*VA_CHAN+ich)]/SIGMA_STEP)]++;	 
      }
      sigma_mean[iva]=GetCleanedSigma(sigma_dist);
    }
    
    double isgaus[NCHAN];

  
    
    for(int evt=0;evt<NCALIBRUN;evt++)
    common_noise_clean(data2_calib[evt],mean2_calib,sigma2_calib,sigma_mean,&CN_matrix_calib[evt][0],IsAlive);
      
  
  for(int iva=0;iva<N_VA;++iva){
    for(int ievt=0;ievt<NCALIBRUN;++ievt){
      co_no_calib[ievt]=CN_matrix_calib[ievt][iva];
      
    }
    for(int j=0;j<VA_CHAN;j++){
      sigma3_calib[iva*VA_CHAN+j]=sigma_clean(data_calib[iva*VA_CHAN+j],co_no_calib,mean2_calib[iva*VA_CHAN+j],(SIGMA_CUT_SIGMA*sigma2_calib[iva*VA_CHAN+j]),NCALIBRUN);
      
      //Is the channel gaussian??
      isgaus[iva*VA_CHAN+j]=IsGauss(data_calib[iva*VA_CHAN+j],mean2_calib[iva*VA_CHAN+j],sigma3_calib[iva*VA_CHAN+j],co_no_calib,NCALIBRUN);
     
    }
  }

 std:: ofstream outfile (outputname.c_str(),std::ofstream::out);
     outfile << "***************  LIMADOU calibration file  **************"<<std::endl;
     outfile<< "From file: " <<namefile_calib<<std::endl
	    <<  "*********************************************************"<<std::endl
	    <<"N_chan \t Pedestal \t Sigma2 \t Sigma3 \t Gaus_index"<<std::endl;
     for(int ichan=0;ichan<NCHAN;++ichan)
       outfile << ichan <<"\t"<<mean2_calib[ichan]<<"\t"<<sigma2_calib[ichan]<<"\t"<<sigma3_calib[ichan]<<"\t"<<isgaus[ichan]<<std::endl;
     /*
     TH1D *gaus_index=new TH1D("gausindex","gausindex;chan;gausindex",NCHAN,0,NCHAN);

     for(int ichan=0;ichan<NCHAN;++ichan)
       gaus_index->Fill(ichan,isgaus[ichan]);

     gaus_index->Draw();
     */
}


void analysis(std::string namefile,std::string calib_file,std::string outputname){
  gROOT->Reset();
  gDirectory->GetList()->Delete();

  //reading and storing calibration values
  std::ifstream infile (calib_file.c_str(),std::ifstream::in);
  std::string buffer;   //avoiding header
  for(int j=0;j<4;++j)
    std::getline(infile,buffer);

  double mean2_calib[NCHAN];
  double sigma2_calib[NCHAN];
  double sigma3_calib[NCHAN];
  double gausindex[NCHAN];
  
  for(int ichan=0;ichan<NCHAN;++ichan){
    mean2_calib[ichan]=0.;
    sigma2_calib[ichan]=0.;
    sigma3_calib[ichan]=0;
    gausindex[ichan]=0.;
  }
  int buf;
 
  for(int ichan=0;ichan<NCHAN;++ichan){
      infile >>  buf >> mean2_calib[ichan] >>sigma2_calib[ichan] >> sigma3_calib[ichan]>>gausindex[ichan];
  }
  //Calculation of the mean value of the distribution of the sigma2 for each channel
  int sigma_dist[BIN_SIGMA];
  double sigma_mean[N_VA];
  for(int iva=0;iva<N_VA;++iva)
    sigma_mean[iva]=0.;
    
  for(int iva=0;iva<N_VA;++iva){
    for(int ibin=0;ibin<BIN_SIGMA;++ibin)
      sigma_dist[ibin]=0;
    
    for(int ich=0;ich<VA_CHAN;++ich){
      if(sigma2_calib[(iva*VA_CHAN+ich)]<MAX_SIGMA)
	sigma_dist[(int)(sigma2_calib[(iva*VA_CHAN+ich)]/SIGMA_STEP)]++;
    }
    sigma_mean[iva]=GetCleanedSigma(sigma_dist);
  }
  bool IsAlive[NCHAN];
  for(int ichan=0;ichan<NCHAN;++ichan){
    if(gausindex[ichan]<=GAUSINDEX_CUT&& sigma3_calib[ichan]<SIGMAHOT) IsAlive[ichan]=1;
    else IsAlive[ichan]=0;
  }
  //Opening data File
  TFile *input=TFile::Open(namefile.c_str());
  TTree *tree=(TTree*)input->Get("T");
  LEvent ev;
  gStyle->SetPalette(1);
  //gROOT->SetBatch(kTRUE);
  SetBranchAddresses(tree,ev);
  tree->SetBranchStatus("*",kFALSE);
  tree->SetBranchStatus("strip[4608]",kTRUE);
  const int MAXEVENTS = tree->GetEntries();
  const int N_PKG=MAXEVENTS/NCALIBEVENTS;

  //Histogram defining
  std::stringstream Stream_root;
  Stream_root<<outputname<<".root";
  
  TFile *output2=new TFile(Stream_root.str().c_str(),"RECREATE");
  
  
  
  TH1D *good_chan_hist=new TH1D("good_chan_hist","",NCHAN,0,NCHAN);
  TH1D *meanmean_hist[N_LADDER];
  for(int ld=0;ld<N_LADDER;++ld)
    meanmean_hist[ld]=new TH1D(Form("Mean_value_ladder_%d",ld),Form("Mean_value_ladder_%d",ld),N_PKG,0,N_PKG);
  
  
  TH2D *mean1_histo=new TH2D("mean1_2d","",NCHAN,0.,NCHAN,3000,0,3000);
  TH2D *sigma1_histo=new TH2D("sigma1_2d","",NCHAN,0,NCHAN,500,0,40.);
  TH2D *mean2_histo=new TH2D("mean2_2d","",NCHAN,0.,NCHAN,3000,0,3000);
  TH2D *sigma2_histo=new TH2D("sigma2_2d","",NCHAN,0,NCHAN,500,0,40.);
  TH2D *sigma3_histo=new TH2D("sigma3_2d","",NCHAN,0.,NCHAN,500,0.,40.);

 
  TH2D *common_noise_total=new TH2D("common_noise_total","",N_VA,0,N_VA,500,-100,100);

  TH2D *total_counts_adc=new TH2D("total_counts_adc","total_counts_adc",NCHAN,0,NCHAN,3000,0,3000);
  TH2D *total_counts_mean=new TH2D("total_counts_mean","total_counts_mean",NCHAN,0,NCHAN,1000,-500,500);
  TH2D *total_counts_clean=new TH2D("total_counts_clear","total_counts_clear",NCHAN,0,NCHAN,3000,-500,500);
 

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

  TH1D *deltax_central_ladder_p_hist=new TH1D("deltax_central_ladder_p_hist","deltax_central_ladder_p;deltax;counts",200,-5,5);
TH1D *deltax_central_ladder_n_hist=new TH1D("deltax_central_ladder_n_hist","deltax_central_ladder_n;deltax;counts",200,-5,5);

TH1D *deltax_central_ladder_noeta_p_hist=new TH1D("deltax_central_ladder_noeta_p_hist","deltax_central_ladder_noeta_p;deltax;counts",200,-5,5);
TH1D *deltax_central_ladder_noeta_n_hist=new TH1D("deltax_central_ladder_noeta_n_hist","deltax_central_ladder_noeta_n;deltax;counts",200,-5,5);
  
  

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

    clusterseed_p_histo[ld]=new TH2D(Form("clusterseed_p_%d",ld),Form("clusterseed_p_%d",ld),500,0,500,400,-200,200);
    clusterseed_n_histo[ld]=new TH2D(Form("clusterseed_n_%d",ld),Form("clusterseed_n_%d",ld),500,0,500,400,-200,200);
    clusterseed_m1_p_histo[ld]=new TH2D(Form("clusterseed_m1_p_%d",ld),Form("clusterseed_m1_p_%d",ld),500,0,500,400,-200,200);
    clusterseed_m1_n_histo[ld]=new TH2D(Form("clusterseed_m1_n_%d",ld),Form("clusterseed_m1_n_%d",ld),500,0,500,400,-200,200);
    clusterseed_p1_p_histo[ld]=new TH2D(Form("clusterseed_p1_p_%d",ld),Form("clusterseed_p1_p_%d",ld),500,0,500,400,-200,200);
clusterseed_p1_n_histo[ld]=new TH2D(Form("clusterseed_p1_n_%d",ld),Form("clusterseed_p1_n_%d",ld),500,0,500,400,-200,200);


    sign_frac_sum_p_histo[ld]=new TH2D(Form("sign_frac_sum_p_histo_%d",ld),Form("sign_frac_sum_p_histo_%d",ld),100,0,100,100,-20,20);
    sign_frac_sum_n_histo[ld]=new TH2D(Form("sign_frac_sum_n_histo_%d",ld),Form("sign_frac_sum_n_histo_%d",ld),100,0,100,100,-20,20);
    sign_frac_m1_p_histo[ld]=new TH2D(Form("sign_frac_m1_p_histo_%d",ld),Form("sign_frac_m1_p_histo_%d",ld),100,0,100,100,-20,20);
    sign_frac_m1_n_histo[ld]=new TH2D(Form("sign_frac_m1_n_histo_%d",ld),Form("sign_frac_m1_n_histo_%d",ld),100,0,100,100,-20,20);
    sign_frac_p1_p_histo[ld]=new TH2D(Form("sign_frac_p1_p_histo_%d",ld),Form("sign_frac_p1_p_histo_%d",ld),100,0,100,100,-20,20);
    sign_frac_p1_n_histo[ld]=new TH2D(Form("sign_frac_p1_n_histo_%d",ld),Form("sign_frac_p1_n_histo_%d",ld),100,0,100,100,-20,20);

    adc_p1_m1_p_histo[ld]=new TH2D(Form("adc_p1_m1_p_%d",ld),Form("adc_p1_m1_p_%d",ld),160,-10,150,160,-10,150);
    adc_p1_m1_n_histo[ld]=new TH2D(Form("adc_p1_m1_n_%d",ld),Form("adc_p1_m1_n_%d",ld),160,-10,150,160,-10,150);

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

  TH1D *calibration_comp_mean=new TH1D("mean_comparison","mean_comparison;chan;ADC",NCHAN,0,NCHAN);
  TH1D *calibration_comp_sigma=new TH1D("sigma_comparison","sigma_comparison;chan;",NCHAN,0,NCHAN);
  
  double *counts_clean[NCHAN];
  for(int ichan=0;ichan<NCHAN;++ichan)
    counts_clean[ichan]=new double[NCALIBEVENTS];
  double *CN_matrix_clean[NCALIBEVENTS];
  for(int iev=0;iev<NCALIBEVENTS;++iev)
    CN_matrix_clean[iev]= new double[N_VA];

  double mean2_data[NCHAN];
  double sigma3_data[NCHAN];
  
  for(int ichan=0;ichan<NCHAN;++ichan){
    mean2_data[ichan]=0.;
    sigma3_data[ichan]=0.;
  }

  short *data[NCHAN];
  for(int ichan=0;ichan<NCHAN;++ichan)
    data[ichan]=new short[NCALIBEVENTS];
  short *data2[NCALIBEVENTS];
  for(int iev=0;iev<NCALIBEVENTS;++iev)
    data2[iev]=new short[NCHAN];

  double step_eta=(double)ETARANGE/(double)ETASTEP;
    double basement[ETASTEP];
    int eta_dist_p[N_LADDER][ETASTEP];
    int eta_dist_n[N_LADDER][ETASTEP];
    for(int ch=0;ch<ETASTEP;++ch){
      basement[ch]=ch*step_eta;
      //std::cout<<ch<<"\t"<<ch*step_eta-1<<std::endl;
      for(int ladder=0;ladder<N_LADDER;++ladder){
	eta_dist_p[ladder][ch]=0;
	eta_dist_n[ladder][ch]=0;
      }
    }
    std::stringstream Stream_out;
    Stream_out<<outputname<<"_events.txt";
    std:: ofstream outputfile(Stream_out.str().c_str());

    std::vector <std::vector <event> > storage;  //good events storage
    //processing events
  for(int ipk=0;ipk<N_PKG;++ipk){
    //std::cout<<"Processing events. Step "<<ipk<<std::endl;
    for(int iev=0; iev<NCALIBEVENTS; ++iev){
      tree->GetEntry(NCALIBEVENTS*ipk+iev);
      for(int ichan=0;ichan<NCHAN;ichan++){
	data[ichan][iev]=ev.strip[ichan];
	//total_counts_adc->Fill(ichan,data[ichan][iev]);
      }
    }
    for(int iev=0;iev<NCALIBEVENTS;++iev){
      for(int ichan=0;ichan<NCHAN;ichan++)
	data2[iev][ichan]=data[ichan][iev];
    }
    //common noise calculation
    double co_no[NCALIBEVENTS];
    for(int iev=0;iev<NCALIBEVENTS;++iev)
      co_no[iev]=0.;

    for(int evt=0;evt<NCALIBEVENTS;evt++){
      common_noise_clean(data2[evt],mean2_calib,sigma2_calib,sigma_mean,&CN_matrix_clean[evt][0],IsAlive);
    }

    for(int iev=0;iev<NCALIBEVENTS;++iev){
      for(int ichan=0;ichan<NCHAN;++ichan)
	counts_clean[ichan][iev]=0.;
    }
    

    for(int ichan=0;ichan<NCHAN;++ichan){
      for(int iev=0;iev<NCALIBEVENTS;++iev){
	counts_clean[ichan][iev]=(data[ichan][iev]-mean2_calib[ichan]-CN_matrix_clean[iev][(int)ichan/VA_CHAN]);  
	total_counts_clean->Fill(ichan,counts_clean[ichan][iev]);
	//total_counts_mean->Fill(ichan,data[ichan][iev]-mean2_calib[ichan]);
       }
     }
     
    /*
    //Calculation of sigma without common noise.
     for(int iva=0;iva<N_VA;++iva){
       for(int ievt=0;ievt<NCALIBEVENTS;++ievt)
	 co_no[ievt]=CN_matrix_clean[ievt][iva];
       for(int j=0;j<VA_CHAN;j++){
	 sigma3_data[iva*VA_CHAN+j]=sigma_clean(data[iva*VA_CHAN+j],co_no,mean2_calib[iva*VA_CHAN+j],(SIGMA_CUT_SIGMA*sigma2_calib[iva*VA_CHAN+j]),NCALIBEVENTS);
	 //sigma3_data[iva*VA_CHAN+j]+=sigma3_data[iva*VA_CHAN+j];
       }
     }
    */
    std::vector< event > clev; // clever OR clean event 
    for(int iev=0;iev<NCALIBEVENTS;++iev){
<<<<<<< HEAD
      std::vector<LTrackerCluster>* clusters=GetClusters(counts_clean[iev],sigma3_calib);
      //std::cout<<"Check. Event N "<<iev<<" Package "<<ipk<<std::endl;
      event myevent;
      for(int ev=0;ev<clusters->size();++ev){
	myevent.cls.push_back(clusters->at(ev));
	int ladder=ChanToLadder(clusters->at(ev).seed);
	double eta=clusters->at(ev).GetEta();
      if(eta>=ETAMIN && eta<=ETAMAX){
      	if(ChanToSide(clusters->at(ev).seed)) //n cases
	   ++eta_dist_n[ladder][(int)(ETASTEP*(eta)/ETARANGE)];	 
        else 
	  ++eta_dist_p[ladder][(int)(ETASTEP*(eta)/ETARANGE)];
	}
      }
      
      clev.push_back(myevent);
      
      /*
=======
>>>>>>> f3ec6811eba387071641107fc1dcc97148ba294c
      event myevent;
      myevent.entry = ipk*NCALIBEVENTS+iev;
      for(int ichan=0;ichan<NCHAN;++ichan){
	if (!IsAlive[ichan]) continue;
	double count=counts_clean[ichan][iev];
	double sigma=sigma3_calib[ichan];
	double sign= count/sigma;
	if(sign<GE_THRESHOLD) continue;	
	good_chan_hist->Fill(ichan);
	// Open window for cluster identification *******
	std::vector< double > window_count;
	std::vector< double > window_sigma;
	std::vector< double > window_sign;
	for(int iw=0; iw<NEV_WINDOW; ++iw) {
	  
	  double mycount = counts_clean[ichan+iw][iev];//(data[ichan+iw][iev]-mean2_calib[ichan+iw]-CN_matrix_clean[iev][(ichan+iw)/VA_CHAN]);
	  double mysigma = sigma3_data[ichan+iw];
	  window_count.push_back(mycount);
	  window_sigma.push_back(mysigma);
	  window_sign.push_back(mysigma!=0. ? mycount/mysigma : 0.);
	  if((ichan+iw+1)%384==0) break;
	}
	int wsize=window_count.size();
	//***********************************************
	int maxind=-999;
	double maxsign=-999.;
	for(int il=0; il<wsize; ++il) {
	  if(window_sign.at(il)>maxsign) {
	    maxind=il;
	    maxsign=window_sign.at(il);
	  }
	}

	int max2ind=999;
	double max2sign=-999.;
	for(int il=0; il<wsize; ++il) {
	  if(window_sign.at(il)>max2sign && il!= maxind) {
	    max2ind=il;
	    max2sign=window_sign.at(il);
	  }
	}

	// Seed your life... *****************************
	/*
	int seed = -999;
	
	// Case for dead channel in between
	if((maxind-max2ind==2 || max2ind-maxind==2) && max2sign>3. )seed=ichan+(maxind+max2ind)/2;//WARNING!  
	else if(abs(maxind-max2ind)==1 ) seed=maxind+ichan;
	else seed=((max2ind<maxind && max2sign>GE_THRESHOLD) ? max2ind+ichan : maxind+ichan);
	*/
	int seed=maxind+ichan;

	// ***********************************************
	// Fill cluster data *****************************

	LTrackerCluster mycluster;
	mycluster.seed = seed;
	for(int iii=0; iii<CLUSTERCHANNELS; ++iii) {
	  int cchan = seed-CLUSTERCHANNELS/2+iii;// remember that the shift depends on clustersize!!
	  bool SameLP_FLAG = SameLadderPlane(seed, cchan);
	  if((IsAlive) && SameLP_FLAG ){
	    //if(counts_clean[cchan][iev]/sigma3_calib[cchan]>3.){
	      mycluster.count[iii] =  counts_clean[cchan][iev];
	      mycluster.sigma[iii] =  sigma3_calib[cchan];
	      mycluster.sn[iii] = counts_clean[cchan][iev]/sigma3_calib[cchan];
	      //}
	  }
	  else{//Is the best choice???
	    mycluster.count[iii] = 0.;
	    mycluster.sigma[iii] = 0.;
	    mycluster.sn[iii] = 0.;
	    
	  }

	}
	//oldeta
	//mycluster.eta=(mycluster.count[3]>mycluster.count[1] ? (mycluster.count[3]-mycluster.count[2])/(mycluster.count[3]+mycluster.count[2]) : (mycluster.count[2]-mycluster.count[1])/(mycluster.count[2]+mycluster.count[1]));

	//symmethrize eta!
	if(mycluster.sn[CLUSTERCHANNELS/2+1]>0. || mycluster.sn[CLUSTERCHANNELS/2-1]>0.)
	  mycluster.eta=(mycluster.count[CLUSTERCHANNELS/2+1]>mycluster.count[CLUSTERCHANNELS/2-1] ? (mycluster.count[CLUSTERCHANNELS/2+1])/(mycluster.count[CLUSTERCHANNELS/2+1]+mycluster.count[CLUSTERCHANNELS/2]) : (mycluster.count[CLUSTERCHANNELS/2])/(mycluster.count[CLUSTERCHANNELS/2]+mycluster.count[CLUSTERCHANNELS/2-1]));
	else  mycluster.eta=-999.;	
	myevent.cls.push_back(mycluster);
	ichan = seed+(CLUSTERCHANNELS); // very important! Restart after the cluster
      }
      
      if(myevent.cls.size()>0) clev.push_back(myevent);
    }
    storage.push_back(clev);
    //}
    
    
    for(int nev=0;nev<(int)(clev.size());++nev){
      //Filling the histos and file
      for(int ncl=0;ncl<(int)(clev.at(nev).cls.size());++ncl){
	int ev_seed=clev.at(nev).cls.at(ncl).seed;
	double ev_adc_seed=clev.at(nev).cls.at(ncl).count[2];
	double ev_adc_p1=clev.at(nev).cls.at(ncl).count[3];
	double ev_adc_m1=clev.at(nev).cls.at(ncl).count[1];
	double ev_sign_seed=clev.at(nev).cls.at(ncl).sn[2];
	double ev_sign_p1=clev.at(nev).cls.at(ncl).sn[3];
	double ev_sign_m1=clev.at(nev).cls.at(ncl).sn[1];
	double eta=clev.at(nev).cls.at(ncl).eta;
	/*
	significativit_histo->Fill(ev_seed, ev_sign_seed);
	clusterseed_histo->Fill(ev_adc_seed,ev_adc_p1+ev_adc_m1);
 	seed_m1_histo->Fill(ev_adc_seed,ev_adc_m1);
	seed_p1_histo->Fill(ev_adc_seed,ev_adc_p1);
	
	clusterseed_corr_histo->Fill(ev_adc_seed,clev.at(nev).cls.at(ncl).GetSides(3.));
	*/
int ladder=(ChanToLadder(ev_seed));
	if(ChanToSide(ev_seed)){ //n cases
	  if(eta>ETAMIN && eta<ETAMAX) ++eta_dist_n[ladder][(int)(ETASTEP*(eta)/ETARANGE)];
	 
	}else{

	  if(eta>ETAMIN && eta<ETAMAX) ++eta_dist_p[ladder][(int)(ETASTEP*(eta)/ETARANGE)];
	}

	//outputfile << nev <<"\t"<< ChanToLadder(ev_seed) <<"\t"<<ChanToSide(ev_seed)<<"\t"<<ev_seed%SIDE_CHAN+eta*CLUSTERCHANNELS/2<<std::endl;
	//WARNING! Plotting sign instead sigma!!!!
	double ev_sigma_seed=clev.at(nev).cls.at(ncl).sn[1];
	double ev_sigma_m1=clev.at(nev).cls.at(ncl).sn[0];
	double ev_sigma_p1=clev.at(nev).cls.at(ncl).sn[2];
	
	if(ev_sigma_seed<MAX_SIGMA2 && ev_sigma_m1<MAX_SIGMA2 && ev_sigma_p1<MAX_SIGMA2){
	  clustersigma_histo->Fill(ev_sigma_seed,sqrt(ev_sigma_m1*ev_sigma_m1+ev_sigma_p1*ev_sigma_p1));
	  clustersigma_m1_histo->Fill(ev_sigma_seed,ev_sigma_m1);
	  clustersigma_p1_histo->Fill(ev_sigma_seed,ev_sigma_p1);
	}
      }
    }
    
  for(int iev=0;iev<NCALIBEVENTS;iev++){
    for(int iva=0;iva<N_VA;iva++){
      common_noise_total->Fill(iva,CN_matrix_clean[iev][iva]);
    }
  }
  }
  double f_eta_p[N_LADDER][ETASTEP];
  double f_eta_n[N_LADDER][ETASTEP];
    for(int ladder=0;ladder<N_LADDER;++ladder){
      f_eta(ETASTEP,basement,eta_dist_n[ladder],f_eta_n[ladder]);
      f_eta(ETASTEP,basement,eta_dist_p[ladder],f_eta_p[ladder]);
      for (int ch=0;ch<ETASTEP;++ch){      
	f_eta_p_hist[ladder]->SetBinContent(ch+1,f_eta_p[ladder][ch]);
	f_eta_n_hist[ladder]->SetBinContent(ch+1,f_eta_n[ladder][ch]);
      }
    }
    //Now I want to correct the cluster using the f(eta). I stored all the vectors of clusters in a vector to avoid a loop on the events. 

    for(int vec=0;vec<(int)storage.size();++vec){
      for(int nev=0;nev<(int)(storage.at(vec).size());++nev){// I want to correct values for f_eta and fill all the histograms here.
	double x_p[N_LADDER]={0.};
	double x_n[N_LADDER]={0.};
	int seed_p[N_LADDER]={0};
	int seed_n[N_LADDER]={0};
	for(int ncl=0;ncl<(int)(storage.at(vec).at(nev).cls.size());++ncl){
	  //double countcluster[CLUSTERCHANNELS];
	  //int ev_seed=clev.at(nev).cls.at(ncl).seed;
	  int ev_seed=storage.at(vec).at(nev).cls.at(ncl).seed;
	  double ev_adc_seed=storage.at(vec).at(nev).cls.at(ncl).count[CLUSTERCHANNELS/2];
	  double ev_adc_p1=storage.at(vec).at(nev).cls.at(ncl).count[CLUSTERCHANNELS/2+1];
	  double ev_adc_m1=storage.at(vec).at(nev).cls.at(ncl).count[CLUSTERCHANNELS/2-1];
	  double ev_sign_seed=storage.at(vec).at(nev).cls.at(ncl).sn[CLUSTERCHANNELS/2];
	  double ev_sign_p1=storage.at(vec).at(nev).cls.at(ncl).sn[CLUSTERCHANNELS/2+1];
	  double ev_sign_m1=storage.at(vec).at(nev).cls.at(ncl).sn[CLUSTERCHANNELS/2-1];
	  //double eta=clev.at(nev).cls.at(ncl).eta;

	  double chargecenter=0.;
	  double totcharge=0.;
	  double eta=storage.at(vec).at(nev).cls.at(ncl).eta;

	  significativit_histo->Fill(ev_seed, ev_sign_seed);
	  clusterseed_histo->Fill(ev_adc_seed,ev_adc_p1+ev_adc_m1);
	  seed_m1_histo->Fill(ev_adc_seed,ev_adc_m1);
	  seed_p1_histo->Fill(ev_adc_seed,ev_adc_p1);
	  clusterseed_corr_histo->Fill(ev_adc_seed,storage.at(vec).at(nev).cls.at(ncl).GetSides(3.));
	  
	  int ladder=ChanToLadder(ev_seed);
	  int side=ChanToSide(ev_seed);
	  if(side){//n side
	    for(int cls=0;cls<CLUSTERCHANNELS;++cls){
	      chargecenter+=(storage.at(vec).at(nev).cls.at(ncl).count[cls]*(cls+f_eta_n[ladder][(int)eta/ETASTEP]*PITCH));
	      totcharge+=storage.at(vec).at(nev).cls.at(ncl).count[cls];
	    }
	    x_n[ladder]=(double)(ev_seed%SIDE_CHAN)+(chargecenter/totcharge-CLUSTERCHANNELS/2);
	    seed_n[ladder]=ev_seed%SIDE_CHAN;
	    outputfile << nev <<"\t"<< ladder <<"\t"<<side<<"\t"<<x_n[ladder] <<std::endl;


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
	eta_ADC_n_histo[ladder]->Fill(eta,storage.at(vec).at(nev).cls.at(ncl).GetCounts(3.));
	real_cluster_pos_n_hist[ladder]->Fill(ev_seed%SIDE_CHAN+eta*CLUSTERCHANNELS/2);
	for(int i=0;i<CLUSTERCHANNELS;++i)
	  cluster_shape_n_histo[ladder]->Fill(i,storage.at(vec).at(nev).cls.at(ncl).sn[i]);
	  }
	  else{//p side
	    for(int cls=0;cls<CLUSTERCHANNELS;++cls){
	      chargecenter+=(storage.at(vec).at(nev).cls.at(ncl).count[cls]*(cls+f_eta_p[ladder][(int)eta/ETASTEP]*PITCH));
	      totcharge+=storage.at(vec).at(nev).cls.at(ncl).count[cls];
	    }
	    
	    x_p[ladder]=(double)(ev_seed%SIDE_CHAN)+(chargecenter/totcharge-CLUSTERCHANNELS/2);
	    seed_p[ladder]=ev_seed%SIDE_CHAN;
	    outputfile << nev <<"\t"<< ladder <<"\t"<<side<<"\t"<<x_p[ladder] <<std::endl;

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
	eta_ADC_p_histo[ladder]->Fill(eta,storage.at(vec).at(nev).cls.at(ncl).GetCounts(3.));
	real_cluster_pos_p_hist[ladder]->Fill(ev_seed%SIDE_CHAN+eta*CLUSTERCHANNELS/2);
	for(int i=0;i<CLUSTERCHANNELS;++i)
	  cluster_shape_p_histo[ladder]->Fill(i,storage.at(vec).at(nev).cls.at(ncl).sn[i]);
	    
	  }
	  if(x_p[2]!=0. && x_p[3]!=0.){
	    deltax_central_ladder_p_hist->Fill(x_p[3]-x_p[2]);
	    deltax_central_ladder_noeta_p_hist->Fill(seed_p[3]-seed_p[2]);
	  }
	    	    
	  if(x_n[2]!=0. && x_n[3]!=0.){
	    deltax_central_ladder_n_hist->Fill(x_n[3]-x_n[2]);
	    deltax_central_ladder_noeta_n_hist->Fill(seed_n[3]-seed_n[2]);
	  }

	}
      }
    }

   
  mean1_histo->Write();
  sigma1_histo->Write();
  mean2_histo->Write();
  sigma2_histo->Write();
  sigma3_histo->Write();
  common_noise_total->Write();
  total_counts_adc->Write();
  total_counts_mean->Write();
  total_counts_clean->Write();
  good_chan_hist->Write();
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

  for(int ld=0;ld<N_LADDER;++ld){
    meanmean_hist[ld]->Write();
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
  //calibration_comp_mean->Write();
  //calibration_comp_sigma->Write();
  

  /*
  seed_frac_sum_histo->Write();
  seed_frac_m1_histo->Write();
  seed_frac_p1_histo->Write();
  sign_frac_sum_histo->Write();
  sign_frac_m1_histo->Write();
  sign_frac_p1_histo->Write();
  */
  std::stringstream Stream2;
  std::stringstream Stream2_open;
  std::stringstream Stream2_close;
  Stream2<<outputname <<"_results.pdf";
  Stream2_open<<outputname<<"_results.pdf[";
  Stream2_close<<outputname<<"_results.pdf]";


  TCanvas *out=new TCanvas();
  
  out->Print(Stream2_open.str().c_str());
  display_ladders2D(common_noise_total,"common_noise",";VA;ADC")->Print(Stream2.str().c_str());
  //display_ladders2D(sigma1_histo,"sigma_1",";chan;")->Print(Stream2.str().c_str());
  //display_ladders2D(sigma2_histo,"sigma_2",";chan;")->Print(Stream2.str().c_str());
  //display_ladders2D(sigma3_histo,"sigma_3",";chan;")->Print(Stream2.str().c_str());
  // display_ladders2D(total_counts_adc,"tot_counts",";chan;ADC")->Print(Stream2.str().c_str());
  //display_ladders2D(mean2_histo,"mean_2",";chan;ADC")->Print(Stream2.str().c_str());
  //display_ladders2D(total_counts_mean,"total_counts_mean",";chan;")->Print(Stream2.str().c_str());
  display_ladders2D(total_counts_clean,"total_counts_clean",";chan;")->Print(Stream2.str().c_str());
  display_ladders2D(significativit_histo,"sign",";chan;")->Print(Stream2.str().c_str());
  display_ladders1D(good_chan_hist,"good_channel",";chan;")->Print(Stream2.str().c_str());
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
  drawing6_1D(sum_adc_p_hist)->Print(Stream3.str().c_str());
  drawing6_1D(sum_adc_n_hist)->Print(Stream3.str().c_str());
  drawing6_1D(sum_adc_p_corr_hist)->Print(Stream3.str().c_str());
  drawing6_1D(sum_adc_n_corr_hist)->Print(Stream3.str().c_str());
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
  drawing1D(deltax_central_ladder_noeta_p_hist,0)->Print(Stream4.str().c_str());
  drawing1D(deltax_central_ladder_noeta_n_hist,0)->Print(Stream4.str().c_str());
  drawing1D(deltax_central_ladder_p_hist,0)->Print(Stream4.str().c_str());
  drawing1D(deltax_central_ladder_n_hist,0)->Print(Stream4.str().c_str());

  out->Print(Stream4_close.str().c_str());

  
  
}

void track_reconstruction(std::string namefile){
  gROOT->Reset();
  gDirectory->GetList()->Delete();

  //reading and storing calibration values
  std::ifstream infile (namefile.c_str(),std::ifstream::in);
  /*  need them later...
      std::vector ladder <int>;
      std::vector plan <int>;
      std::vector position <double>;
  */
  int ladder_index;
  int side_index;
  int ev_index;
  double real_position;

  while(!infile.eof()){
    infile >> ev_index >> ladder_index >> side_index >> real_position;


  }




}
    

void run1(){
  analysis("./Servo-OFF/Run-20161112/20161112-112421-Run_3C_37MeV_SERVO_EASIROC2.root","output_cal.txt","./run_22_12/30MeV");
  analysis("./Servo-OFF/Run-20161112/RUN_3C_51MeV_SERVO_EASIROC2_HOT.root","output_cal.txt","./run_22_12/50MeV");
  analysis("./Servo-OFF/Run-20161111/RUN_POS4XC_70MeV_EASIROC2_TrigMask0_SERVO_HOT.root","calib_70MeV.txt","./run_22_12/70MeV");
  analysis("./Servo-OFF/Run-20161112/RUN_3C_228MeV_SERVO_EASIROC2_HOT.root","output_cal.txt","./run_22_12/228MeV");
 analysis("./Servo-OFF/Run-20161112/RUN_3C_100MeV_SERVO_EASIROC2_HOT.root","output_cal.txt","./run_22_12/100MeV");
  analysis("./Servo-OFF/Run-20161112/RUN_3C_125MeV_SERVO_EASIROC2_HOT.root","output_cal.txt","./run_22_12/125MeV");
analysis("./Servo-OFF/Run-20161112/RUN_3C_154MeV_SERVO_EASIROC2_HOT.root","output_cal.txt","./run_22_12/154MeV");
  analysis("./Servo-OFF/Run-20161112/RUN_3C_174MeV_SERVO_EASIROC2_HOT.root","output_cal.txt","./run_22_12/174MeV");

}
