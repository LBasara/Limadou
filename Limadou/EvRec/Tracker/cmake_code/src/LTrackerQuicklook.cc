#include "detector_const.hh"
#include "analysis_const.hh"

#include "LEvRec0File.hh"
#include "LEvRec0.hh"
#include "LTrackerTools.hh"
#include "LTrackerCalibrationManager.hh"
#include "LTrackerCalibration.hh"
#include "silicon_calibration.hh"
#include "LTrackerQuicklook.hh"
#include "LTrackerMask.hh"

#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>

#include "TCanvas.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLine.h"



TCanvas *drawing6_chan_2D(TH2D *ladder[N_LADDER]){
  TCanvas *display=new TCanvas();
  TLine *line_adc[N_LADDER];
  TLine *line_va[N_LADDER][SIDE_VA];
  display->Divide(3,2);
  int pad;
  double ymin=ladder[0]->GetYaxis()->GetXmin();
  double ymax=ladder[0]->GetYaxis()->GetXmax();
  
  for(int i=0;i<N_LADDER;++i){
    if (i==0) pad=1;
    else if(i==1) pad=4;
    else if(i==2) pad=2;
    else if(i==3) pad=5;
    else if(i==4) pad=3;
    else if(i==5) pad=6;
    
    line_adc[i]=new TLine((double)ADC_CHAN,ymin,(double)ADC_CHAN,ymax);
    line_adc[i]->SetLineColor(1);
    line_adc[i]->SetLineStyle(1);
    
    display->cd(pad);
    gPad->cd(pad);
    gPad->SetLogz();
    ladder[i]->SetStats(0);
    ladder[i]->Draw("zcol");
    line_adc[i]->Draw("same");
    for(int iva=0;iva<SIDE_VA;++iva){
      line_va[i][iva]=new TLine(VA_CHAN+iva*VA_CHAN,ymin,VA_CHAN+iva*VA_CHAN,ymax);
      line_va[i][iva]->SetLineColor(20);
      line_va[i][iva]->SetLineStyle(2);

      line_va[i][iva]->Draw("same");
    }
  }
  return display;
}

TCanvas *drawing6_va_2D(TH2D *ladder[N_LADDER]){
  TCanvas *display=new TCanvas();
  TLine *line_adc[N_LADDER];
  TLine *line_va[N_LADDER][SIDE_VA];
  display->Divide(3,2);
  int pad;
  double ymin=ladder[0]->GetYaxis()->GetXmin();
  double ymax=ladder[0]->GetYaxis()->GetXmax();
  
  for(int i=0;i<N_LADDER;++i){
    if (i==0) pad=1;
    else if(i==1) pad=4;
    else if(i==2) pad=2;
    else if(i==3) pad=5;
    else if(i==4) pad=3;
    else if(i==5) pad=6;
    
    line_adc[i]=new TLine((double)SIDE_VA/2,ymin,(double)SIDE_VA/2,ymax);
    line_adc[i]->SetLineColor(1);
    line_adc[i]->SetLineStyle(1);
    
    display->cd(pad);
    gPad->cd(pad);
    gPad->SetLogz();
    ladder[i]->SetStats(0);
    ladder[i]->Draw("zcol");
    line_adc[i]->Draw("same");
    for(int iva=0;iva<SIDE_VA;++iva){
      line_va[i][iva]=new TLine(iva,ymin,iva,ymax);
      line_va[i][iva]->SetLineColor(20);
      line_va[i][iva]->SetLineStyle(2);

      line_va[i][iva]->Draw("same");
    }
  }
  return display;
}


TCanvas *drawing12_2D(TH2D *pside[N_LADDER],TH2D *nside[N_LADDER]){
   TCanvas *display=new TCanvas();
   display->Divide(3,4);
   int pad1;
   int pad2;
   for(int i=0;i<N_LADDER;++i){
     if(i==0){pad1=1; pad2=4;}
     else if(i==1){pad1=7; pad2=10;}
     else if(i==2){pad1=2; pad2=5;}
     else if(i==3){pad1=8; pad2=11;}
     else if(i==4){pad1=3; pad2=6;}
     else if(i==5){pad1=9; pad2=12;}
     display->cd(pad1);
     gPad->cd(pad1);
     gPad->SetLogz();
     pside[i]->SetStats(0);
     pside[i]->Draw("zcol");
     display->cd(pad2);
     gPad->cd(pad2);
     gPad->SetLogz();
     nside[i]->SetStats(0);
     nside[i]->Draw("zcol");
     
   }
   return display;
}

void QuickLook(std::string namefile){
  gROOT->Reset();
  gDirectory->GetList()->Delete();
  gStyle->SetPalette(1);
  
  TH2D *pedestal[N_SIDES][N_LADDER];
  TH2D *sigma[N_SIDES][N_LADDER];
  TH2D *NGindex[N_SIDES][N_LADDER];
  TH2D *CN[N_SIDES][N_LADDER];

  TH2D *counts_clean[N_SIDES][N_LADDER];
  TH2D *signal_noise[N_SIDES][N_LADDER];

  TH2D *temp_map[N_SIDES][N_LADDER];
  TH2D *gaus_map[N_SIDES][N_LADDER];
  
  for(int side=0;side<N_SIDES;++side){
    char plan=(side ? 'n' : 'p');
    for(int ld=0;ld<N_LADDER;++ld){
      pedestal[side][ld]=new TH2D(Form("pedestal_%c_%d",plan,ld),Form("pedestal_%c_%d;chan;ADC",plan,ld),SIDE_CHAN,0,SIDE_CHAN,100,1000,1500);
      sigma[side][ld]=new TH2D(Form("sigma_%c_%d",plan,ld),Form("sigma_%c_%d;chan;ADC",plan,ld),SIDE_CHAN,0,SIDE_CHAN,100,0,50);
      NGindex[side][ld]=new TH2D(Form("Non_gaussianty_%c_%d",plan,ld),Form("Non_gaussianty_%c_%d;chan;NG index",plan,ld),SIDE_CHAN,0,SIDE_CHAN,50,-10,10);
      CN[side][ld]=new TH2D(Form("Common_Noise_%c_%d",plan,ld),Form("Common_Noise_%c_%d;VA;ADC",plan,ld),SIDE_VA,0,SIDE_VA,200,-100,100);
      
      counts_clean[side][ld]=new TH2D(Form("Counts_clean_%c_%d",plan,ld),Form("Counts_clean_%c_%d;chan;ADC",plan,ld),SIDE_CHAN,0,SIDE_CHAN,200,-200,500);
     signal_noise[side][ld]=new TH2D(Form("signal_noise_%c_%d",plan,ld),Form("signal_noise_%c_%d;chan;S/N",plan,ld),SIDE_CHAN,0,SIDE_CHAN,100,-100,100);

     temp_map[side][ld]=new TH2D(Form("Hotness_map_%c_%d",plan,ld),Form("Hotness_map_%c_%d;channel;VA",plan,ld),VA_CHAN,0,VA_CHAN,SIDE_VA,0,SIDE_VA);
     gaus_map[side][ld]=new TH2D(Form("non_gaussianity_map_%c_%d",plan,ld),Form("non_gaussianity_map_%c_%d;channel;VA",plan,ld),VA_CHAN,0,VA_CHAN,SIDE_VA,0,SIDE_VA);
    }
  }
   LEvRec0File input(namefile.c_str());
  LEvRec0 ev;
  input.SetTheEventPointer(ev);
  const int MAXEVENTS = input.GetEntries();
  const int N_PKG=MAXEVENTS/NCALIBEVENTS_QL;

  //Calibration on file
  LTrackerCalibrationManager::GetInstance().LoadRun(namefile.c_str());
  LTrackerCalibration *cal=LTrackerCalibrationManager::GetInstance().Calibrate(NCALIBEVENTS_QL,0);
  //int slots=cal->GetNSlots();
  
  for(int ipk=0;ipk<N_PKG;++ipk){
    //std::cout<<"Processing events "<< ipk<<std::endl;
    double *sigma_chan=cal->GetSigma(ipk);
    double *mean_chan=cal->GetPedestal(ipk);
    double *NGindex_chan=cal->GetNGIndex(ipk);
    bool *CNmask=cal->GetCNMask(ipk);
    LTrackerMask hotchan_tmp=cal->GetMaskOnSigma(ipk,-999.,HOTCHANNELTHRESHOLD);
    LTrackerMask coldchan_tmp=cal->GetMaskOnSigma(ipk,COLDCHANNELTHRESHOLD,999.);
    LTrackerMask nongauschan_tmp=cal->GetMaskOnNGI(ipk,-999.,GAUSCHANNELTHRESHOLD);
    
    short data[NCHAN];
    double CN_va[N_VA];
    double counts_clean_chan[NCHAN];
    double sign_noise[NCHAN];
    
    
    for(int iev=0;iev<NCALIBEVENTS_QL;++iev){
      input.GetEntry(NCALIBEVENTS_QL*ipk+iev);
      for(int ichan=0;ichan<NCHAN;ichan++){
	data[ichan]=ev.strip[ichan];
      }
      ComputeCN(data,mean_chan,CNmask,CN_va);
      
      for(int ichan=0;ichan<NCHAN;++ichan){
	int chanva=ChanToVA(ichan);
	int chanside=ChanToSide(ichan);
	int chanladder=ChanToLadder(ichan);
	int chanonside=ChanToSideChan(ichan);
	int chanVAside=ChanToVA(ichan)%SIDE_VA;
	int chanonVA=ChanToVAChan(ichan);
  
	counts_clean_chan[ichan]=data[ichan]-CN_va[chanva]-mean_chan[ichan];
	sign_noise[ichan]=counts_clean_chan[ichan]/sigma_chan[ichan];
	counts_clean[chanside][chanladder]->Fill(chanonside,counts_clean_chan[ichan]);
	signal_noise[chanside][chanladder]->Fill(chanonside,sign_noise[ichan]);
      }
      for(int ld=0;ld<N_LADDER;++ld){
      for(int sd=0;sd<N_SIDES;++sd){
	for(int vaside=0;vaside<SIDE_VA;++vaside)
	  CN[sd][ld]->Fill(vaside,CN_va[vaside+sd*SIDE_VA+ld*LADDER_VA]);
      }
    }
    }
    bool *hot_bool=hotchan_tmp.GetBool();
    bool *cold_bool=coldchan_tmp.GetBool();
    bool *nongaus_bool=nongauschan_tmp.GetBool();
    /*
    for(int ichan=0;ichan<NCHAN;++ichan){
      int chanside=ChanToSide(ichan);
      int chanladder=ChanToLadder(ichan);
      int chanonside=ChanToSideChan(ichan);
      int chanVAside=ChanToVA(ichan)%SIDE_VA;
      int chanonVA=ChanToVAChan(ichan);
      if(!hot_bool[ichan]) temp_map[chanside][chanladder]->Fill(chanonVA,chanVAside,+10);
      //if(!cold_bool[ichan])temp_map[chanside][chanladder]->Fill(chanonVA,chanVAside,-10.);
      if(!nongaus_bool[ichan])gaus_map[chanside][chanladder]->Fill(chanonVA,chanVAside);
    }
    */
    

    for(int ichan=0;ichan<NCHAN;++ichan){
      int chanside=ChanToSide(ichan);
      int chanladder=ChanToLadder(ichan);
      int chanonside=ChanToSideChan(ichan);
      int chanVAside=ChanToVA(ichan)%SIDE_VA;
      int chanonVA=ChanToVAChan(ichan);
      pedestal[chanside][chanladder]->Fill(chanonside,mean_chan[ichan]);
      sigma[chanside][chanladder]->Fill(chanonside,sigma_chan[ichan]);
      NGindex[chanside][chanladder]->Fill(chanonside,NGindex_chan[ichan]);
      if(!hot_bool[ichan]) temp_map[chanside][chanladder]->Fill(chanonVA,chanVAside,+10);
      if(!cold_bool[ichan])temp_map[chanside][chanladder]->Fill(chanonVA,chanVAside,-10.);
      if(!nongaus_bool[ichan])gaus_map[chanside][chanladder]->Fill(chanonVA,chanVAside);
    }
  }

  /*
  std::stringstream OutName;
  std::stringstream OutName_open;
  std::stringstream OutName_close;
  OutName<<"QuickLook_"<<namefile<<".pdf";
  OutName_open<<"QuickLook_"<<namefile<<".pdf[";
  OutName_close<<"QuickLook_"<<namefile<<".pdf]";
  */
  TCanvas *output=new TCanvas();
  output->Print("QuickLook.pdf[","pdf");
 
  drawing6_chan_2D(pedestal[0])->Print("QuickLook.pdf","pdf");
    drawing6_chan_2D(pedestal[1])->Print("QuickLook.pdf","pdf");
    drawing6_chan_2D(sigma[0])->Print("QuickLook.pdf","pdf");
    drawing6_chan_2D(sigma[1])->Print("QuickLook.pdf","pdf");
    drawing6_va_2D(CN[0])->Print("QuickLook.pdf","pdf");
    drawing6_va_2D(CN[1])->Print("QuickLook.pdf","pdf");
    drawing6_chan_2D(NGindex[0])->Print("QuickLook.pdf","pdf");
    drawing6_chan_2D(NGindex[1])->Print("QuickLook.pdf");
    drawing6_chan_2D(counts_clean[0])->Print("QuickLook.pdf","pdf");
    drawing6_chan_2D(counts_clean[1])->Print("QuickLook.pdf","pdf");
    drawing6_chan_2D(signal_noise[0])->Print("QuickLook.pdf","pdf");
    drawing6_chan_2D(signal_noise[1])->Print("QuickLook.pdf","pdf");
    drawing12_2D(temp_map[0],temp_map[1])->Print("QuickLook.pdf","pdf");
    drawing12_2D(gaus_map[0],gaus_map[1])->Print("QuickLook.pdf","pdf");
    /*
    drawing6_2D(gaus_map[0])->Print("QuickLook.pdf");
    drawing6_2D(gaus_map[1])->Print("QuickLook.pdf");
    drawing6_2D(temp_map[0])->Print("QuickLook.pdf");
    drawing6_2D(temp_map[1])->Print("QuickLook.pdf");
    */
  output->Print("QuickLook.pdf]","pdf");

}
