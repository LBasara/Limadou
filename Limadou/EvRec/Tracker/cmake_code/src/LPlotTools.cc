#include "LEvRec0.hh"
#include "LTrackerTools.hh"
#include "LTrackerCluster.hh"
#include "analysis_follega.hh"

#include <iostream>
#include "TFile.h"
#include "TLatex.h"
#include "TH1.h"
#include "TMath.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH2.h"
#include "TPad.h"
#include "string.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLine.h"
#include <string>
#include "TPaveLabel.h"
#include <time.h>
#include <vector>
#include <stdlib.h>
#include <algorithm>
#include <fstream>
#include <string>
#include "string.h"

TCanvas * Plot_2d(TH2F * histo, string xaxis, string yaxis){
  TCanvas * c2 = new TCanvas();
  histo->GetXaxis()->SetTitle(xaxis.c_str());
  histo->GetYaxis()->SetTitle(yaxis.c_str());
  gStyle->SetPalette(1);
  gPad->SetLogz();
  gPad -> SetRightMargin(0.15);
  histo->Draw("COLZ");
  return c2;

}

TCanvas * Plot_1d(TH1F * histo, string xaxis, string yaxis){
  TCanvas * c2 = new TCanvas();
  histo->GetXaxis()->SetTitle(xaxis.c_str());
  histo->GetYaxis()->SetTitle(yaxis.c_str());
  histo->Draw();
  return c2;

}

TCanvas * Plot6_2d(TH2F * histo, string name, string xaxis, string yaxis, float min, float max, string cond, string cond2){
    TCanvas * c2 = new TCanvas();
    c2->SetTitle(name.c_str());
    c2->Divide(3,2);
    histo -> GetXaxis()->SetTitle(xaxis.c_str());
    histo -> GetYaxis()->SetTitle(yaxis.c_str());
    int binx = histo->GetNbinsX();
  int biny = histo->GetNbinsY();
  TH2F * histo_six2[6];
  TLine *line[6];
  TLine *line2[6];
  string title, name_0;
  name_0 = name;
  float pos[6] = {0,3,1,4,2,5};
  for ( int w=0; w<6; w++){
    name = name_0;
    name = name + "_ladder_"+to_string(w);
    title = name;
    histo_six2[w] = new TH2F(name.c_str(),title.c_str(),binx/6,0,binx/6,biny,min,max);
    line[w] = new TLine(binx/12,min,binx/12,max);
    line2[w] = new TLine(0,5.,4608,5.);
    for (int ix = 0; ix<binx/6; ix++){
      for (int iy = 0; iy<biny; iy++){
      histo_six2[w]->SetBinContent(ix+1, iy+1, histo->GetBinContent(binx/6*w+ix+1,iy+1));
      }
    }
    histo_six2[w] -> GetXaxis()->SetTitle(xaxis.c_str());
    histo_six2[w] -> GetYaxis()->SetTitle(yaxis.c_str());
    c2 -> cd(pos[w]+1);
    gPad -> cd(pos[w]+1);
      gPad -> SetRightMargin(0.15);
      gStyle -> SetOptStat(0);
      gPad->SetLogz();
      gStyle->SetPalette(1);
      histo_six2[w]->Draw("COLZ");
      if (cond=="yes"){ 
        line[w]->SetLineColor(2);
        line[w]->SetLineStyle(2);
        line[w]->Draw("same");
    }
    if (cond2=="5sign"){
      line2[w]->SetLineColor(2);
        line2[w]->SetLineStyle(2);
        line2[w]->Draw("same");
    }
  }
    return c2;
}
  
TCanvas * Plot6_1d(TH1F * histo, string name, string xaxis, string yaxis ,float min, float max, string cond){
  TCanvas * c20 = new TCanvas();
  c20->SetTitle(name.c_str());
  c20->Divide(3,2);
  int bin = histo->GetNbinsX();
  TH1F * histo_six[6];
  TLine *line[6];
  string title, name_0;
  name_0 = name;
  float pos[6] = {0,3,1,4,2,5};
  for ( int w=0; w<6; w++){
    name = name_0;
    name = name + "ladder_"+to_string(w);
    title = name;
    histo_six[w] = new TH1F(name.c_str(),title.c_str(),bin/6,0,bin/6);
    line[w] = new TLine(bin/12,min,bin/12,max);
    for (int i = 0; i<bin/6; i++){
      histo_six[w]->SetBinContent(i+1,histo->GetBinContent(bin/6*w+i+1));
    }
    histo_six[w] -> GetXaxis()->SetTitle(xaxis.c_str());
    histo_six[w] -> GetYaxis()->SetTitle(yaxis.c_str());
    histo_six[w]->SetMaximum(max);
      histo_six[w]->SetMinimum(min);
    c20 -> cd(pos[w]+1);
      gStyle -> SetOptStat(0);
      histo_six[w]->Draw();
      if (cond=="yes"){ 
        line[w]->SetLineColor(2);
        line[w]->SetLineStyle(2);
        line[w]->Draw("same");
    }
  }

  return c20;
}

TCanvas * Plot_6histo_1d(TH1F * histo[6], string xaxis, string yaxis, string log){
  TCanvas * c6 = new TCanvas();
  c6->Divide(3,2);
  float pos[6] = {0,3,1,4,2,5};
  for ( int w=0; w<6; w++){
    c6 -> cd(pos[w]+1);
    histo[w]->GetXaxis()->SetTitle(xaxis.c_str());
    histo[w]->GetYaxis()->SetTitle(yaxis.c_str());
    //gStyle -> SetOptStat(0);
    if(log=="log") gPad->SetLogy();
    histo[w]->Draw();
  }
  return c6;
}
TCanvas * Plot_6histo_profile(TProfile * histo[6], string xaxis, string yaxis, string log){
  TCanvas * c6 = new TCanvas();
  c6->Divide(3,2);
  float pos[6] = {0,3,1,4,2,5};
  for ( int w=0; w<6; w++){
    c6 -> cd(pos[w]+1);
    histo[w]->GetXaxis()->SetTitle(xaxis.c_str());
    histo[w]->GetYaxis()->SetTitle(yaxis.c_str());
    //gStyle -> SetOptStat(0);
    if(log=="log") gPad->SetLogy();
    histo[w]->Draw();
  }
  return c6;
}

TCanvas * Plot_6histo_2d(TH2F * histo[6], string xaxis, string yaxis){
  TCanvas * c6 = new TCanvas();
  c6->Divide(3,2);
  float pos[6] = {0,3,1,4,2,5};
  for ( int w=0; w<6; w++){
    c6 -> cd(pos[w]+1);
    histo[w]->GetXaxis()->SetTitle(xaxis.c_str());
    histo[w]->GetYaxis()->SetTitle(yaxis.c_str());
    gStyle -> SetOptStat(0);
    gPad->SetLogz();
    histo[w]->Draw("COLZ");
  }
  return c6;
}