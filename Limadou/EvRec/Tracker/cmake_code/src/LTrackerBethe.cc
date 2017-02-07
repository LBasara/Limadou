#include <iostream>
#include "string.h"
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include "TMath.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "TGraphErrors.h"

#include "LTrackerBethe.hh"



double Ftau(Double_t * x, double mass){
  double beta2 = (x[0]*x[0])/(1+x[0]*x[0]);
  double KinEn = mass*(-1+sqrt(x[0]*x[0]+1));
  double facttau = (KinEn*KinEn/8-(2*KinEn+1)*TMath::Log(2))/pow(KinEn+1,2);
  double result = 1-beta2+facttau;
  return result;
}

void Dummy(){
  Double_t eta=0.1;
  double *punt=&eta;
  while( eta<10){
    std::cout<<"eta = "<<eta <<" delta="<<Ftau(punt,0.511)<<std::endl;
    eta=eta+0.1;
  }
}

float Conv(double x, double mass){
  double toverm = x/mass;
  return sqrt(toverm*toverm+2*toverm);
}

Double_t bethefuncion(Double_t *x, Double_t *par){
  //x[0] = Conv(x[0],par[0]); // USE if kin. en. needed!!!!!!!!!!
  double K = 0.307075; //MeV*g^-1*cm^2 for A= 1 g/mol 
  double Z = 14; // atomic number of Silicon
  //double Z = 29; // atomic number of Copper
  double z2 = par[1]*par[1]; // charge squared of the incident particle
  double mass_e = 0.511; //MeV
  //double rho = 8.92;// copper g cm^-3
  double rho = 2.33 ;// silicon g cm^-3
  double A = 28; // mass number of Silicon
  //double A = 63.5; // mass number of Copper
  double r_e = 2.817*pow(10,-15);//radius electron in meters
  double moverM = mass_e/par[0];
  //double MeanExEn = 12.2 * pow(10,-6) *Z;// Mean exitation energy Silicon Z = 14 in eV
  double MeanExEn = 11.1 * pow(10,-6) *Z;// Mean exitation energy Silicon Z = 29 in eV
  double numT_max = 2*mass_e*x[0]*x[0];
  double denTmax = 1+2.*sqrt(x[0]*x[0]+1)*moverM+(moverM*moverM);
  double beta2inv =((1+x[0]*x[0])/(x[0]*x[0]));
  double fact1 = 0.5*TMath::Log(2*mass_e*x[0]*x[0]*numT_max/(MeanExEn*MeanExEn*denTmax));
  double fact2 = 1./beta2inv;
  double Enplasma = sqrt(rho*Z/A)*28.816*pow(10,-6);
  //double fact3 = 0.5*delta(x,-0.0254,3.2792,4.4190,0.1434,2.9044,0.08);//x_0 par_0,x_1 par_1,C par_2,a par_3,m par_4, delta_0 par_5
  double fact3 = 0.5*delta(x,0.2014,2.8715,4.4351,0.1492,3.2546,0.058);
  double factel=0;
  if(par[0]<1.){
    factel = Ftau(x,par[0]);
  }
  //double fact3 =0.;// TMath::Log(Enplasma/MeanExEn)+TMath::Log(x[0])-0.5;
  double result = K/A*Z*z2*beta2inv*(fact1+factel-fact2-fact3);
  return result;
}

void bethe(){
  double mass_proton=938.27; // mass of the incident particle
  double mass_muon=105.6;
  double mass_elec=0.511;
  double mass_helium = 3727.4;
  double charge_proton=1.;
  double charge_electron=1.;
  double charge_helium=2.;
  TCanvas * c1 = new TCanvas();
  //c1->SetLogx();
  c1->SetLogy();
  double minXaxis=30.;
  double maxXaxis=300.;
  TH1F *bg = (TH1F*)c1->DrawFrame(0.9*minXaxis,1.,maxXaxis*1.1,31.17);
  bg->GetXaxis()->SetTitle("#beta#gamma");
  bg->GetYaxis()->SetTitle("-#frac{dE}{dx} [MeV g^{-1} cm^{2}]");

  TF1 * bethefunc_e = new TF1("bethefunc_e",bethefuncion,minXaxis, maxXaxis, 2);
  bethefunc_e->SetParameter(0,mass_elec);
  bethefunc_e->SetParameter(1,charge_electron);
  bethefunc_e->SetLineColor(3);
  bethefunc_e->GetXaxis()->SetTitle("#gamma #beta");
  bethefunc_e->GetYaxis()->SetTitle("#frac{dE}{dx} [MeV/g cm]");
  bethefunc_e->SetMaximum(9.);
  bethefunc_e->SetMinimum(1.0);
  bethefunc_e->Draw("same");
  
  TF1 * bethefunc_p = new TF1("bethefunc_p",bethefuncion,minXaxis, maxXaxis, 2);
  bethefunc_p->SetParameter(0,mass_proton);
  bethefunc_p->SetParameter(1,charge_proton);
  bethefunc_p->SetLineColor(2);
  bethefunc_p->Draw("same");
  /*
  TF1 * bethefunc_mu = new TF1("bethefunc_mu",bethefuncion,3minXaxis, maxXaxis, 2);
  bethefunc_mu->SetParameter(0,mass_muon);
  bethefunc_mu->SetParameter(1,charge_electron);
  bethefunc_mu->SetLineColor(4);
  bethefunc_mu->Draw("same");
  */
  
  TF1 * bethefunc_helium = new TF1("bethefunc_helium",bethefuncion,minXaxis, maxXaxis, 2);
  bethefunc_helium->SetParameter(0,mass_helium);
  bethefunc_helium->SetParameter(1,charge_helium);
  bethefunc_helium->SetLineColor(4);
  bethefunc_helium->Draw("same");
  


  TLegend * leg = new TLegend(0.7,0.70,0.9,0.9);
  leg->SetHeader("Silicon","l");
  leg->AddEntry("bethefunc_e","e^{-}","l");
    leg->AddEntry("bethefunc_p","p","l");
    leg->AddEntry("bethefunc_helium","He","l");
    gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendTextSize(0.039);
    //leg->AddEntry("bethefunc_e","electrons","l");
    leg->Draw("same");
    c1->Print("Bethes.pdf");

  
}

double delta(Double_t * x, double x_0, double x_1, double Cbar, double a, double k, double delta_0 ){//x_0 par_0,x_1 par_1,C par_2,a par_3,k par_4, delta_0 par_5
  double result;
  double Logeta = TMath::Log10(x[0]);
  if(Logeta>=x_1){result = 2*TMath::Log(10)*Logeta-Cbar;}
  else if((Logeta>x_0)&&(Logeta<x_1)){result = 2*TMath::Log(10)*Logeta-Cbar+a*pow(x_1-Logeta,k);}
  else {result = delta_0*pow(10,2*(Logeta-x_0));}//for conductor, 0 for semiconduction
  return result;
}

Double_t Bethe_p(Double_t *x, Double_t *par){
  double MASS = 938.27;// MeV
  double K = 0.307075; //MeV*g^-1*cm^2 for A= 1 g/mol 
  double Z = 14; // atomic number of Silicon
  //double Z = 29; // atomic number of Copper
  double z2 = 1; // charge squared of the incident particle
  double mass_e = 0.511; //MeV
  //double rho = 8.92;// copper g cm^-3
  double rho = 2.33 ;// silicon g cm^-3
  double A = 28; // mass number of Silicon
  //double A = 63.5; // mass number of Copper
  double r_e = 2.817*pow(10,-15);//radius electron in meters
  double moverM = mass_e/MASS;
  //double MeanExEn = 12.2 * pow(10,-6) *Z;// Mean exitation energy Silicon Z = 14 in eV
  double MeanExEn = 11.1 * pow(10,-6) *Z;// Mean exitation energy Silicon Z = 29 in eV
  double numT_max = 2*mass_e*x[0]*x[0];
  double denTmax = 1+2.*sqrt(x[0]*x[0]+1)*moverM+(moverM*moverM);
  double beta2inv =((1+x[0]*x[0])/(x[0]*x[0]));
  double fact1 = 0.5*TMath::Log(2*mass_e*x[0]*x[0]*numT_max/(MeanExEn*MeanExEn*denTmax));
  double fact2 = 1./beta2inv;
  double Enplasma = sqrt(rho*Z/A)*28.816*pow(10,-6);
  //double fact3 = 0.5*delta(x,-0.0254,3.2792,4.4190,0.1434,2.9044,0.08);//x_0 par_0,x_1 par_1,C par_2,a par_3,m par_4, delta_0 par_5
  double fact3 = 0.5*delta(x,0.2014,2.8715,4.4351,0.1492,3.2546,0.058);
  //double fact3 =0.;// TMath::Log(Enplasma/MeanExEn)+TMath::Log(x[0])-0.5;
  double result = par[0]*K/A*Z*z2*beta2inv*(fact1-fact2-fact3);
  return result;
}

TCanvas * Plotter(TGraphErrors * graph){
  TCanvas * c1 = new TCanvas();
  double minbg=0.01;
  double maxbg=40.;
  //TH1 *bg=c1->DrawFrame(minbg,0.,maxbg,400.);
  graph->GetXaxis()->SetTitle("#beta#gamma");
  graph->GetYaxis()->SetTitle("Landau Mean");
  graph->SetMarkerStyle(7);
  // x is the kinetic energy of the protons, m_p=938 proton mass
  //TF1 * bbfunc = new TF1("bbfunc","[0]/x*TMath::Log([1]*(x*x+2*938*x))",0,300);
  TF1 * bbfunc = new TF1("bbfunc","Bethe_p",0,40,1);
  bbfunc->SetParameter(0,1000.);
  bbfunc->SetParameter(1,1.);
  //bbfunc->SetParameter(2,-2.1); //hv=16.9eV I=140eV
  gStyle->SetOptFit();
  c1->SetLogx();
  graph->SetMinimum(30.);
  graph->SetMaximum(400.);
  graph->Fit("bbfunc","R");
  graph->Draw("AP");

  return c1;
}


void langau_fit(){
  // the point at 70 MeV is not coherent with the others, and so is omitted
  //float energies[n_points]={37,51,100,125,154,174,228}; // MeV
  double mass_p = 938.; //MeV
  double mass_mu = 105.6; //MeV
  float energies[n_points]={Conv(37,mass_p),Conv(51,mass_p),Conv(100,mass_p),Conv(125,mass_p),Conv(154,mass_p),Conv(174,mass_p),Conv(228,mass_p)};//,Conv(4000,mass_mu)-float(5.)};//-5. is the correction muons->protons
  float zeros[n_points]={0};

  std::string title = "";
  std::ifstream langau_param("langau_param.txt");
  std::ofstream langau_param_william;
  langau_param_william.open("langau_param_william.txt");
  float MPV_read[6][2][n_points];
  float sigma_read[6][2][n_points];
  for(int ene=0; ene<n_points; ene++){
    langau_param_william <<"protons_"<<energies[ene]<< "MeV" <<std::endl;
    for (int o=0; o<6; o++){
      //reads parameter from "langau_param.txt"...pay attention deleting the 70 MeV lines!!
      langau_param >> MPV_read[o][0][ene] >> sigma_read[o][0][ene] >> MPV_read[o][1][ene] >> sigma_read[o][1][ene]; 
      // writes parameters in a more understanble txt ("langau_param_william.txt")
      langau_param_william << "ladder_"<< o <<" MPV_p =" << MPV_read[o][0][ene] << " sigma_p =" << sigma_read[o][0][ene] << " MPV_n =" << MPV_read[o][1][ene] <<  " sigma_n =" <<sigma_read[o][1][ene]<<std::endl;
    }

  }
/*  check paramters read by the lines above
  for(int ene=0; ene<n_points; ene++){
    for (int o=0; o<6; o++){
      std::cout<<MPV_read[o][0][ene] << " "<< sigma_read[o][0][ene]<< " " << MPV_read[o][1][ene] <<" "<< sigma_read[o][1][ene]<< std::endl;
    }
  }
*/
  //graphs for each side/ladder for all the energies
  TGraphErrors * global_langau[2][6];
  for(int i=0; i<6; i++){
    for(int j=0; j<2; j++){
      global_langau[j][i] = new TGraphErrors(n_points,energies,MPV_read[i][j],zeros,sigma_read[i][j]);

      if(j==0){
        title = "Ladder "+std::to_string(i)+" side p";
      }
      else{
        title = "Ladder "+std::to_string(i)+" side n";
      }   
      global_langau[j][i] -> SetTitle(title.c_str());
    }
  }

  
  TCanvas * plot= new TCanvas();
  plot->Print("langau_fit.pdf[");
  for(int i=0; i<6; i++){
    for(int j=0; j<2; j++){
      Plotter(global_langau[j][i])->Print("langau_fit.pdf");
    }
  }
  
}
