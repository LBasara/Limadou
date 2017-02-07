#include "LEvRec0.hh"
#include "LTrackerTools.hh"
#include "LTrackerCluster.hh"
#include "LTrackerFittingTools.hh"
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

//GET THE MEAN OVER 1000 EVENTS
float GetMean( int * channel, float inf, float sup){
	float sum=0;
	int cont=0;
	for (int j=0; j<n_ev; j++){
		if ((channel[j]>=inf)&&(channel[j]<=sup)){			//mean_1 range = 0,5000, mean_2 range = +/- 3 sigma
			sum+=channel[j];
			cont++;
		}
	}
	return sum/cont;
}

//first attampt of correction function
double Correction_1(double eta, double par0, double par1, double par2){
	double result=-999;
	double fitting = par0+par1+par2;
	result = (par0/fitting)+(par1/fitting)*eta+(par2/fitting)*eta*eta;
	return 1./result; 
}
// final correction function, based on the parabolic fit parameters that are readed from "fitting_parameters... .txt"
double Correction_2(TProfile * histo, double x){
	double pos_x=-999;
	double Nbin = histo->GetNbinsX();
	double max = histo->GetBinContent(1);
	double min = histo->GetBinContent(Nbin+1);
	double first_corr = (max+min)/2.;
	histo->Scale(1/first_corr);
	pos_x = histo->GetBin(x);
	return 1/(histo->GetBinContent(pos_x));
}


//GET THE SIGMA OVER 1000 EVENTS
float GetSigma(int * channel, float * noise, float mean, float inf, float sup){
	float sum=0;
	int cont=0;
	for (int j=0; j<n_ev; j++){
			if ((channel[j]>=inf)&&(channel[j]<=sup)){				//sigma_1 range = 0,5000, sigma_2, sigma_3 range = +/- 3 sigma
			sum +=((channel[j]-mean-noise[j])*(channel[j]-mean-noise[j]));
			cont++;
		}
	}
	return sqrt(sum/cont);
}

//GET THE COMMONNOISE
void CommonNoise( int * channel, float *mean, float *sigma, float *max_sigma, float *noise, int * GEcounter/*, int * deadchannel*/){
	int va_index;
	for(int i=0; i<n_VA; ++i) {
		GEcounter[i]=0;					
		noise[i]=0.;
	}
	for (int i = 0; i<n_chann; ++i){
		// line for excluding non gaussian channels
		//if(deadchannel[i]!=0) continue;
		va_index = i/VA_chan;
		if ((channel[i]>=mean[i]-4*sigma[i])&&(channel[i]<=mean[i]+4*sigma[i])){ //selection on the channel +/- 3*sigma_2
			if ((sigma[i]>=(max_sigma[va_index]-sigmarange))&&(sigma[i]<=(max_sigma[va_index]+sigmarange))&&(sigma[i]>=6.)&&(sigma[i]<=30.)){ //selection on the sigma
				noise[va_index] += (channel[i]-mean[i]);
				++GEcounter[va_index];
			}
		}
	}
	for(int i=0; i<n_VA; ++i) {
		if(GEcounter[i]<1) 	noise[i]=0;   //NOT NEEDED. REDUNDANT
		else noise[i] /= GEcounter[i];					
	}
	return;
}
// it returns the sigma for the common noise calculation
double GetCleanedSigma(TH1F *h) {
    // Put the histo into something workable on
    int nb=h->GetNbinsX();

    double max=MIN_NCHANSPERBIN_4CN-1.; // minimum number of chans we want for CN calculation: careful! If split into more bins you cannot demand too much...
    double xmax=99999.;
    for(int ib=1; ib<nb+1; ++ib) {
        double x=h->GetBinCenter(ib);
        double y=h->GetBinContent(ib);
           
        if(x<=MIN_SIGMA_CHANNEL_OFF||x>=MAX_SIGMA_CHANNEL_OFF) continue; // dead channels or noisy channel. Go ahead.
        if(x>xmax+MAX_SIGMA_NOISE_LEVEL) break; // you have gone too far. The first 'structure' is the one you are interested in
        if(y>max) { // update max infos
            max=y;
            xmax=x;
        }
    }
    
    return xmax;
}
// it fits the counts vs. eta profile and it gets the parabola parameters
void Fitting_func(TProfile *histo[6], double par0[6], double par1[6], double par2[6]){
	for (int i=0; i< 6; i++){
		TF1 * myfit= new TF1("myfit","[0]+[1]*x+[2]*x*x", -0.1, 1.1);
		histo[i]->Fit("myfit");
		par0[i] = myfit->GetParameter("p0");
		par1[i] = myfit->GetParameter("p1");
		par2[i] = myfit->GetParameter("p2");
	}
	return;
}
// langaufunction taken from root examples
Double_t langaufun(Double_t *x, Double_t *par) {
      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;

      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}

Double_t langaufunmean(Double_t *x, Double_t *par) {
      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;

      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (x[0]*par[2] * step * sum * invsq2pi / par[3]);
}
// langau fit "ad hoc" for protons
void Langau(TH1F *histo[6], double langau_Mean[6], double langau_Sigma[6], double chi2red[6], double inf, double sup){
	for (int i=0; i< 6; i++){
		if ((i!=2)&&(i!=3)){
			langau_Mean[i]=0;
			langau_Sigma[i]=0;
			chi2red[i]=0;
		} else{
		float max_bin = histo[i]->GetMaximumBin();
		float max_center = histo[i]->GetBinCenter(max_bin);
		float integral = histo[i]->Integral(inf,sup);
		TF1 *langau = new TF1("langau",langaufun,inf,sup,4);
		Double_t values[4];
		values[0]=16; 
		values[1]=max_center; 
		values[2]=integral;
		values[3]=6.;
		langau->SetParameters(values);
		histo[i]->Fit("langau","","",inf,sup);
		values[0] = langau->GetParameter("p0");
		values[1] = langau->GetParameter("p1");
		values[2] = langau->GetParameter("p2");
		values[3] = langau->GetParameter("p3");
		TF1 *langaumean = new TF1("langaumean",langaufunmean,inf,sup,4);
		langaumean->SetParameters(values);
		double fitintegral = langau->Integral(inf,sup);
		double fitmean = langaumean->Integral(inf,sup);

		std::cout<<"------------------------>"
				 <<fitmean/fitintegral<<std::endl;
		//Test line
		langau_Mean[i] = fitmean/fitintegral;
		langau_Sigma[i]= GetMeanUncertanty(langau_Mean[i],langau->GetParameter(1),langau->GetParError(1),langau->GetParameter(0),langau->GetParError(0),langau->GetParameter(3),langau->GetParError(3),langau->GetParameter(2),langau->GetParError(2), inf, sup);//sqrt(pow(langau->GetParameter("p0"),2)+pow(langau->GetParameter("p3"),2))/sqrt(histo[i]->GetEntries());
		chi2red[i]=langau->GetChisquare();
		std::cout<<"**********ladder_"<<i<<" "<<langau_Mean[i]<<std::endl;
	}
	}
	return;
}

double GetMeanUncertanty(double mean, double MPV, double err_MPV, double sig_l, double err_sig_l, double sig_g, double err_sig_g, double area, double err_area, double inf, double sup){
	std::vector<double> resid;
	int dim1=10;
	int dim2=10;
	double MPV_L[dim1];
	double SIGMA_L[dim1];
	double SIGMA_G[dim2];
	double AREA_L[dim2];
	int cont = 0;
	while(cont<dim1){
		MPV_L[cont] = (MPV - err_MPV + 2*err_MPV/(dim1-1)*cont);
		SIGMA_L[cont] = (sig_l - err_sig_l + (2*err_sig_l/(dim1-1))*cont);
		++cont;
	}
	cont=0;
	while(cont<dim2){
		SIGMA_G[cont] = (sig_g - err_sig_g + (2*err_sig_g/(dim2-1))*cont);
		AREA_L[cont] = (area - err_area + (2*err_area/(dim2-1))*cont);
		++cont;
	}
	Double_t values[4];
	for(int i_mpv=0; i_mpv<dim1; ++i_mpv){
		values[1]=MPV_L[i_mpv];
		for(int i_sl=0; i_sl<dim1; ++i_sl){
			values[0]=SIGMA_L[i_sl];
			for(int i_sg=0; i_sg<dim2; ++i_sg){
				values[3]=SIGMA_G[i_sg];
				for(int i_area = 0; i_area<dim2; ++ i_area){
					values[2]=AREA_L[i_area];
					TF1 * mylangau = new TF1("mylangau",langaufun,inf,sup,4);
					mylangau->SetParameters(values);
					TF1 * mylangaumean = new TF1("mylangaumean",langaufunmean,inf,sup,4);
					mylangaumean->SetParameters(values);
					double myfitintegral = mylangau->Integral(inf,sup,1.e-3);
					double myfitmean = mylangaumean->Integral(inf,sup,1.e-3);
					double tmp = myfitmean/myfitintegral;
					double unc = tmp - mean;
					resid.push_back(unc);
					mylangaumean->Delete();
					mylangau->Delete();
				}
			}
		}
		std::cout<<float(i_mpv)/float(dim1)<< " " <<resid.back()<<std::endl;
	}
	double result = *std::max_element(resid.begin(),resid.end());
	std::cout<<"err "<<result<<std::endl;
	return result;
}


