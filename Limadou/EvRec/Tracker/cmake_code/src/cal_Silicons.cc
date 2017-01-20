#include <iostream>
#include "TFile.h"
#include "TLatex.h"
#include "TH1.h"
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


const int n_chann = 4608; //number of channels 
const int n_steps = 1; // number of steps of 1000 events
int n_ev = 0; //number of events per step, default value. I assign it later in the code
const int VA_chan = 64; // number of channel in a VA
const int n_VA = 72; // total number of VA
const int ADC_CHAN = 4*VA_chan;
const int LADDER_CHAN = 3*ADC_CHAN;
const double sigmarange = 1.5;
const double MIN_NCHANSPERBIN_4CN=5;
const double MIN_SIGMA_CHANNEL_OFF=6;
const double MAX_SIGMA_CHANNEL_OFF=30;
const double MAX_NOISE_LEVEL=4; 
const double sigma_cut_gauss = 3.;

bool IsAlive(double sigma){
	bool deadoralive = 1;
	if (sigma<3.5) deadoralive = 0;
	return deadoralive;

}
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
void CommonNoise( int * channel, float *mean, float *sigma, float *max_sigma, float *noise, int * GEcounter){
	int va_index;
	//int *GEcounter = new int[n_VA];

	for(int i=0; i<n_VA; ++i) {
		GEcounter[i]=0;					
		noise[i]=0.;
	}
	for (int i = 0; i<n_chann; ++i){
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

double GetCleanedSigma(TH1F *h) {
    // Put the histo into something workable on
    int nb=h->GetNbinsX();

    double max=MIN_NCHANSPERBIN_4CN-1.; // minimum number of chans we want for CN calculation: careful! If split into more bins you cannot demand too much...
    double xmax=99999.;
    for(int ib=1; ib<nb+1; ++ib) {
        double x=h->GetBinCenter(ib);
        double y=h->GetBinContent(ib);
           
        if(x<=MIN_SIGMA_CHANNEL_OFF||x>=MAX_SIGMA_CHANNEL_OFF) continue; // dead channels or noisy channel. Go ahead.
        if(x>xmax+MAX_NOISE_LEVEL) break; // you have gone too far. The first 'structure' is the one you are interested in
        if(y>max) { // update max infos
            max=y;
            xmax=x;
        }
    }
    
    return xmax;
}


				//==================================================================================================
				//========================================== MAIN =================================================
				//==================================================================================================


void cal_Silicons(string name_root_file, string typeofparticle){

		//================= DECLARE AND INITIALIZE====================
	gDirectory->GetList()->Delete();
	TFile * data = new TFile( name_root_file.c_str(),"READ");
	TFile * data2 = new TFile("nongauss.root","RECREATE");
	TTree * tree = (TTree*)data->Get("T");// take the tree
	n_ev = tree->GetEntries();
	std::cout<<tree->GetEntries()<<std::endl;
	std::cout<<"Welcome to the silicon layers calibration!"<<std::endl;
	string name;
	string nameside;
	string title;

	TH1F * h_non_gauss = new TH1F("","",4608,0,4608);
	TH1F * significance = new TH1F("","",200,-10,10);

	int * ADC[n_chann]; // 4608 pointers to the ADC array !
	int * ADC2[n_ev]; // 1000 pointers to the ADC array !
	short strip[n_chann]; // pointer to the strip !
	float mean_1[n_chann]; // mean of the channels !
	float sigma_1[n_chann]; // sigma of the channels !
	float mean_2[n_chann]; // mean of the channels in +/- 3 sigma !
	float sigma_2[n_chann]; // sigma of the channels in +/- 3 sigma !
	float sigma_3[n_chann]; // sigma after comnoise subtraction !
	float * comnoise_real[n_ev]; // 1000 pointers to the comnoise array  !
	float * comnoise_real2[n_VA]; // 72 pointers to the comnoise array  !
	float * max[n_steps]; // max per step
	float mean_final[n_chann];
	float sigma_2_final[n_chann];
	float sigma_3_final[n_chann];
	float non_gauss[n_chann];
	float non_gauss_index[n_chann];

	for (int i=0; i<n_chann; ++i){
		mean_final[i]=0;
		sigma_2_final[i]=0.;
		sigma_3_final[i]=0.;
		non_gauss[i]=0.;
		non_gauss_index[i]=0.;
	}

	for (int s=0; s<n_steps; ++s){ max[s] = new float[n_VA];} // initialize the max
	for(int i=0; i<n_chann; ++i) { ADC[i] = new int[n_ev];} // initialize ADC one way
	for(int j=0; j<n_ev; j++){ADC2[j] = new int[n_chann];} // initialize ADC other way
	for (int j=0; j<n_ev; j++){ comnoise_real[j] = new float[n_VA];} // initialize comnoise_real one way
	for (int i=0; i<n_VA; i++){ comnoise_real2[i] = new float[n_ev];} // initialize comnoise_real other way
	for (int i=0; i<n_chann; i++){ 		//initialize all arrays
		strip[i]=0;
		mean_1[i]=0;
		mean_2[i]=0;
		sigma_1[i]=0;
		sigma_2[i]=0;
		sigma_3[i]=0;
	}

				//*****************************************************************
				//*****************************EXECUTE*****************************
				//*****************************************************************


	tree->SetBranchAddress("strip[4608]",&strip);// choose the branch 
	tree->SetBranchStatus("*",kFALSE);// close the other branch
	tree->SetBranchStatus("strip[4608]",kTRUE); // close the other branch

	std::cout<<"Calculation in process..."<<std::endl;

	for (int s=0; s<n_steps;++s){
		clock_t start = clock();
		for (int j=0; j<n_ev; j++){ 		// event index
			tree->GetEntry(s*n_ev+j);				// all channels for a fixed event
			for (int i=0; i<n_chann; i++){ 			// channel index
				mean_1[i]=mean_2[i]=0;
				sigma_1[i]=sigma_2[i]=sigma_3[i]=0;	
				ADC[i][j] = strip[i];	// filling the vector with channel/event ADCs
				ADC2[j][i] = strip[i];  // filling the vector with event/channel ADCs
				comnoise_real2[i/64][j]=0; //initialize to 0 common noise every 
				comnoise_real[j][i/64]=0; //initialize to 0 common noise every 
			}
		}

		//=========PEDESTAL_1==========
		for (int i=0; i<n_chann; i++){
			mean_1[i] = GetMean(ADC[i],0,5000); // all the allowed range
		}
		//=========SIGMA_1===========
		for (int i=0; i<n_chann; i++){ 
			sigma_1[i] = GetSigma(ADC[i],comnoise_real2[i/64],mean_1[i],0,5000);
		}

		//=========PEDESTAL_2==========
		for (int i=0; i<n_chann; i++){
			mean_2[i] = GetMean(ADC[i],mean_1[i]-3*sigma_1[i],mean_1[i]+3*sigma_1[i]);
		}
		//=========SIGMA_2===========
		for (int i=0; i<n_chann; i++){ 
			sigma_2[i] = GetSigma(ADC[i],comnoise_real2[i/64],mean_2[i],mean_2[i]-3*sigma_1[i],mean_2[i]+3*sigma_1[i]);
			
		}

		//==============SIGMA_2_VA==============	
		TH1F * h_sigma_VA;
		for (int t=0; t<n_VA; ++t){ //loop over all the VA
			name = "sigma_VA_"+to_string(t);
			name = title = "sigma_VA_"+to_string(t)+"_"+to_string(s);
			h_sigma_VA = new TH1F(name.c_str(),title.c_str(),200,0,100); // POSSIBLE BIAS
			for (int k=0; k<VA_chan; ++k){		//loop over the channel in a single VA
				h_sigma_VA->Fill(sigma_2[t*VA_chan+k]);
			}
			
			max[s][t] = GetCleanedSigma(h_sigma_VA);
		}
		//==================CALCULATION OF THE COMMON NOISE=====================
		int GEcounter[100][4608];
		for (int j=0; j<n_ev; j++){			// loop over the events
			CommonNoise(ADC2[j], mean_2, sigma_2, max[s], &comnoise_real[j][0], &GEcounter[s][0]);		
			for(int t=0; t<n_VA; ++t){
				comnoise_real2[t][j] = comnoise_real[j][t];			
			}

		}

		//==================CALCULATION OF THE SIGMA AFTER CN SUBTRACTION=====================
		for (int i=0; i<n_chann; i++){
			sigma_3[i] = GetSigma(ADC[i],comnoise_real2[i/64],mean_2[i],mean_2[i]-3*sigma_2[i],mean_2[i]+3*sigma_2[i]);
			
		}
		// time each step
	    clock_t stop = clock();
		std::cout<<"step [" <<int(s+1)<<"/"<<n_steps<< "] time = " <<((double)stop - (double)start)/CLOCKS_PER_SEC<<std::endl;

		//It sums all the sigmas and means each steps
		for(int i=0; i<n_chann; i++){
			mean_final[i] += mean_2[i];
			sigma_2_final[i] += sigma_2[i];
			sigma_3_final[i] += sigma_3[i];
		}


	}
	//It looks for non gaussian channels!!
	for (int i=0; i<n_chann; i++){
		non_gauss[i]=0;
		for(int j=0; j<n_ev; j++){
			double cont = (ADC[i][j]-mean_2[i]-comnoise_real2[int(i/VA_chan)][j])/sigma_3[i];
			// non gaussian index +1 if is out 3 sigma from the pedestal
			if(abs(cont)>sigma_cut_gauss){
				non_gauss[i]++;
			}
			//only p side of the first ladder
			if((i<768/2) && (IsAlive(sigma_3[i]))) significance->Fill(cont);
		}
		//out 3 sigma 0.0027 probability 
		non_gauss_index[i]=(non_gauss[i]-0.0027*n_ev)/(sqrt(0.0027*n_ev));
	}
	// nice output!!
	std::ofstream calib;
	string calib_name;
	calib_name = "calib_"+typeofparticle+".txt";
	calib.open(calib_name.c_str());
	calib <<"\t \t **************************************************************"<<endl
		  <<"\t \t ********* WELCOME TO SILICON LAYERS CALIBRATION **************"<<endl
		  <<"\t \t ********* run "<<name_root_file<<" ********************"<<endl
		  <<"\t \t ********* "<<typeofparticle<<" ************************"<<endl
		  <<"\t \t *************************************************************"<<std::endl<<std::endl;
	 
	for (int i=0; i<4608; i++){
		calib <<"channel_"<< i<<"_mean_2="<<mean_final[i]/n_steps<<" sigma_2="<<sigma_2_final[i]/n_steps<<" sigma_3="<<sigma_3_final[i]/n_steps<<" non_gauss_index="<<non_gauss_index[i]<<" "<<std::endl;
		//spectrum of non gaussian channels
		h_non_gauss->SetBinContent(i+1,non_gauss_index[i]);

	}
	calib.close();

	//***********************PLOTTING AND WRITING****************************

	TCanvas * c1 = new TCanvas();
	string name_non_gauss = "";
	name_non_gauss = "non_gauss_"+ typeofparticle +"_prova.pdf";
	h_non_gauss->Draw();
	c1->Print(name_non_gauss.c_str());

	TCanvas * c2 = new TCanvas();
	name_non_gauss = "significance_"+ typeofparticle +"_prova.pdf";
	significance->Draw();
	c2->SetLogy();
	c2->Print(name_non_gauss.c_str());
	data->Close();

	data2->cd();
	h_non_gauss->Write();
	significance->Write();
	data2->Close();
	

}

void Series(){
	cal_Silicons("FakeCalib-20161112-101327_SERVO_EASIROC2_HOT.root","protons_37MeV");
	cal_Silicons("FakeCalib-20161112-101327_SERVO_EASIROC2_HOT.root","protons_50MeV");
	cal_Silicons("20161111-213008-Run_48059_4660.root","protons_70MeV");
	cal_Silicons("FakeCalib-20161112-101327_SERVO_EASIROC2_HOT.root","protons_100MeV");
	cal_Silicons("FakeCalib-20161112-101327_SERVO_EASIROC2_HOT.root","protons_125MeV");
	cal_Silicons("FakeCalib-20161112-101327_SERVO_EASIROC2_HOT.root","protons_154MeV");
	cal_Silicons("FakeCalib-20161112-101327_SERVO_EASIROC2_HOT.root","protons_174MeV");
	cal_Silicons("FakeCalib-20161112-101327_SERVO_EASIROC2_HOT.root","protons_228MeV");
	cal_Silicons("FakeCalib-20161112-101327_SERVO_EASIROC2_HOT.root","protons_beamrotated.txt");
}



