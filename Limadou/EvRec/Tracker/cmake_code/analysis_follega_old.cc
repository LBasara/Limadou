#include "LEvRec0.hh"
#include "LTrackerTools.hh"
#include "LTrackerCluster.hh"
#include "analysis_follega.hh"
#include "LPlotTools.hh"
#include "LTrackerFittingTools.hh"

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

/*
const int n_chann = 4608; //number of channels 
int n_steps = 0; // number of steps of 1000 events
const int n_ev = 1000; //number of events per step
const int VA_chan = 64; // number of channel in a VA
const int n_VA = 72; // total number of VA
const int ADC_CHAN = 3*VA_chan;
const int LADDER_CHAN = 4*ADC_CHAN;
const double sigmarange = 1.5;
const double MIN_NCHANSPERBIN_4CN=5;
const double MIN_SIGMA_CHANNEL_OFF=6;
const double MAX_SIGMA_CHANNEL_OFF=30;
const double MAX_SIGMA_NOISE_LEVEL=4; 
const double sigma_cut=6.5;
const double gaussian_threshold = 3.;
*/


struct event {
  int entry;
  std::vector<LTrackerCluster> cls;
};


				//==================================================================================================
				//========================================== MAIN =================================================
				//==================================================================================================


void analysis_follega(string root_data_file, string typeofparticle, string file_calib, string file_fit_par){

		//================= DECLARE AND INITIALIZE====================
	//command for deleting all the objects, usfull if you want to write 
	gDirectory->GetList()->Delete();
	std::cout<<"Welcome to the silicon layers calibration!"<<std::endl;

	//********************* Reader ********************
	// it reads from the file analysis ... _2.root, do not delete it

	string name;
	string nameside;
	string title;
	string file_name = "analysis_"+typeofparticle+"_2.root";
	
	TFile * getter = new TFile(file_name.c_str());
	TProfile * getmyprofile[2][6];
	for(int side = 0; side < 2; side++){
		for(int w = 0; w <6; w++){
			nameside = "p";
			if (side==1) nameside = "n";
			name = "h_myprofile_"+to_string(w)+"_side_"+nameside;
			getmyprofile[side][w] = (TProfile*)getter->Get(name.c_str());
		}
	}

	TFile * getter2 = new TFile(file_name.c_str());
	TProfile * getmyprofile2[2][6];
	for(int side = 0; side < 2; side++){
		for(int w = 0; w <6; w++){
			nameside = "p";
			if (side==1) nameside = "n";
			name = "h_sumcontcluster_eta_prof"+to_string(w)+"_side_"+nameside;
			getmyprofile2[side][w] = (TProfile*)getter2->Get(name.c_str());
		}
	}
	
	//********************* End Reader ********************

	TFile * data = new TFile(root_data_file.c_str(),"READ");
	TTree * tree = (TTree*)data->Get("T");// take the tree
	//n_steps = int(tree->GetEntries()/1000);

	TH1F * h_mean_1 = new TH1F("h_mean_1","Pedestal",4608,0,4608); 
	TH1F * h_sigma_1 = new TH1F("h_sigma_1","Sigma_raw",4608,0,4608);   
	TH2F * h_comnoise = new TH2F("h_comnoise","Common Noise",72,0,72,500,-100,100);
	TH1F * h_sumcontcluster[2][6];
	TH1F * h_sumcontcluster_uncorrected[2][6];
	TH1F * h_sumcontcluster_corr2[2][6];
	TH1F * h_sumcontcluster_corr3[2][6];
	TH1F * h_sumsigncluster[2][6];
	TH2F * h_ratiocontcluster_plus[2][6];
	TH2F * h_ratiocontcluster_minus[2][6];
	TH2F * h_ratiocontcluster[2][6];
	TH2F * h_ratiosigncluster_plus[2][6];
	TH2F * h_ratiosigncluster_minus[2][6];
	TH2F * h_ratiosigncluster[2][6];
	TProfile * h_sumcontcluster_eta_prof[2][6];
	TH1F * h_countseed[2][6];
	TH1F * h_signseed[2][6];
	TH1F * h_clustersize[2][6];
	TH1F * h_cls_barycenter[2][6];
	TH1F * h_cls_barycenter_abs[2][6];
	TH1F * h_eta_distribution[2][6];
	TH2F * h_sumcontcluster_eta[2][6];	
	TProfile * h_myprofile[2][6];
	TH2F * h_correlation_adc[2][6];
	TH2F * h_corr_adc1chan[24];
	//checking correlation adc
	for (int w=0; w<24; w++){
		name = "h_corr_adc1chan_"+to_string(w);
		h_corr_adc1chan[w] = new TH2F(name.c_str(),name.c_str(),1000,-150,300,1000,-150,300);	
	}		


	

	// customized bins of a TProfile, for an experiment on binning
	Double_t edges_myprofile[11] = {-1.,-0.9,-0.8,-0.3,-0.15,0.,0.15,0.3,0.8,0.9,1.};


	// declaration of all the histograms, side(2)-ladder(6)
	for (int side = 0; side < 2; side ++){
		nameside = "p";
		if (side==1) nameside = "n";
		for (int w=0; w<6; w++){
			name = "h_sumcontcluster_"+to_string(w)+"_side_"+nameside;
			h_sumcontcluster[side][w] = new TH1F(name.c_str(),name.c_str(),520,-20,600);
			name = "h_sumcontcluster_uncorrected_"+to_string(w)+"_side_"+nameside;
			h_sumcontcluster_uncorrected[side][w] = new TH1F(name.c_str(),name.c_str(),520,-20,600);
			name = "h_sumcontcluster_corr2_"+to_string(w)+"_side_"+nameside;
			h_sumcontcluster_corr2[side][w] = new TH1F(name.c_str(),name.c_str(),520,-20,600);
			name = "h_sumcontcluster_corr3_"+to_string(w)+"_side_"+nameside;
			h_sumcontcluster_corr3[side][w] = new TH1F(name.c_str(),name.c_str(),520,-20,600);
			name = "h_sumcontcluster_eta"+to_string(w)+"_side_"+nameside;
			h_sumcontcluster_eta[side][w] = new TH2F(name.c_str(),name.c_str(),30,-1,1,600,0,300);
			name = "h_sumcontcluster_eta_prof"+to_string(w)+"_side_"+nameside;
			h_sumcontcluster_eta_prof[side][w] = new TProfile(name.c_str(),name.c_str(),30,-1,1);
			name = "h_sumsigncluster_"+to_string(w)+"_side_"+nameside;
			h_sumsigncluster[side][w] = new TH1F(name.c_str(),name.c_str(),300,0,300);
			name = "h_ratiocontcluster_plus_"+to_string(w)+"_side_"+nameside;
			h_ratiocontcluster_plus[side][w] = new TH2F(name.c_str(),name.c_str(),400,0,100,200,-20,20);
			name = "h_ratiocontcluster_minus_"+to_string(w)+"_side_"+nameside;
			h_ratiocontcluster_minus[side][w] = new TH2F(name.c_str(),name.c_str(),400,0,100,200,-20,20);
			name = "h_ratiocontcluster_"+to_string(w)+"_side_"+nameside;
			h_ratiocontcluster[side][w] = new TH2F(name.c_str(),name.c_str(),400,0,100,200,-20,20);
			name = "h_ratiosigncluster_plus_"+to_string(w)+"_side_"+nameside;
			h_ratiosigncluster_plus[side][w] = new TH2F(name.c_str(),name.c_str(),400,0,100,200,-20,20);
			name = "h_ratiosigncluster_minus_"+to_string(w)+"_side_"+nameside;
			h_ratiosigncluster_minus[side][w] = new TH2F(name.c_str(),name.c_str(),400,0,100,200,-20,20);
			name = "h_ratiosigncluster_"+to_string(w)+"_side_"+nameside;
			h_ratiosigncluster[side][w] = new TH2F(name.c_str(),name.c_str(),400,0,100,200,-20,20);
			name = "h_clustersize_"+to_string(w)+"_side_"+nameside;
			h_clustersize[side][w] = new TH1F(name.c_str(),name.c_str(),6,0,6);
			name = "h_cls_barycenter_"+to_string(w)+"_side_"+nameside;
			h_cls_barycenter[side][w] = new TH1F(name.c_str(),name.c_str(),40,0,5);
			name = "h_cls_barycenter_abs_"+to_string(w)+"_side_"+nameside;
			h_cls_barycenter_abs[side][w] = new TH1F(name.c_str(),name.c_str(),40,0,5);
			name = "eta_distribution_"+to_string(w)+"_side_"+nameside;
			h_eta_distribution[side][w] = new TH1F(name.c_str(),name.c_str(),30,-1,1);
			name = "h_myprofile_"+to_string(w)+"_side_"+nameside;
			h_myprofile[side][w] = new TProfile(name.c_str(),name.c_str(),10,edges_myprofile);
			name = "h_correlation_adc_"+to_string(w)+"_side_"+nameside;
			h_correlation_adc[side][w] = new TH2F(name.c_str(),name.c_str(),1000,-150,300,1000,-150,300);
		}
	}
	
	for (int side = 0; side < 2; side ++){
		nameside = "p";
		if (side==1) nameside = "n";
		for (int w=0; w<6; w++){
			name = "h_countseed_"+to_string(w)+"_side_"+nameside;
			h_countseed[side][w] = new TH1F(name.c_str(),name.c_str(),520,-20,500);
			name = "h_signseed_"+to_string(w)+"_side_"+nameside;
			h_signseed[side][w] = new TH1F(name.c_str(),name.c_str(),520,-20,500);
		}
	}
	// declaration of all the other histograms
	TH2F * h_count = new TH2F("h_count","ADC - Pedestal - CommonNoise",4608,0,4608,500,-500,500); 
	TH2F * h_sign = new TH2F("h_sign","Significance",4608,0,4608,200,-50,50); 
	TH1F * h_good_events = new TH1F("h_good_events","Events above 5 #sigma",4608,0,4608);
	TH1F * h_mean_2 = new TH1F("h_mean_2","Pedestal only within +/- sigma",4608,0,4608);
	TH1F * h_sigma_2 = new TH1F("h_sigma_2","Sigma_raw only within +/- 3sigma",4608,0,4608);
	TH1F * h_sigma_3 = new TH1F("h_sigma_3","Sigma_raw only within +/- 3sigma, considering common noise",4608,0,4608);
	TH2F * h_sigma_tot = new TH2F("h_sigma_tot","Sigma_raw only within +/- 3sigma, considering common noise, all events",4608,0,4608,400,0,40);
	TH2F * h_sigma_tot_prime = new TH2F("h_sigma_tot_prime","Sigma_raw all events",4608,0,4608,400,0,200);
	TH2F * h_sigma_tot_second = new TH2F("h_sigma_tot_second","Sigma_raw only within +/- 3sigma all events",4608,0,4608,400,0,40);
	TH2F * h_ADC_tot = new TH2F("h_ADC_tot","ADC all events",4608,0,4608,500,0,3000);
	TH2F * h_mean_tot = new TH2F("h_mean_tot","mean all events",4608,0,4608,500,0,3000);
	TH2F * h_ADC_mean_tot = new TH2F("h_ADC_mean_tot","ADC-mean all events",4608,0,4608,500,-500,500);
	TH2F * h_correlation = new TH2F("h_correlation","h_correlation", 1000, 0, 1000, 700, -200, 500);
	TH2F * h_correlation_plus_minus = new TH2F("h_correlation_plus_minus","h_correlation_plus_minus", 1000, 0, 1000, 700, -200, 500);
	TH2F * h_correlation_minus = new TH2F("h_correlation_minus","correlation_plot_minus, seed-1", 1000, 0, 1000, 700, -200, 500);
	TH2F * h_correlation_plus = new TH2F("h_correlation_plus","correlation_plot_plus, seed+1", 1000, 0, 1000, 700, -200, 500);
	TH2F * h_correlation_sign = new TH2F("h_correlation_sign","h_correlation_sign", 100, 0, 100, 100, -20, 100);
	TH2F * h_correlation_sign_plus_minus = new TH2F("h_correlation_sign_plus_minus","h_correlation_sign_plus_minus", 100, 0, 100, 100, -20, 100);
	TH2F * h_correlation_sign_minus = new TH2F("h_correlation_sign_minus","correlation_sign_plot, seed-1", 100, 0, 100, 100, -20, 100);
	TH2F * h_correlation_sign_plus = new TH2F("h_correlation_sign_plus","correlation_sign_plot, seed+1", 100, 0, 100, 100, -20, 100);
	TH1F * h_distribution_seed_5 = new TH1F("h_distribution_seed_5","h_distribution_seed_5 #sigma",100, 0, 100);
	TH1F * h_distribution_seed_4 = new TH1F("h_distribution_seed_4","h_distribution_seed_4 #sigma",100, 0, 100);
	TH1F * h_distribution_seed_35 = new TH1F("h_distribution_seed_35","h_distribution_seed_3.5 #sigma",100, 0, 100);
	TH1F * h_distribution_seed_3 = new TH1F("h_distribution_seed_3","h_distribution_seed_3 #sigma",100, 0, 100);
	TH1F * sigma_mean = new TH1F("sigma_mean","sigma_mean",4608,0,4608);
	TH1F * mean_mean = new TH1F("mean_mean","mean_mean",4608,0,4608);
	TH1F * h_non_gauss = new TH1F("h_non_gauss","h_non_gauss",4608,0,4608);
	TH2F * h_checkdeadchannel = new TH2F("h_checkdeadchannel","h_checkdeadchannel",4608,0,4608,400,0,200);
	TH2F * checkRobCode = new TH2F("checkRobCode","h_correlation",1000, 0, 1000, 700, -200, 500);

	//declaration of all the histograms on the ladder
	TH1F * checkped[6];
	string checkpedname;
	for(int ilad = 0; ilad<6; ++ilad){
		checkpedname = "checkpedname_ladder_"+to_string(ilad);
		checkped[ilad] = new TH1F(checkpedname.c_str(),checkpedname.c_str(),110,0,110);
	}

	TFile * analysis2 = new TFile(file_name.c_str(),"RECREATE");
	file_name = "analysis_"+typeofparticle+"_3.root";
	TFile * analysis3 = new TFile(file_name.c_str(),"RECREATE");
	for(int t=0; t<n_VA; ++t){
			name = "sigma_VA_"+to_string(t);
			analysis3->mkdir(name.c_str());
	}

	int * counter;
	counter = new int[6];
	float content_plus,content_minus,content;
	float significance_plus,significance_minus,significance,signif;
	int * ADC[n_chann]; // 4608 pointers to the ADC array !
	int * ADC2[n_ev]; // 1000 pointers to the ADC array !
	short strip[n_chann]; // pointer to the strip !
	UInt_t event_index; // event index in the Tree
	float mean_1[n_chann]; // mean of the channels !
	float sigma_1[n_chann]; // sigma of the channels !
	float mean_2[n_chann]; // mean of the channels in +/- 3 sigma !
	float sigma_2[n_chann]; // sigma of the channels in +/- 3 sigma !
	double sigma_3[n_chann]; // sigma after comnoise subtraction !
	float * comnoise_real[n_ev]; // 1000 pointers to the comnoise array  !
	float * comnoise_real2[n_VA]; // 72 pointers to the comnoise array  !
	float * max[n_steps]; // max per step
	double * clearchann[n_ev];
	double * clearchann2[n_ev];
	float mean_final[n_chann];
	float sigma_final[n_chann];
	float non_gauss[n_chann];
	float clearchanntot[n_chann];
	float clearchanntot2[n_chann];
	int MASK[n_chann];
	for(int i=0; i<n_chann; i++){
		if(int(i*6/n_chann)==0) MASK[i]=3;
		if(int(i*6/n_chann)==1) MASK[i]=5;
		if(int(i*6/n_chann)==2) MASK[i]=3;
		if(int(i*6/n_chann)==3) MASK[i]=5;
		if(int(i*6/n_chann)==4) MASK[i]=3;
		if(int(i*6/n_chann)==5) MASK[i]=10;
	}

	for(int i=0; i<n_ev; ++i) { clearchann[i] = new double[n_chann];} 
	for(int i=0; i<n_ev; ++i) { clearchann2[i] = new double[n_chann];} 
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
		clearchanntot[i]=0;
		clearchanntot2[i]=0;
	}
				//*****************************************************************
				//*****************************EXECUTE*****************************
				//*****************************************************************

	std::cout<<"Calculation in process..."<<typeofparticle<<std::endl;

	tree->SetBranchAddress("strip[4608]",&strip);// choose the branch 
	tree->SetBranchAddress("event_index",&event_index);
	tree->SetBranchStatus("*",kFALSE);// close the other branch
	tree->SetBranchStatus("strip[4608]",kTRUE); // close the other branch

	

	//********************************************
	//********* Reading Fit Parameters ***********
	//********************************************

	// read fit parameters from fitting_parameters files, do NOT delete it
	double p0_read[2][6];
	double p1_read[2][6];
	double p2_read[2][6];

	std::ifstream fitting_param(file_fit_par.c_str());
	for (int o=0; o<6; o++){
		fitting_param >> p0_read[0][o] >> p1_read[0][o] >> p2_read[0][o] >> p0_read[1][o] >> p1_read[1][o] >> p2_read[1][o];
	
	}
	
	//********************************************
	//******** End Reading Fit Parameters ********
	//********************************************

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
		// reads the calibration pedestals, sigma_3s, and non_gaussian_indicies
		std::ifstream calib(file_calib.c_str());
    	std::string linecalib,mean_2_str,sigma_2_str,sigma_3_str, non_gauss_str;
    	int filler = 0;
    	int counter_reading = 0;
    	while (std::getline(calib, linecalib)){
    		counter_reading++;
    		if(counter_reading<7) continue;
    		mean_2_str = "";
    		sigma_2_str = "";
    		sigma_3_str = "";
    		non_gauss_str = "";
    		int dim=0;
    		int contator=0;
    		for( std::string::iterator it = linecalib.begin(); it!=linecalib.end() ; it++){
      			dim ++;
    		}
    		for(int i =0; i<dim-1; i++){
    			int t=0;
    			if (linecalib[i]=='=' && contator==0){ 
    				while (linecalib[i+1+t]!=' '){
	  					mean_2_str+=linecalib[i+1+t] ;
	  					t++;
	  				}
	  				contator++;
	  				i=i+t+1;
	  			}

	  			t=0;
	  			if (linecalib[i]=='=' && contator==1){ 
    				while (linecalib[i+1+t]!=' '){
	  					sigma_2_str+=linecalib[i+1+t] ;
	  					t++;
	  				}
	  				contator++;
	  				i=i+t+1;
	  			}

	  			t=0;
	  			if (linecalib[i]=='=' && contator==2){ 
    				while (linecalib[i+1+t]!=' '){
	  					sigma_3_str+=linecalib[i+1+t] ;
	  					t++;
	  				}
	  				contator++;
	  				i=i+t+1;
	  			}
	  			t=0;
	  			if (linecalib[i]=='=' && contator==3){ 
    				while (linecalib[i+1+t]!=' '){
	  					non_gauss_str+=linecalib[i+1+t] ;
	  					t++;
	  				}
	  				contator++;
	  				i=i+t+1;
	  			}
	  			
    		}

    		mean_2[filler] = float(std::stof( mean_2_str ));
    		sigma_2[filler] = float(std::stof( sigma_2_str));
    		sigma_3[filler] = float(std::stof( sigma_3_str));
    		non_gauss[filler] = float(std::stof( non_gauss_str));
			filler++;
			//if(non_gauss[filler]>gaussian_threshold) deadchannel[filler] = 1;
    	}

		//==============SIGMA_2_VA==============	
		analysis3->cd();
		TH1F * h_sigma_VA;
		for (int t=0; t<n_VA; ++t){ //loop over all the VA
			name = "sigma_VA_"+to_string(t);
			analysis3->cd(name.c_str()); //print histograms in the different folders
			name = title = "sigma_VA_"+to_string(t)+"_"+to_string(s);
			h_sigma_VA = new TH1F(name.c_str(),title.c_str(),200,0,100); // POSSIBLE BIAS
			h_sigma_VA->SetDirectory(analysis3);
			for (int k=0; k<VA_chan; ++k){		//loop over the channel in a single VA
				h_sigma_VA->Fill(sigma_2[t*VA_chan+k]);
			}
			
   			//int bmax = h_sigma_VA->GetMaximumBin(); //you get the number of the maximum bit
   			//max[s][t] = h_sigma_VA->GetBinCenter(bmax); // coordinate of the maximum
			h_sigma_VA->Write();
			max[s][t] = GetCleanedSigma(h_sigma_VA);
			analysis3->cd("/../");
		}
		//?????????????????????????????????____________????????????????????????????????
		int GEcounter[n_steps][n_VA];
		//==================CALCULATION OF THE COMMON NOISE=====================
		for (int j=0; j<n_ev; j++){			// loop over the events
			CommonNoise(ADC2[j], mean_2, sigma_2, max[s], &comnoise_real[j][0], &GEcounter[s][0]/*, deadchannel*/);
			
			for(int t=0; t<n_VA; ++t){
				h_comnoise->Fill(t,comnoise_real[j][t]);
				comnoise_real2[t][j] = comnoise_real[j][t];
			}

		}
		/*
		//==================CALCULATION OF THE SIGMA AFTER CN SUBTRACTION=====================
		for (int i=0; i<n_chann; i++){
			sigma_3[i] = GetSigma(ADC[i],comnoise_real2[i/64],mean_2[i],mean_2[i]-3*sigma_2[i],mean_2[i]+3*sigma_2[i]);
			//h_sigma_3->SetBinContent(i+1,sigma_3[i]);
			h_sigma_tot->Fill(i,sigma_3[i]);
		}
		*/
		
		// distribution of number of seeds INCOMPLETE, TO BE IMPROVED
		/*
		for(int j = 0; j<n_ev; ++j){	
			N_seed(ADC2[j], mean_2, comnoise_real[j], sigma_3, 5., &counter[0];
			for (int w = 0; w < 6; w++){
				h_distribution_seed_5->Fill(counter[w]+20*w);
			}
			N_seed(ADC2[j], mean_2, comnoise_real[j], sigma_3, 4., &counter[0]);
			for (int w = 0; w < 6; w++){
				h_distribution_seed_4->Fill(counter[w]+20*w);
			}
			N_seed(ADC2[j], mean_2, comnoise_real[j], sigma_3, 3., &counter[0];
			for (int w = 0; w < 6; w++){
				h_distribution_seed_3->Fill(counter[w]+20*w);
			}
			N_seed(ADC2[j], mean_2, comnoise_real[j], sigma_3, 3.5, &counter[0]);
			for (int w = 0; w < 6; w++){
				h_distribution_seed_35->Fill(counter[w]+20*w);
			}
		}
		*/

		for(int iev= 0 ; iev<n_ev; ++iev){
			for(int ichan=0;ichan<n_chann;++ichan){
				//if(deadchannel[ichan]!=0) continue;
				h_ADC_tot->Fill(ichan,ADC[ichan][iev]); //<---------- With this lines the processing time increases a lot(x3)
				h_ADC_mean_tot->Fill(ichan,ADC[ichan][iev]-mean_2[ichan]);//<---------- With this lines the processing time increases a lot(x3)
				clearchann[iev][ichan]=ADC[ichan][iev]-mean_2[ichan]-comnoise_real[iev][int(ichan/VA_chan)];
				clearchann2[iev][ichan]=ADC[ichan][iev]-mean_2[ichan];
				h_count->Fill(ichan, clearchann[iev][ichan]);
				signif = clearchann[iev][ichan]/sigma_3[ichan];//<---------- With this lines the processing time increases a lot(x3)
				h_sign->Fill(ichan, signif);//<---------- With this lines the processing time increases a lot (x3)
			}
		}
		for(int iev=0; iev<n_ev; ++iev){
			for(int n_adc=0; n_adc<24; n_adc+=2){
				for(int chan_adc=0; chan_adc<ADC_CHAN; ++chan_adc){
					int truechann = n_adc*ADC_CHAN+chan_adc;
					//std::cout<< truechann<<" "<< truechann+ADC_CHAN << std::endl;
					int side = 0;
	    			if((truechann)%(n_chann/6)>383) side = 1;
					h_correlation_adc[side][int(truechann*6/n_chann)]->Fill(clearchann[iev][truechann],clearchann[iev][truechann+ADC_CHAN]);
				}
			}
		}
		for(int iev=0; iev<n_ev; ++iev){
			for(int n_adc=0; n_adc<24; n_adc++){
				int truechann = n_adc*ADC_CHAN;
				h_corr_adc1chan[n_adc]->Fill(clearchann[iev][0],clearchann[iev][truechann]);
			}
		}
		//cluster finding code
		std::vector< event > clev; 
    	for(int iev=0;iev<NCALIBEVENTS;++iev){
     		std::vector<LTrackerCluster> *clusters=GetClusters(clearchann[iev],sigma_3);
    		event myevent;
     		for(int ev=0;ev<clusters->size();++ev){
				myevent.cls.push_back(clusters->at(ev));
      		}
      	clev.push_back(myevent);
      	}
		/*
		//_____________________ROBERTO CODE__________________________
		std::vector< event > clev; // clever OR clean event
	    for(int iev=0;iev<n_ev;++iev){
	    	event myevent;
	    	myevent.entry = s*n_ev+iev;
	    	//myevent.entry = event_index;
	    	for(int ichan=0;ichan<n_chann;++ichan){
	    		//if(deadchannel[ichan]!=0) continue;
	    		double count=(ADC[ichan][iev]-mean_2[ichan]-comnoise_real[iev][ichan/VA_chan]);
	    		double sigma=sigma_3[ichan];
	   			double sign= count/sigma;
	    		if(sign>sigma_cut){
		 		    // Open window for cluster identification *******
		    		std::vector< double > window_count;
		    		std::vector< double > window_sigma;
		    		std::vector< double > window_sign;
		    		for(int iw=0; iw<7; ++iw) {
		     			double mycount = (ADC[ichan+iw][iev]-mean_2[ichan+iw]-comnoise_real[iev][(ichan+iw)/VA_chan]);
		      			double mysigma = sigma_3[ichan+iw];
		     			window_count.push_back(mycount);
		     			window_sigma.push_back(mysigma);
		  			    window_sign.push_back(mysigma != 0. ? mycount/mysigma : 0.);
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

		    		int max2ind=-999;
		    		double max2sign=-999.;
		    		if(wsize!=1){
			    		for(int il=0; il<wsize; ++il) {
			    			if(window_sign.at(il)>max2sign && il!= maxind) {
			        			max2ind=il;
			        			max2sign=window_sign.at(il);
			    			}
			    		}
			    	}
		    		else {
		    			max2ind=maxind;
		    		}
		    		// Seed your life... *****************************
		   			int seed = -999;
		    		// Case for dead channel in between
		   			cluster mycluster;

		   			//Previus selection of seed, but incloplete or wrong
		   			
				    if(abs(max2ind-maxind)==2 && max2sign>3.) {
				    	//if(window_sign.at((maxind+max2ind)/2)<3.) mycluster.goodorbad = 1;
				    	seed=ichan+(maxind+max2ind)/2;
				    }
				    else if(abs(maxind-max2ind)==1) seed = maxind+ichan;
				    else seed=(max2ind < maxind ? max2ind+ichan : maxind+ichan); 
				    

				    seed = maxind + ichan;

		    		// ***********************************************
		    		// Fill cluster data *****************************
				    mycluster.seed = seed;
		    		for(int iii=0; iii<5; ++iii) {
		    			int cchan = seed-2+iii;//seed-ichan-1+iii;
		    			bool SameLP_FLAG = SameLadderPlane(seed, cchan);
		    			
		    			if(deadchannel[cchan]!=0){
		    				mycluster.count[iii] = 0.;
		    				mycluster.sign[iii] = 0.;
		    				mycluster.sigma[iii] = 0.;
		    			}
		    			
		    			//else{
		    				mycluster.count[iii] = ( SameLP_FLAG ? clearchann[cchan][iev] : 0.); //window_count.at(cchan) : 0.);
		    				mycluster.sign[iii] = ( SameLP_FLAG ? clearchann[cchan][iev]/sigma_3[cchan] : 0.); //window_sign.at(cchan) : 0.);
		    				mycluster.sigma[iii] = ( SameLP_FLAG ? sigma_3[cchan] : 0.); //window_sigma.at(cchan) : 0.);
		    			//}
		    		}

		    		myevent.cls.push_back(mycluster);

		    		ichan = seed+5; // very important! Restart after the cluster
	    		}
	      	}
			
			calib.close();
	      	if(myevent.cls.size()>0) clev.push_back(myevent);
	    }
		//________________________.:End Roberto Code:._______________________________
	    */
	    for(int i=0; i<int(clev.size()); ++i){
	    	for(int h=0; h<int(clev.at(i).cls.size()); ++h){
    			LTrackerCluster mycl = clev.at(i).cls.at(h);
    			int side = 0;
    			if(mycl.seed%(n_chann/6)>383) side = 1;
	    		h_good_events->Fill(mycl.seed);
	    		h_correlation->Fill(mycl.count[2],mycl.count[1]+mycl.count[3]);
	    		h_correlation_plus_minus->Fill(mycl.count[1],mycl.count[3]);
	    		h_correlation_minus->Fill(mycl.count[2],mycl.count[1]);
	    		h_correlation_plus->Fill(mycl.count[2],mycl.count[3]);
	    		h_correlation_sign->Fill(mycl.sn[2], sqrt(pow(mycl.sn[1],2)+pow(mycl.sn[3],2)));
	    		h_correlation_sign_plus_minus->Fill(mycl.sn[1],mycl.sn[3]);
	    		h_correlation_sign_minus->Fill(mycl.sn[2], mycl.sn[1]);
				h_correlation_sign_plus->Fill(mycl.sn[2], mycl.sn[3]);
				h_countseed[side][int(mycl.seed*6/n_chann)]->Fill(mycl.count[2]);
				h_signseed[side][int(mycl.seed*6/n_chann)]->Fill(mycl.sn[2]);	
				//h_sumcontcluster[side][int(mycl.seed*6/n_chann)]->Fill(mycl.count[1]+mycl.count[2]+mycl.count[3]);
				h_sumsigncluster[side][int(mycl.seed*6/n_chann)]->Fill(sqrt(pow(mycl.sn[1],2)+pow(mycl.sn[2],2)+pow(mycl.sn[3],2)));
				h_ratiocontcluster[side][int(mycl.seed*6/n_chann)]->Fill(mycl.count[2],(mycl.count[1]+mycl.count[3])/mycl.count[2]);
				h_ratiocontcluster_minus[side][int(mycl.seed*6/n_chann)]->Fill(mycl.count[2],mycl.count[1]/mycl.count[2]);
				h_ratiocontcluster_plus[side][int(mycl.seed*6/n_chann)]->Fill(mycl.count[2],mycl.count[3]/mycl.count[2]);
				h_ratiosigncluster[side][int(mycl.seed*6/n_chann)]->Fill(mycl.sn[2],sqrt(pow(mycl.sn[1],2)+pow(mycl.sn[3],2))/mycl.sn[2]);
				h_ratiosigncluster_minus[side][int(mycl.seed*6/n_chann)]->Fill(mycl.sn[2],mycl.sn[1]/mycl.sn[2]);
				h_ratiosigncluster_plus[side][int(mycl.seed*6/n_chann)]->Fill(mycl.sn[2],mycl.sn[3]/mycl.sn[2]);
				
				double barycenter=0.;
				double barycenter_abs=0.;
				int counter_bar=0;
				double sum_count_bar=0.;
				double eta = -999.;
				for (int l=0; l<5; l++){
					if (mycl.sn[l]>3.){
						barycenter += l*mycl.count[l];
						barycenter_abs += (l+mycl.seed)*mycl.count[l];
						sum_count_bar += mycl.count[l];
						counter_bar++;
					}
				}
				//if(mycl.goodorbad==0){
				eta = clev.at(i).cls.at(h).GetEta();
					//eta = ( mycl.count[1]>mycl.count[3] ? (mycl.count[2]-mycl.count[1])/(mycl.count[2]+mycl.count[1]) : (mycl.count[3]-mycl.count[2])/(mycl.count[3]+mycl.count[2]));
				//}
				h_clustersize[side][int(mycl.seed*6/n_chann)]->Fill(counter_bar);
				h_sumcontcluster_uncorrected[side][int(mycl.seed*6/n_chann)]->Fill(sum_count_bar);
				h_sumcontcluster[side][int(mycl.seed*6/n_chann)]->Fill(sum_count_bar*Correction_1(eta,p0_read[side][int(mycl.seed*6/n_chann)],p1_read[side][int(mycl.seed*6/n_chann)],p2_read[side][int(mycl.seed*6/n_chann)]));
				h_sumcontcluster_corr2[side][int(mycl.seed*6/n_chann)]->Fill(sum_count_bar*Correction_2(getmyprofile[side][int(mycl.seed*6/n_chann)],eta));
				h_sumcontcluster_corr3[side][int(mycl.seed*6/n_chann)]->Fill(sum_count_bar*Correction_2(getmyprofile2[side][int(mycl.seed*6/n_chann)],eta));
				barycenter/=sum_count_bar;
				barycenter_abs/=sum_count_bar;
				h_cls_barycenter[side][int(mycl.seed*6/n_chann)]->Fill(barycenter);
				//h_cls_barycenter_abs[side][int(mycl.seed*6/n_chann)]->Fill(barycenter_abs);
				h_eta_distribution[side][int(mycl.seed*6/n_chann)]->Fill(eta);
				h_sumcontcluster_eta[side][int(mycl.seed*6/n_chann)]->Fill(eta,sum_count_bar);
				h_sumcontcluster_eta_prof[side][int(mycl.seed*6/n_chann)]->Fill(eta,sum_count_bar);
				h_myprofile[side][int(mycl.seed*6/n_chann)]->Fill(eta,sum_count_bar);


	    	}
	    }
	   
		// calculate time each step
	    clock_t stop = clock();
		std::cout<<"step [" <<int(s+1)<<"/"<<n_steps<< "] time = " <<((double)stop - (double)start)/CLOCKS_PER_SEC<<std::endl;
		
		//It sums all the sigmas and means each steps
		for(int i=0; i<n_chann; i++){
			mean_final[i] += mean_2[i];
			sigma_final[i] += sigma_3[i];
		}

		for(int iev=0; iev<n_ev; iev++){
			for(int i=0; i<n_chann; i++){
				clearchanntot2[i]+=clearchann2[iev][i];
			}
		}
		for(int iev=0; iev<n_ev; iev++){
			for(int i=0; i<n_chann; i++){
				clearchanntot[i]+=clearchann[iev][i];
			}
		}

		std::ofstream sigmatxt;
		name = "sigma"+typeofparticle+".txt";
		sigmatxt.open(name);
		for(int i=0; i<n_chann;i++){
			sigmatxt<<sigma_3[i]<<std::endl;
		}
	}

	
	analysis3->Close();
	//sigma_mean and mean_mean for all the steps
	for(int i=0; i<n_chann; i++){
		sigma_mean->SetBinContent(i+1,sigma_final[i]/n_steps);
		mean_mean->SetBinContent(i+1,mean_final[i]/n_steps);
		h_non_gauss->SetBinContent(i+1,non_gauss[i]);
	}
	for(int side=0; side<2; side ++){
		for(int w=0; w<6; w++){
			double nent_cluster_size = h_clustersize[side][w]->GetEntries();
			h_clustersize[side][w]->Scale(1/nent_cluster_size);
		}
	}
	

				//*****************************************************************
				//**********************Get Fit Parameters ************************   //ATTENCTION
				//*****************************************************************
	
	double p0[2][6];
	double p1[2][6];
	double p2[2][6];
	Fitting_func(h_sumcontcluster_eta_prof[0], &p0[0][0], &p1[0][0], &p2[0][0]);
	Fitting_func(h_sumcontcluster_eta_prof[1], &p0[1][0], &p1[1][0], &p2[1][0]);

	double langau_MPV[2][6];
	double langau_sigma[2][6];
	Langau(h_sumcontcluster[0], &langau_MPV[0][0], &langau_sigma[0][0]);
	Langau(h_sumcontcluster[1], &langau_MPV[1][0], &langau_sigma[1][0]);
	/*
	std::ofstream fitting_paremters;
	name = "fitting_parameters_"+typeofparticle+".txt";
	fitting_paremters.open(name);

	for(int i_par=0; i_par<6; i_par++){
		fitting_paremters<<p0[0][i_par]<<" "<<p1[0][i_par]<<" "<<p2[0][i_par]<<" "<<p0[1][i_par]<<" "<<p1[1][i_par]<<" "<<p2[1][i_par]<<std::endl;
	}
	*/
	
	std::ofstream langau_paremters;
	name = "langau_param.txt";
	langau_paremters.open(name, fstream::app);

	for(int i_par=0; i_par<6; i_par++){
		langau_paremters<<langau_MPV[0][i_par]<<" "<<langau_sigma[0][i_par]<<" "<<langau_MPV[1][i_par]<<" "<<langau_sigma[1][i_par]<<std::endl;
	}
	std::ofstream furier;
	name = "furier_"+typeofparticle+".txt";
	furier.open(name);
	for(int i=0; i<n_chann; i++){
		if(non_gauss[i]<MASK[i]){
			furier<</*"channel = "<<i<<" clearchanncount = " <<*/clearchanntot[i]/*<<" non gaus index "<< non_gauss[i]<<" MASK[i] "<<MASK[i]*/<<std::endl;
		}
		else{
			furier<</*"channel = "<<i<<" clearchanncount = " <<*/0./*<<" non gaus index "<< non_gauss[i]<<" MASK[i] "<<MASK[i]*/<<std::endl;
		}
	}
	std::ofstream furierprime;
	name = "furierprime_"+typeofparticle+".txt";
	furierprime.open(name);
	for(int i=0; i<n_chann; i++){
		if(non_gauss[i]<MASK[i]){
			furierprime<</*"channel = "<<i<<" clearchanncount = " <<*/clearchanntot2[i]/*<<" non gaus index "<< non_gauss[i]<<" MASK[i] "<<MASK[i]*/<<std::endl;
		}
		else{
			furierprime<</*"channel = "<<i<<" clearchanncount = " <<*/0./*<<" non gaus index "<< non_gauss[i]<<" MASK[i] "<<MASK[i]*/<<std::endl;
		}
	}

	std::ofstream masktxt;
	name = "mask.txt";
	masktxt.open(name);
	for(int i=0; i<n_chann; i++){
		if(non_gauss[i]<MASK[i]){
			masktxt<<1<<std::endl;
		}
		else{
			masktxt<<0<<std::endl;
		}
	}
	


					//========================================================================
					//=====================WRITING & PRINTING RESULTS=========================
					//========================================================================

	
	analysis2->cd();
	h_comnoise->Write();
	h_sigma_tot_prime->Write();
	h_sigma_tot_second->Write();
	h_sigma_tot->Write();
	h_ADC_tot->Write();
	h_mean_tot->Write();
	h_ADC_mean_tot->Write();
	h_good_events->Write();
	h_sign->Write();
	h_count->Write();
	h_correlation->Write();
	h_correlation_minus->Write();
	h_correlation_plus->Write();
	h_correlation_sign->Write();
	sigma_mean->Write();
	mean_mean->Write();
	for (int ilad = 0; ilad<6; ilad++){
		checkped[ilad]->Write();
		for(int side = 0; side < 2; side++){
			h_sumsigncluster[side][ilad]->Write();
			h_sumcontcluster[side][ilad]->Write();
			h_sumcontcluster_uncorrected[side][ilad]->Write();
			h_sumcontcluster_corr2[side][ilad]->Write();
			h_sumcontcluster_corr3[side][ilad]->Write();
			h_ratiocontcluster[side][ilad]->Write();
			h_ratiocontcluster_plus[side][ilad]->Write();
			h_ratiocontcluster_minus[side][ilad]->Write();
			h_ratiosigncluster[side][ilad]->Write();
			h_ratiosigncluster_plus[side][ilad]->Write();
			h_ratiosigncluster_minus[side][ilad]->Write();
			h_countseed[side][ilad]->Write();
			h_signseed[side][ilad]->Write();
			h_clustersize[side][ilad]->Write();
			h_myprofile[side][ilad]->Write();
			h_eta_distribution[side][ilad]->Write();
			h_sumcontcluster_eta[side][ilad]->Write();
			h_sumcontcluster_eta_prof[side][ilad]->Write();
			h_correlation_adc[side][ilad]->Write();
		}
	}
	h_correlation_sign_plus_minus->Write();
	h_correlation_plus_minus->Write();
	h_non_gauss->Write();
	h_correlation_sign_minus->Write();
	h_correlation_sign_plus->Write();
	/*
	h_distribution_seed_5->Write();
	h_distribution_seed_4->Write();
	h_distribution_seed_35->Write();
	h_distribution_seed_3->Write();
	*/
	h_checkdeadchannel->Write();
	checkRobCode->Write();
	analysis2->Close();

	// Number of seed distributions overlapped
	gROOT->SetBatch(kTRUE);
	/*
	TCanvas * c5 = new TCanvas();
	c5->Print("NSeedDist_calibration_and_data.pdf[");
		Plot6_1d(h_distribution_seed_5,"h_distribution_seed_5","number of seed with 5 #sigma cut","entries",0,60000,"no")->Print("NSeedDist_calibration_and_data.pdf");
		Plot6_1d(h_distribution_seed_4,"h_distribution_seed_4","number of seed with 4 #sigma cut","entries",0,60000,"no")->Print("NSeedDist_calibration_and_data.pdf");		
		Plot6_1d(h_distribution_seed_35,"h_distribution_seed_35","number of seed with 3.5 #sigma cut","entries",0,60000,"no")->Print("NSeedDist_calibration_and_data.pdf");
		Plot6_1d(h_distribution_seed_3,"h_distribution_seed_3","number of seed with 3 #sigma cut","entries",0,60000,"no")->Print("NSeedDist_calibration_and_data.pdf");
	c5->Print("NSeedDist_calibration_and_data.pdf]");
	*/
	string name_plot1 = typeofparticle+"_calibration_and_data_1.pdf[";
	string name_plot11 = typeofparticle+"_calibration_and_data_1.pdf";
	string name_plot111 = typeofparticle+"_calibration_and_data_1.pdf]";
	TCanvas * c1 = new TCanvas();
	c1->Print(name_plot1.c_str());
		Plot6_1d(h_non_gauss,"h_non_gauss","channel","non_gauss_index",-5,10,"yes")->Print(name_plot11.c_str());
		Plot6_2d(h_comnoise,"h_comnoise","n_VA","ComNoise",-100,100,"yes","")->Print(name_plot11.c_str());
		//Plot6_2d(h_sigma_tot_prime,"h_sigma_tot_prime","channel","#sigma_raw",0,200,"yes","")->Print(name_plot11.c_str());
		//Plot6_2d(h_sigma_tot_second,"h_sigma_tot_second","channel","#sigma_raw",0,200,"yes","")->Print(name_plot11.c_str());
		//Plot6_2d(h_sigma_tot,"h_sigma_tot","channel","#sigma_raw",0,40,"yes","")->Print(name_plot11.c_str());
		Plot6_2d(h_ADC_tot,"h_ADC_tot","channel","ADC",0,3000,"yes","")->Print(name_plot11.c_str());
		//Plot6_2d(h_mean_tot,"h_mean_tot","channel","Ped",0,3000,"yes","")->Print(name_plot11.c_str());
		Plot6_2d(h_ADC_mean_tot,"h_ADC_mean_tot","channel","ADC-Ped",-500,500,"yes","")->Print(name_plot11.c_str());
		Plot6_2d(h_count,"h_count","channel","events",-500,500,"yes","")->Print(name_plot11.c_str());
		Plot6_2d(h_sign,"h_sign","channel","sign",-50,50,"yes","5sign")->Print(name_plot11.c_str());
		Plot6_1d(h_good_events,"h_good_events","channel","events",0,400,"yes")->Print(name_plot11.c_str());
		Plot_2d(h_correlation,"ADC_{seed}","ADC_{seed-1}+ADC_{seed+1}")->Print(name_plot11.c_str());
		Plot_2d(h_correlation_plus,"ADC_{seed}","ADC_{seed+1}")->Print(name_plot11.c_str());
		Plot_2d(h_correlation_plus_minus,"ADC_{seed-1}","ADC_{seed+1}")->Print(name_plot11.c_str());
		Plot_2d(h_correlation_sign_plus_minus,"sign_{seed-1}","sign_{seed+1}")->Print(name_plot11.c_str());
		Plot_2d(h_correlation_minus,"ADC_{seed}","ADC_{seed-1}")->Print(name_plot11.c_str());
		Plot_2d(h_correlation_sign_minus,"sign_{seed}","sign_{seed-1}")->Print(name_plot11.c_str());
		Plot_2d(h_correlation_sign_plus,"sign_{seed}","sign_{seed+1}")->Print(name_plot11.c_str());
		Plot_2d(h_correlation_sign,"sign_{seed}","#sqrt(sign_{seed+1}^2+sign_{seed-1}^2)")->Print(name_plot11.c_str());
		Plot_2d(h_checkdeadchannel,"h_checkdeadchannel","channel")->Print(name_plot11.c_str());
		for(int ilad = 0; ilad<6; ++ilad){
			checkpedname = "checkpedname_ladder_"+to_string(ilad);
			Plot_1d(checkped[ilad],"step","mean_ped")->Print(name_plot11.c_str());
		}

	c1->Print(name_plot111.c_str());
	string name_plot2 = typeofparticle+"_calibration_and_data_2.pdf[";
	string name_plot22 = typeofparticle+"_calibration_and_data_2.pdf";
	string name_plot222 = typeofparticle+"_calibration_and_data_2.pdf]";
	TCanvas * c3 = new TCanvas();
	c3->Print(name_plot2.c_str());
		for (int side = 0; side <2; side++){
			Plot_6histo_1d(h_clustersize[side],"clustersize","entries","")->Print(name_plot22.c_str());
			Plot_6histo_1d(h_cls_barycenter[side],"barycenter","entries","")->Print(name_plot22.c_str());
			Plot_6histo_1d(h_countseed[side],"count seed","entries","log")->Print(name_plot22.c_str());
			Plot_6histo_1d(h_signseed[side],"sign seed","entries","log")->Print(name_plot22.c_str());
			Plot_6histo_1d(h_eta_distribution[side],"eta","count","")->Print(name_plot22.c_str());
			Plot_6histo_2d(h_sumcontcluster_eta[side],"sum counts","entries")->Print(name_plot22.c_str());
			Plot_6histo_profile(h_sumcontcluster_eta_prof[side],"eta","count","")->Print(name_plot22.c_str());
			Plot_6histo_profile(getmyprofile[side],"eta","count","")->Print(name_plot22.c_str());
			Plot_6histo_profile(getmyprofile2[side],"eta","count","")->Print(name_plot22.c_str());
			Plot_6histo_profile(h_myprofile[side],"eta","count","")->Print(name_plot22.c_str());
			Plot_6histo_1d(h_sumcontcluster_uncorrected[side],"sum counts","entries","")->Print(name_plot22.c_str());
			Plot_6histo_1d(h_sumcontcluster[side],"sum counts","entries","")->Print(name_plot22.c_str());
			Plot_6histo_1d(h_sumcontcluster_corr2[side],"sum counts","entries","")->Print(name_plot22.c_str());
			Plot_6histo_1d(h_sumcontcluster_corr3[side],"sum counts","entries","")->Print(name_plot22.c_str());
			Plot_6histo_1d(h_sumsigncluster[side],"sum sign","entries","log")->Print(name_plot22.c_str());
			Plot_6histo_2d(h_ratiocontcluster[side],"ADC_{seed}","(ADC_{seed-1}+ADC_{seed+1})/ADC_{seed}")->Print(name_plot22.c_str());
			Plot_6histo_2d(h_ratiocontcluster_minus[side],"ADC_{seed}","ADC_{seed-1}/ADC_{seed}")->Print(name_plot22.c_str());
			Plot_6histo_2d(h_ratiocontcluster_plus[side],"ADC_{seed}","ADC_{seed+1}/ADC_{seed}")->Print(name_plot22.c_str());
			Plot_6histo_2d(h_ratiosigncluster[side],"sign_{seed}","(sign_{seed-1}+sign_{seed+1})/sign_{seed}")->Print(name_plot22.c_str());
			Plot_6histo_2d(h_ratiosigncluster_minus[side],"sign_{seed}","sign_{seed-1}/sign_{seed}")->Print(name_plot22.c_str());
			Plot_6histo_2d(h_ratiosigncluster_plus[side],"sign_{seed}","sign_{seed+1}/sign_{seed}")->Print(name_plot22.c_str());
			Plot_6histo_2d(h_correlation_adc[side],"count_{ADC}","count_{ADC+1}")->Print(name_plot22.c_str());
		}
	c3->Print(name_plot222.c_str());

	string name_plot3 = typeofparticle+"_calibration_and_data_3.pdf[";
	string name_plot33 = typeofparticle+"_calibration_and_data_3.pdf";
	string name_plot333 = typeofparticle+"_calibration_and_data_3.pdf]";
	TCanvas * c4 = new TCanvas();
	c4->Print(name_plot3.c_str());
	for (int w = 0; w <24; w++){
		Plot_2d(h_corr_adc1chan[w],"1st chan ADC","1st chan ADC+1")->Print(name_plot33.c_str());
	}
	c4->Print(name_plot333.c_str());
}



int main(){
	analysis_follega("20161112-112421-Run_3C_37MeV_SERVO_EASIROC2.root","protons37MeV_c5","calib_protons_37MeV.txt","fitting_parameters_protons37MeV_c5.txt");
	analysis_follega("RUN_3C_51MeV_SERVO_EASIROC2_HOT.root","protons50MeV_c5","calib_protons_50MeV.txt","fitting_parameters_protons50MeV_c5.txt");
	analysis_follega("data_proton_11Nov_4C.root","protons70MeV_c5","calib_protons_70MeV.txt","fitting_parameters_protons70MeV_c5.txt");
	analysis_follega("RUN_3C_100MeV_SERVO_EASIROC2_HOT.root","protons100MeV_c5","calib_protons_100MeV.txt","fitting_parameters_protons100MeV_c5.txt");
	analysis_follega("RUN_3C_125MeV_SERVO_EASIROC2_HOT.root","protons125MeV_c5","calib_protons_125MeV.txt","fitting_parameters_protons125MeV_c5.txt");
	analysis_follega("RUN_3C_154MeV_SERVO_EASIROC2_HOT.root","protons154MeV_c5","calib_protons_154MeV.txt","fitting_parameters_protons154MeV_c5.txt");
	analysis_follega("RUN_3C_174MeV_SERVO_EASIROC2_HOT.root","protons174MeV_c5","calib_protons_174MeV.txt","fitting_parameters_protons174MeV_c5.txt");
	analysis_follega("data_proton_12Nov_3C.root","protons228MeV_c5","calib_protons_228MeV.txt","fitting_parameters_protons228MeV_c5.txt");
	//Silicon_analysis_cls5("RUN_ROTAZIONE_70MeV_SERVO_EASIROC2_HOT.root","protons_beamrotated","calib_protons_beamrotated.txt","fitting_parameters_beamrotated_c5.txt");

	return 0;
}


		/*
		//==========================================
		//================= OLD CODE ===============
		//==========================================
		float cont = 0;
		string name;
		string title;
		string address;
	
		address = "strip["+to_string(i)+"]";
		====================Drawing histograms=====================
		TH1S * hist1D[4608]; 
		TFile * analysis = new TFile("analysis.root","RECREATE");
		std::cout<<"channel : "<<i<<" event : "<<j<<" "<<strip[i]<<std::endl;
		name = "1step_histo" + to_string(i);
		title = "1step_histo" + to_string(i);
		hist1D[i] = new TH1S(name.c_str(),title.c_str(),4096,0,4096);
		address = "strip["+to_string(i)+"]>>histo(4096,0,4096)";
		tree->Draw(address.c_str(),"event_index>=1000&&event_index<=2000");
		hist1D[i] = (TH1S*)gDirectory->Get("histo");
		analysis->cd();
		hist1D[i]->Write(name.c_str());
		std::cout<<"percentage drawing ["<<int(i*100/4608.)<<"%]"<<std::endl;
		
		//=====================ANALYSIS_STEP1========================
		hist1D[i]->Fit("gaus","q");
		mean_1[i] = hist1D[i]->GetMean();
		sigma_1[i] = hist1D[i]->GetStdDev();
		h_mean_1->SetBinContent(i,mean_1[i]);
		h_sigma_1->SetBinContent(i,sigma_1[i]);
		
		//=====================ANALYSIS_STEP2========================
		hist1D[i]->Fit("gaus","q","",mean_1[i]-3*sigma_1[i],mean_1[i]+3*sigma_1[i]);
		mean_2[i] = hist1D[i]->GetMean();
		sigma_2[i] = hist1D[i]->GetStdDev();
		h_mean_2->SetBinContent(i,mean_2[i]);
		h_sigma_2->SetBinContent(i,sigma_2[i]);	

		std::cout<<"percentage drawing ["<<int(100)<<"%]"<<std::endl;
		
		h_mean_1->Write();
		h_sigma_1->Write();
		h_mean_2->Write();
		h_sigma_2->Write();
		


	
		//===============SELECTION WITH THE THRESHOLD==============
		if (sigma_1[i]>threshold){
			std::cout<<"above"<<std::endl;
			goodorbad[i]=1;				//---> BAD CHANNELS
		}

			
		for (int i=0; i<4608; i++){
			sigma_2[i]+=(strip[i]-mean_1[i]-comnoise_real[i])*(strip[i]-mean_1[i]-comnoise_real[i]);
		}
		//std::cout<<sigma_2[i]<<std::endl;
		
		//std::cout<<"percentage of caluculation of the common noise ["<<int(cont*100/(4608000*2.))<<"%]"<<std::endl;
	}
	
	for (int i = 0; i < 4608; i++){
		sigma_2[i]=(1/n_ev)*sigma_2[i];
		h_sigma_2->SetBinContent(i,sigma_2[i]);
	}
	

	
		if (s==0){
			analysis->cd();
			h_mean_1->Write();
			h_sigma_1->Write();
			h_mean_2->Write();
			h_sigma_2->Write();
			h_sigma_3->Write();
			analysis->Close();
		}
		
		
		for (int i=0;i<n_chann;++i){
			for (int j=0;j<1000;++j){
				content = content_minus = content_plus = 0;
				h_ADC_tot->Fill(i,ADC[i][j]);
				h_ADC_mean_tot->Fill(i,ADC[i][j]-mean_2[i]);
				content = (ADC[i][j] - mean_2[i] - comnoise_real[j][i/64]);
				clearchann[i][j] = content;
				//std::cout<<clearchann[i][j]<<std::endl;
				content_minus = (i>0 ? ADC[i-1][j] - mean_2[i-1] - comnoise_real[j][(i-1)/64] : 0.);
				content_plus = (i<4607 ? ADC[i+1][j] - mean_2[i+1] - comnoise_real[j][(i+1)/64] : 0.);
				h_count->Fill(i, content);
				significance = content/sigma_3[i];
				significance_minus = 0.;
				significance_plus = 0.;
				if ((i+1)%384!=0 || i==0 ) significance_plus = content_plus/sigma_3[i+1];
				if (i%384!=0) significance_minus = content_minus/sigma_3[i-1];
				h_sign->Fill(i, significance);
				if (significance>5.){
					h_good_events->Fill(i);
					h_correlation_minus->Fill(content,content_minus);
					h_correlation_sign_minus->Fill(significance, significance_minus);
					h_correlation_sign_plus->Fill(significance, significance_plus);
					h_correlation->Fill(content,content_plus+content_minus);
					h_correlation_plus->Fill(content,content_plus);
					h_correlation_sign->Fill(significance,sqrt(significance_minus*significance_minus+significance_plus*significance_plus));

				}
			}
		}




// experiment of tracking, now ester does it better
	std::ofstream printout;
	string printout_name;
	printout_name = "printout"+typeofparticle+"_golden.txt";
	printout.open(printout_name);
	std::ofstream printout_triple;
	string printout_name_triple;
	printout_name_triple = "printout"+typeofparticle+"_triple.txt";
	printout_triple.open(printout_name_triple);
	std::ofstream printout_double;
	string printout_name_double;
	printout_name_double = "printout"+typeofparticle+"_double.txt";
	printout_double.open(printout_name_double);
	





 
		for(int i=0; i<int(clev.size()); ++i){

	    	int lladder = 4;
	    	int sside = 0;
	    	int maxmax_pos = -999;
	    	float maxcount = 0.;
	    	float maxmax = 0.;
	    	int checker = 0;
	    	for(int h=0; h<int(clev.at(i).cls.size()); ++h){
	    		int pos_seed = clev.at(i).cls.at(h).seed;
	    		
	    		if((clev.at(i).cls.at(h).sign[2] > maxmax)&&(ChanToSide(pos_seed)==sside && ChanToLadder(pos_seed)==lladder)){
	    			maxmax_pos = clev.at(i).cls.at(h).seed;
	    			maxmax = clev.at(i).cls.at(h).sign[2];
	    			maxcount = clev.at(i).cls.at(h).count[2];
	    		}
	    	}
			int seed_p_4 = maxmax_pos;
			if (seed_p_4==-999) checker ++;
			int sign_p_4 = maxmax;
			int count_p_4 = maxcount;

			lladder = 4;
		   	sside = 1;
	    	maxmax_pos = -999;
	    	maxmax = 0.;
	    	maxcount = 0.;
	    	for(int h=0; h<int(clev.at(i).cls.size()); ++h){
	    		int pos_seed = clev.at(i).cls.at(h).seed;
	    		
	    		if((clev.at(i).cls.at(h).sign[2] > maxmax)&&(ChanToSide(pos_seed)==sside && ChanToLadder(pos_seed)==lladder)){
	    			maxmax_pos = clev.at(i).cls.at(h).seed;
	    			maxmax = clev.at(i).cls.at(h).sign[2];
	    			maxcount = clev.at(i).cls.at(h).count[2];
	    		}
	    	}
			int seed_n_4 = maxmax_pos;
			if (seed_n_4==-999) checker ++;
			int sign_n_4 = maxmax;
			int count_n_4 = maxcount;

			lladder = 5;
	   		sside = 0;
	    	maxmax_pos = -999;
	    	maxmax = 0.;
	    	maxcount = 0.;
	    	for(int h=0; h<int(clev.at(i).cls.size()); ++h){
	    		int pos_seed = clev.at(i).cls.at(h).seed;
	    		
	    		if((clev.at(i).cls.at(h).sign[2] > maxmax)&&(ChanToSide(pos_seed)==sside && ChanToLadder(pos_seed)==lladder)){
	    			maxmax_pos = clev.at(i).cls.at(h).seed;
	    			maxmax = clev.at(i).cls.at(h).sign[2];
	    			maxcount = clev.at(i).cls.at(h).count[2];
	    		}
	    	}
			int seed_p_5 = maxmax_pos;
			if (seed_p_5==-999) checker ++;
			int sign_p_5 = maxmax;
			int count_p_5 = maxcount;

			lladder = 5;
	   		sside = 1;
	   	 	maxmax_pos = -999;
	    	maxmax = 0.;
	    	maxcount = 0.;
	    	for(int h=0; h<int(clev.at(i).cls.size()); ++h){
	    		int pos_seed = clev.at(i).cls.at(h).seed;
	    		
	    		if((clev.at(i).cls.at(h).sign[2] > maxmax)&&(ChanToSide(pos_seed)==sside && ChanToLadder(pos_seed)==lladder)){
	    			maxmax_pos = clev.at(i).cls.at(h).seed;
	    			maxmax = clev.at(i).cls.at(h).sign[2];
	    			maxcount = clev.at(i).cls.at(h).count[2];
	    		}
	    	}
			int seed_n_5 = maxmax_pos;
			if (seed_n_5==-999) checker ++;
			int sign_n_5 = maxmax;
			int count_n_5 = maxcount;
			
			if(checker==0){
				printout<<clev.at(i).entry << " " << seed_p_4 << " " << seed_n_4 << " " << seed_p_5 << " " << seed_n_5 <<
										      " " << sign_p_4 << " " << sign_n_4 << " " << sign_p_5 << " " << sign_n_5 << 
										      " " << count_p_4 << " " << count_n_4 << " " << count_p_5 << " " << count_n_5 <<std::endl;
								   		  			   		 
			}
			if(checker==1){
				printout_triple<<clev.at(i).entry << " " << seed_p_4 << " " << seed_n_4 << " " << seed_p_5 << " " << seed_n_5 <<
										      " " << sign_p_4 << " " << sign_n_4 << " " << sign_p_5 << " " << sign_n_5 << 
										      " " << count_p_4 << " " << count_n_4 << " " << count_p_5 << " " << count_n_5 <<std::endl;
								   		  			   		 
			}
			if(checker==2){
				printout_double<<clev.at(i).entry << " " << seed_p_4 << " " << seed_n_4 << " " << seed_p_5 << " " << seed_n_5 <<
										      " " << sign_p_4 << " " << sign_n_4 << " " << sign_p_5 << " " << sign_n_5 << 
										      " " << count_p_4 << " " << count_n_4 << " " << count_p_5 << " " << count_n_5 <<std::endl;
								   		  			   		 
			}

	    
		}
		*/






		




