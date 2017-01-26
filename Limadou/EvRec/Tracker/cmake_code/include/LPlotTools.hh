#ifndef __LPLOTTOOLS__
#define __LPLOTTOOLS__ 1

#include "TCanvas.h"
#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include <string>
#include "string.h"

TCanvas * Plot_2d(TH2F * histo, string xaxis, string yaxis);
TCanvas * Plot_1d(TH1F * histo, string xaxis, string yaxis);
TCanvas * Plot6_2d(TH2F * histo, string name, string xaxis, string yaxis, float min, float max, string cond, string cond2);
TCanvas * Plot6_1d(TH1F * histo, string name, string xaxis, string yaxis ,float min, float max, string cond);
TCanvas * Plot_6histo_1d(TH1F * histo[6], string xaxis, string yaxis, string log);
TCanvas * Plot_6histo_profile(TProfile * histo[6], string xaxis, string yaxis, string log);
TCanvas * Plot_6histo_2d(TH2F * histo[6], string xaxis, string yaxis);


#endif

