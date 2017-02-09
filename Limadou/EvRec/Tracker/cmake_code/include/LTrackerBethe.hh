#ifndef __LTRACKERBETHE__
#define __LTRACKERBETHE__ 1

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


const int n_points = 7;

double delta(Double_t * x, double x_0, double x_1, double Cbar, double a, double k, double delta_0 );
double Ftau(Double_t * x, double mass);
void Dummy();
float Conv(double x, double mass);
Double_t bethefuncion(Double_t *x, Double_t *par);
void bethe();
double delta(Double_t * x, double x_0, double x_1, double Cbar, double a, double k, double delta_0 );
Double_t Bethe_p(Double_t *x, Double_t *par);
TCanvas * Plotter(TGraphErrors * graph);
void langau_fit();

#endif
