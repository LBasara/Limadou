#ifndef LIMADOU_CONST
#define LIMADOU_CONST 1

const int VA_CHAN=64;
const int ADC_CHAN=3*VA_CHAN;
const int LADDER_CHAN=4*ADC_CHAN;
const int SIDE_CHAN=LADDER_CHAN/2;
const int NCHAN = 6*LADDER_CHAN;
const int NADC = 4096;
const int NPMT = 64;
const int NRATEMETER = 9; // ?????????????
const int NTRIGSCINT=6;
const int N_VA=NCHAN/VA_CHAN;
const int SIDE_VA=6;
const int LADDER_VA=SIDE_VA*2;

const int N_SIDES=2;
const int N_LADDER=6;
const int LADDER_BIN=NCHAN/N_LADDER;
const double PITCH=182.e-6;

#endif
