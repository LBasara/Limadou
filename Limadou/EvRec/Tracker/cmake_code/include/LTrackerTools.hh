#ifndef __LTRACKERTOOLS__
#define __LTRACKERTOOLS__ 1

int ChanToLadder(int nStrip);
int ChanToADC(int nStrip);
int ChanToVA(int nStrip);
int ChanToSide(int nStrip);  // 0 p - 1 n
int ChanToPlane(int nStrip); // 0 external - 1 internal
int ChanToLadderPlane (int nChan); // return 0,1,....11
bool SameLadderPlane(int Chan1, int Chan2);
int ChanToLadderChan(int Chan); // channel inside the ladder

#endif
