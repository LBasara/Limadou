#ifndef RootTrackerHit_h
#define RootTrackerHit_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include <vector>

#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

////////////////////////////////////////////////////////////////////////////////
//
class RootTrackerHit : public TObject
{
public:
  RootTrackerHit ();
  RootTrackerHit (TVector3 aEntry,TVector3 aExit,float aMom,float aToF,float aEloss,
		  int aPDG ,unsigned int aID,unsigned int atkID,float aTheta,float aPhi);
  
  ~RootTrackerHit ();
  RootTrackerHit (const RootTrackerHit&);
  
  inline TVector3 GetEntryPoint(){return theEntryPoint;}
  inline TVector3 GetExitPoint(){return theExitPoint;}
  inline float GetMomentum(){return theMomentum;}
  inline float GetToF(){return theTimeOfFlight;}
  inline float GetELoss(){return theEnergyLoss;}
  inline int GetParticleType(){return theParticeType;}
  inline unsigned int GetDetectorId(){return theDetectorId;}
  inline unsigned int GetTrackId(){return theTrackIdInHit;}
  inline float GetThetaAtEntry(){return theThetaAtEntry;}
  inline float GetPhiAtEntry(){return thePhiAtEntry;}

  ClassDef(RootTrackerHit,1);
  
private:
  TVector3 theEntryPoint;
  TVector3 theExitPoint;
  float theMomentum;
  float theTimeOfFlight;
  float theEnergyLoss;
  int    theParticeType;
  unsigned int theDetectorId;
  unsigned int theTrackIdInHit;
  float theThetaAtEntry;
  float thePhiAtEntry;
};

#endif
