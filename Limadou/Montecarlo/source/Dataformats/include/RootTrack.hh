#ifndef RootTrack_h
#define RootTrack_h 1
////////////////////////////////////////////////////////////////////////////////
//

#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

////////////////////////////////////////////////////////////////////////////////
//
class RootTrack : public TObject
{
public:
  RootTrack ();
  RootTrack (int aTkId,int aPDG,TVector3 aPos,TVector3 aDir,float aKinE,int aPTkId);
  
  ~RootTrack ();
  RootTrack (const RootTrack&);
  
  
  inline int      GetTrackID(){return theTrackId;}
  inline int      GetPDG(){return thePDG;}
  inline TVector3 GetDirection(){return theDirection;}
  inline TVector3 GetPosition(){return thePosition;}
  inline float    GetKinEnergy(){return theKinEnergy;}
  inline int      GetParentTrackID(){return theParentTrackId;}

  ClassDef(RootTrack,2);
  
private:
  int          theTrackId;	 	 
  int          thePDG;	 		 
  TVector3     thePosition;	 
  TVector3     theDirection;		 
  float        theKinEnergy;
  int          theParentTrackId;	 
};

#endif
