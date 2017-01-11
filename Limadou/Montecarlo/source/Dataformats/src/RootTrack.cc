////////////////////////////////////////////////////////////////////////////////
//
#include "RootTrack.hh"

ClassImp(RootTrack)

////////////////////////////////////////////////////////////////////////////////
//
RootTrack::RootTrack():theTrackId(0),thePDG(0),thePosition(0,0,0),theDirection(0,0,0),theKinEnergy(0),theParentTrackId(0)
{
  ;
}
////////////////////////////////////////////////////////////////////////////////
//
RootTrack::RootTrack (int aTkId,int aPDG,TVector3 aPos,TVector3 aDir,float aKinE,int aPTkId)
{   
  theTrackId = aTkId;	 	 
  thePDG = aPDG;	 		 
  thePosition = aPos;	 
  theDirection = aDir;	 
  theKinEnergy = aKinE;		 
  theParentTrackId = aPTkId;	 
}
////////////////////////////////////////////////////////////////////////////////
//
RootTrack::~RootTrack ()
{}
////////////////////////////////////////////////////////////////////////////////
//
RootTrack::RootTrack (const RootTrack& right) : TObject()
{
  theTrackId = right.theTrackId;
  thePDG = right.thePDG;
  thePosition = right.thePosition;	 
  theDirection = right.theDirection;	 
  theKinEnergy = right.theKinEnergy;		 
  theParentTrackId = right.theParentTrackId;
}
////////////////////////////////////////////////////////////////////////////////
//
