////////////////////////////////////////////////////////////////////////////////
//
#include "RootVertex.hh"

ClassImp(RootVertex)

////////////////////////////////////////////////////////////////////////////////
//
RootVertex::RootVertex ():thePosition(0,0,0)
{
  theVolumeName="";
  theInelastic=false;
  theQuasielastic=false;
}

////////////////////////////////////////////////////////////////////////////////
//
RootVertex::RootVertex (bool aQuasielastic,bool aInelastic, TString aVolumeName,
		TVector3 aPosition)
{ 
  theVolumeName          = aVolumeName;
  thePosition            = aPosition;	 
  theQuasielastic        = aQuasielastic;
  theInelastic           = aInelastic;
}
////////////////////////////////////////////////////////////////////////////////
//
RootVertex::~RootVertex ()
{}
////////////////////////////////////////////////////////////////////////////////
//
RootVertex::RootVertex (const RootVertex& right) : TObject()
{
  theVolumeName          = right.theVolumeName;
  thePosition            = right.thePosition;	 
  theQuasielastic        = right.theQuasielastic;
  theInelastic           = right.theInelastic;
}
