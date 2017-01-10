#ifndef RootVertex_h
#define RootVertex_h 1
////////////////////////////////////////////////////////////////////////////////
//

#include "TObject.h"
#include "TString.h"
#include "TVector3.h"


////////////////////////////////////////////////////////////////////////////////
//
class RootVertex : public TObject
{
public:
  RootVertex ();
  RootVertex (bool aQuaiselastic, bool aInelastic,
	  TString aVolumeName,TVector3 aPosition);
  
  ~RootVertex ();
  RootVertex (const RootVertex&);

  inline TString         GetVolumeName(){return theVolumeName;}
  inline TVector3        GetPosition(){return thePosition;}
  inline bool            IsQuasiElastic(){return theQuasielastic;} 
  inline bool            IsInelastic(){return theInelastic;} 
  
  ClassDef(RootVertex,2);
  
private:
  TString          theVolumeName;
  TVector3         thePosition;	 
  bool             theInelastic;
  bool             theQuasielastic;
};

#endif
