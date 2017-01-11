#ifndef RootCaloHit_h
#define RootCaloHit_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include <vector>
#include <map>

#include "TObject.h"
#include "TString.h"

////////////////////////////////////////////////////////////////////////////////
//
class RootCaloHit : public TObject
{
public:
  RootCaloHit ();
  //  RootCaloHit (TString aVolume,double aTotalEdep,std::map<int,double> aEdep);
  RootCaloHit (int aDetID,double aTotalEdep,std::map<int,double> aEdep);
  
  ~RootCaloHit ();
  RootCaloHit (const RootCaloHit&);
  
  
  inline double      GetTotalEdep(){return theTotalEdep;}
  //  inline TString     GetVolume(){return theVolume;}
  inline int GetVolume(){return theDet; }
  inline std::vector<int> GetParticleList(){return theParticleList;}
  double      GetEdep(int aTkID);

  ClassDef(RootCaloHit,2);
  
private:
  //  TString         theVolume;  
  int    theDet;
  double   theTotalEdep;	 
  std::map<int,double> theEdep;
  std::vector<int>     theParticleList;
};

#endif
