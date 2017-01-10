////////////////////////////////////////////////////////////////////////////////
//
#include "RootCaloHit.hh"
#include <iostream>
ClassImp(RootCaloHit)

////////////////////////////////////////////////////////////////////////////////
//
RootCaloHit::RootCaloHit():theTotalEdep(0),theParticleList(0)
{
  //  theVolume="";
  theDet=-1000;
}
////////////////////////////////////////////////////////////////////////////////
//
//RootCaloHit::RootCaloHit(TString aVolume,double aTotalEdep,std::map<int,double> aEdep) 
RootCaloHit::RootCaloHit(int aDetID, double aTotalEdep,std::map<int,double> aEdep) 
{ 
  //  theVolume     = aVolume;
  theDet     = aDetID;
  theTotalEdep  = aTotalEdep;	  
  theEdep       = aEdep;
  for(std::map<int,double>::iterator i = theEdep.begin();i!=theEdep.end();i++){
    theParticleList.push_back(i->first);
  }
}
////////////////////////////////////////////////////////////////////////////////
//
RootCaloHit::~RootCaloHit()
{;}
////////////////////////////////////////////////////////////////////////////////
//
RootCaloHit::RootCaloHit(const RootCaloHit& right) : TObject()
{
  //  theVolume       = right.theVolume;
  theDet       = right.theDet;
  theTotalEdep    = right.theTotalEdep;	 
  theEdep         = right.theEdep;
  theParticleList = right.theParticleList;
}
////////////////////////////////////////////////////////////////////////////////
//
double RootCaloHit::GetEdep(int aTkID)
{
  if(theEdep.find(aTkID)==theEdep.end()){
    std::cout<<"You are asking the wrong track ID = "<<aTkID<<" try with one of these ";
      for(unsigned int i=0;i<theParticleList.size();i++)
	std::cout<<theParticleList[i]<<"; ";
    std::cout<<" "<<std::endl;
    return 0;
  }else{
    return theEdep[aTkID];
  }
}
