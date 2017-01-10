#include "HEPDSWProducerManager.hh"
#include "HEPDSWProducerMessenger.hh"

//Analyser
#include "G4RunManager.hh"
#include "G4ParticleTable.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4PrimaryParticle.hh"

#include "G4Step.hh"
#include "G4Timer.hh"
#include "Randomize.hh"

#include "RootCaloHit.hh"
#include "CalorimeterSD.hh"
#include "VetoSD.hh"
#include <map>
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VTrajectory.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4TrajectoryContainer.hh"
#include <iterator>

HEPDSWProducerManager* HEPDSWProducerManager::instance = 0;

////////////////////////////////////////////////////////////////////////////////
//
HEPDSWProducerManager::HEPDSWProducerManager():theEvent(0),theRootFile(0),theEventTree(0),thePathDir(0)
{
  theMessenger = new HEPDSWProducerMessenger(this);
  vetoHitCollID=-1;	
  trackerHitCollID=-1;	
  trackCollID=-1;	
  vertexCollID=-1;	
  theAutoSaveLimit = 1000;
  eventID=0;
  verboseLevel = 0;
  saveTracker=true;
  saveMCTruth=true;
  saveCalo=true;
  theEvent = new RootEvent();
}
////////////////////////////////////////////////////////////////////////////////
//  
HEPDSWProducerManager::~HEPDSWProducerManager() 
{
  delete theMessenger;
  delete theEventTree;
  delete theEvent;
  delete theRootFile;
}
////////////////////////////////////////////////////////////////////////////////
//

void HEPDSWProducerManager::SetRootFile(G4String aFileName,G4String aDirName) 
{
  theRootFile = new TFile(aFileName,"RECREATE","Storing of HEPD info");
  
  G4String aNameDir = aDirName;
  thePathDir = theRootFile->mkdir(aNameDir);

  thePathDir->cd();

  theEventTree = new TTree("EventTree","The Tree with the variable used to performe the calculation of energy deposition on the HEPD detector");
  theEventTree->Branch("Event","RootEvent",&theEvent);
  theEventTree->SetAutoSave(theAutoSaveLimit);
}


HEPDSWProducerManager* HEPDSWProducerManager::GetInstance()
{
  if (instance == 0) instance = new HEPDSWProducerManager;
  return instance;
}
////////////////////////////////////////////////////////////////////////////////
//
void HEPDSWProducerManager::BeginOfEventAction(const G4Event*)
{
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  if(caloHitCollID<0||vetoHitCollID<0||trackerHitCollID<0||trackCollID<0){
    if(saveCalo){
      caloHitCollID = SDman->GetCollectionID("caloCollection");
      vetoHitCollID = SDman->GetCollectionID("vetoCollection");
    }
    if(saveTracker)
      trackerHitCollID = SDman->GetCollectionID("trackerHitCollection");
    if(saveMCTruth){
      trackCollID = SDman->GetCollectionID("trackCollection");
      vertexCollID = SDman->GetCollectionID("vertexCollection");
    }
  } 
}
  


void HEPDSWProducerManager::StoreRootFile()
{
  theRootFile->Write(0,TObject::kOverwrite);
  theRootFile->Close();
}
void HEPDSWProducerManager::BeginOfRunAction(const G4Run*)
{
  eventID=0;
}
///////////////////////////////////////////////////////////////////////////////
//
void HEPDSWProducerManager::EndOfEventAction(const G4Event* evt)
{ 
  if(verboseLevel>0)
    G4cout << "entering in EndOfEventAction..." << G4endl;

  eventID = evt->GetEventID();

  
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  
  CaloHitsCollection * caloHC  = 0;
  CaloHitsCollection * vetoHC  = 0;
  TrackerHitsCollection * trackerHC   = 0;
  TracksCollection * trackHC   = 0;
  VertexsCollection * vertexHC   = 0;

  if(theCaloHitContainer.size())
    theCaloHitContainer.clear();
  if(theVetoHitContainer.size())
    theVetoHitContainer.clear();
  if(theTrackerHitContainer.size())
    theTrackerHitContainer.clear();
  if(theTrackContainer.size())
    theTrackContainer.clear();
  if(theVertexContainer.size())
    theVertexContainer.clear();

  if(HCE){
    if(verboseLevel>0)
      std::cout<<"Evento # "<<eventID<<std::endl;
    if(!(trackerHitCollID<0)){
      trackerHC = (TrackerHitsCollection*)(HCE->GetHC(trackerHitCollID));
      for(int i=0;i<trackerHC->entries();i++){
	TVector3 Entry((*trackerHC)[i]->GetEntryPoint().getX(),(*trackerHC)[i]->GetEntryPoint().getY(),(*trackerHC)[i]->GetEntryPoint().getZ());
	TVector3 Exit((*trackerHC)[i]->GetExitPoint().getX(),(*trackerHC)[i]->GetExitPoint().getY(),(*trackerHC)[i]->GetExitPoint().getZ());
	theTrackerHitContainer.push_back(RootTrackerHit(Entry,Exit,
							(*trackerHC)[i]->GetMomentum(),
							(*trackerHC)[i]->GetToF(),
							(*trackerHC)[i]->GetELoss(),
							(*trackerHC)[i]->GetParticleType(),
							(*trackerHC)[i]->GetDetectorId(),
							(*trackerHC)[i]->GetTrackId(),
							(*trackerHC)[i]->GetThetaAtEntry(),
							(*trackerHC)[i]->GetPhiAtEntry()));
	if(verboseLevel>0)
	  std::cout<<"TrackerHit  # "<<i<<" ; Edep = "<<(*trackerHC)[i]->GetELoss()<<" MeV"<<std::endl;
      }
    }
    if(!(caloHitCollID<0)){
      caloHC = (CaloHitsCollection*)(HCE->GetHC(caloHitCollID));
      for(int i=0;i<caloHC->entries();i++){
	theCaloHitContainer.push_back(RootCaloHit((*caloHC)[i]->GetVolume(),(*caloHC)[i]->GetTotalEdep()/MeV,(*caloHC)[i]->GetEdepMap()));
	if(verboseLevel>0)
	  std::cout<<"CaloHit  # "<<i<<" ; Volume = "<<(*caloHC)[i]->GetVolume()<<" ; Edep = "<<(*caloHC)[i]->GetTotalEdep()/MeV<<" MeV"<<std::endl;
      }
    }
    if(!(vetoHitCollID<0)){
      vetoHC = (CaloHitsCollection*)(HCE->GetHC(vetoHitCollID));
      for(int i=0;i<vetoHC->entries();i++){
	theVetoHitContainer.push_back(RootCaloHit((*vetoHC)[i]->GetVolume(),(*vetoHC)[i]->GetTotalEdep()/MeV,(*vetoHC)[i]->GetEdepMap()));
	if(verboseLevel>0)
	  std::cout<<"VetoHit  # "<<i<<" ; Volume = "<<(*vetoHC)[i]->GetVolume()<<" ; Edep = "<<(*vetoHC)[i]->GetTotalEdep()/MeV<<" MeV"<<std::endl;
      }
    }
    
    if(!(trackCollID<0)){
      trackHC = (TracksCollection*)(HCE->GetHC(trackCollID));
      for(int i=0;i<trackHC->entries();i++){
	TVector3 Pos((*trackHC)[i]->GetPosition().getX(),(*trackHC)[i]->GetPosition().getY(),(*trackHC)[i]->GetPosition().getZ());
	TVector3 Dir((*trackHC)[i]->GetMomentumDirection().getX(),(*trackHC)[i]->GetMomentumDirection().getY(),(*trackHC)[i]->GetMomentumDirection().getZ());
	theTrackContainer.push_back(RootTrack((*trackHC)[i]->GetTrackId(),
					      (*trackHC)[i]->GetPDGCode(),
					      Pos,Dir,
					      (*trackHC)[i]->GetKinEnergy(),
					      (*trackHC)[i]->GetMotherTrackId()));
      }
    }
    if(!(vertexCollID<0)){
      vertexHC = (VertexsCollection*)(HCE->GetHC(vertexCollID));
      if(vertexHC->entries()){
	for(int i=0;i<vertexHC->entries();i++){
	  TString volumeName((std::string)(*vertexHC)[i]->GetVolumeName());
	  TVector3 Pos((*vertexHC)[i]->GetPosition().getX(),(*vertexHC)[i]->GetPosition().getY(),(*vertexHC)[i]->GetPosition().getZ());
	  theVertexContainer.push_back(RootVertex((*vertexHC)[i]->IsQuasielastic(),
						  (*vertexHC)[i]->IsInelastic(),
						  volumeName,
						  Pos));
	}
      }else{
	theVertexContainer.push_back(RootVertex(false,false,"NA",TVector3(0,0,0)));
      }
    }
    

    if(verboseLevel>0)
      std::cout<<"=============================================================================================="<<std::endl;
  }


  if(verboseLevel>0)
    G4cout << "theEvent: " << theEvent << G4endl;
  theEvent->SetEventID(eventID);
  if(verboseLevel>0)
    G4cout << "saveCalo: " << saveCalo <<  "\tsaveTracker: " << saveTracker <<"\tsaveMCTruth: " << saveMCTruth << G4endl;
  if(saveCalo){
    theEvent->SetCaloHit(theCaloHitContainer);
    theEvent->SetVetoHit(theVetoHitContainer);
  }
  if(saveTracker)
    theEvent->SetTrackerHit(theTrackerHitContainer);
  if(saveMCTruth){
    theEvent->SetTracks(theTrackContainer);
    theEvent->SetVertex(theVertexContainer);
  }
 if(verboseLevel>0)
   G4cout << "theEventTree: 0x" << theEventTree << G4endl;
  theEventTree->Fill();
  if(verboseLevel>0)
    G4cout << "leaving in EndOfEventAction..." << G4endl;

}

