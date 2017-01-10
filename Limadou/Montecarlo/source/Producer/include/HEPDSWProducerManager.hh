#ifndef HEPDSWProducerManager_HH
#define HEPDSWProducerManager_HH

#include"G4ios.hh"
#include"G4strstreambuf.hh"
#include <vector>
#include"globals.hh"
#include <fstream>
#include"G4ThreeVector.hh"
#include"G4Event.hh"
#include"G4Run.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TDirectory.h"
#include "TVector3.h"

#include "RootEvent.hh"
#include "RootCaloHit.hh"
#include "RootTrackerHit.hh"
#include "RootTrack.hh"
#include "RootVertex.hh"
#include "CaloHit.hh"
#include "TrackerHit.hh"
#include "Track.hh"
#include "Vertex.hh"
#include <string>

class G4Step;
class HEPDSWProducerMessenger;
class HEPDSWProducerManager
{
public:
  
  ~HEPDSWProducerManager();
  static HEPDSWProducerManager* GetInstance();
   
  //Event actions
  //-------------
  void BeginOfRunAction(const G4Run*);
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
    
  void SetRootFile(G4String aFileName,G4String aDirName);
  void StoreRootFile();
  
  inline void SaveTrackerInfo(G4bool aVal){saveTracker=aVal;}
  inline void SaveMCTruthInfo(G4bool aVal){saveMCTruth=aVal;}
  inline void SaveCalorimeterInfo(G4bool aVal){saveCalo=aVal;}

  inline void SetAutoSaveLimit(G4int aAutoSave){theAutoSaveLimit=aAutoSave;}

private:
  static HEPDSWProducerManager* instance;

private:
  HEPDSWProducerManager(); 
 
private:
  HEPDSWProducerMessenger* theMessenger;
  
  G4int verboseLevel;
  
  G4bool saveTracker;
  G4bool saveMCTruth;
  G4bool saveCalo;

  G4int caloHitCollID;
  G4int vetoHitCollID;
  G4int trackerHitCollID;
  G4int trackCollID;
  G4int vertexCollID;

  G4int theAutoSaveLimit;

  G4int eventID;

  std::vector<RootCaloHit> theCaloHitContainer;
  std::vector<RootCaloHit> theVetoHitContainer;
  std::vector<RootTrackerHit> theTrackerHitContainer;
  std::vector<RootTrack> theTrackContainer;
  std::vector<RootVertex> theVertexContainer;

  RootEvent* theEvent;
  
  TFile* theRootFile;
  TTree* theEventTree;
  TDirectory* thePathDir;
  
};

#endif




