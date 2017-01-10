#include "HEPDSWProducerMessenger.hh"
#include "HEPDSWProducerManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4UserStackingAction.hh"
#include "G4RunManager.hh"
#include "G4UIparameter.hh"


HEPDSWProducerMessenger::HEPDSWProducerMessenger
(HEPDSWProducerManager* theManager)
  :pProducerManager(theManager)
{

  G4UIparameter* param;
  //command directory 
  //-------------------
  producerDir = new G4UIdirectory("/hepd/producer/");
  producerDir->SetGuidance("Producer control");
  
  //commands
  //--------
  SetRootFileCmd =  new G4UIcommand("/hepd/producer/SetRootFile",this);
  SetRootFileCmd->SetGuidance("Save the detector information in a ROOT file");
  param = new G4UIparameter("FileName",'s',false);
  SetRootFileCmd->SetParameter(param);
  param = new G4UIparameter("DirectoryName",'s',false);
  SetRootFileCmd->SetParameter(param);
  SetRootFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  StoreRootFileCmd =  new G4UIcmdWithoutParameter("/hepd/producer/StoreRootFile",this);
  StoreRootFileCmd->SetGuidance("Save the ROOT file");
  StoreRootFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetAutoSaveLimitCmd = new G4UIcommand("/hepd/producer/SetAutoSaveLimit",this);
  SetAutoSaveLimitCmd->SetGuidance("Set the AutoSave limit and mode you can choose between Event and Memory");
  param = new G4UIparameter("Limit",'i',false);
  SetAutoSaveLimitCmd->SetParameter(param);
  param = new G4UIparameter("Mode",'s',false);
  SetAutoSaveLimitCmd->SetParameter(param);
  SetAutoSaveLimitCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fCaloActivateCmd = new G4UIcmdWithABool("/hepd/producer/StoreCalorimeterInfo",this);
  fCaloActivateCmd->SetGuidance("Enable or disable the calorimeter");
  fCaloActivateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fTrackerActivateCmd = new G4UIcmdWithABool("/hepd/producer/StoreTrackerInfo",this);
  fTrackerActivateCmd->SetGuidance("Enable or disable the tracker");
  fTrackerActivateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

//   fScintillatorActivateCmd = new G4UIcmdWithABool("/hepd/producer/StoreScintillatorInfo",this);
//   fScintillatorActivateCmd->SetGuidance("Enable or disable the scintillator");
//   fScintillatorActivateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fMCTruthActivateCmd = new G4UIcmdWithABool("/hepd/producer/StoreMCTruthInfo",this);
  fMCTruthActivateCmd->SetGuidance("Enable or disable the storing of MCTruth info");
  fMCTruthActivateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

HEPDSWProducerMessenger::~HEPDSWProducerMessenger()
{  
  delete producerDir;
  delete fCaloActivateCmd;
  delete fTrackerActivateCmd;
  //  delete fScintillatorActivateCmd;
  delete fMCTruthActivateCmd;
  delete SetAutoSaveLimitCmd;
  delete StoreRootFileCmd;
  delete SetRootFileCmd;
}

void HEPDSWProducerMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ 
  if(command == StoreRootFileCmd)
    pProducerManager->StoreRootFile();
  if(command == SetRootFileCmd){
    const char* paramString=newValues;
    G4String fileName,dirName;
    std::istringstream is((char*)paramString);
    is >>fileName>>dirName;
    pProducerManager->SetRootFile(fileName,dirName);
  }
  if(command == SetAutoSaveLimitCmd){
    const char* paramString=newValues;
    G4String autosaveMode;
    G4int number;
    std::istringstream is((char*)paramString);
    is >>number>>autosaveMode;
    if(autosaveMode=="Event")
      pProducerManager->SetAutoSaveLimit(number);
    else if(autosaveMode=="Memory")
      pProducerManager->SetAutoSaveLimit(-1*number);
    else{
      G4cout<<"You have select the wrong AutoSave mode "<<autosaveMode<<";"<<std::endl;
      G4cout<<"The possibility is Event or Memeory by default Event has been chose"<<std::endl;
      pProducerManager->SetAutoSaveLimit(number);
    }
  }
  if (command == fCaloActivateCmd)
    pProducerManager->SaveCalorimeterInfo(fCaloActivateCmd->GetNewBoolValue(newValues));
  if (command == fTrackerActivateCmd)
    pProducerManager->SaveTrackerInfo(fTrackerActivateCmd->GetNewBoolValue(newValues));
  //  if (command == fScintillatorActivateCmd)
  //   pProducerManager->SaveScintillatorInfo(fScintillatorActivateCmd->GetNewBoolValue(newValues));
  if (command == fMCTruthActivateCmd)
    pProducerManager->SaveMCTruthInfo(fMCTruthActivateCmd->GetNewBoolValue(newValues));

}

