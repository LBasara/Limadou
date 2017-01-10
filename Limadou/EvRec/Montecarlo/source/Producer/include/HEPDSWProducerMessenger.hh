#ifndef HEPDSWProducerMessenger_h
#define HEPDSWProducerMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class HEPDSWProducerManager;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;


class HEPDSWProducerMessenger: public G4UImessenger
{
public:
  HEPDSWProducerMessenger(HEPDSWProducerManager * theManager);
  ~HEPDSWProducerMessenger();
    
  void SetNewValue(G4UIcommand * command,G4String newValues);
    
private:
  HEPDSWProducerManager*  pProducerManager;
  G4UIdirectory*          producerDir;
   
  G4UIcmdWithoutParameter*  StoreRootFileCmd;
  
  G4UIcommand* SetRootFileCmd;
  G4UIcommand* SetAutoSaveLimitCmd;

  
  G4UIcmdWithABool* fCaloActivateCmd;
  G4UIcmdWithABool* fTrackerActivateCmd;
  //  G4UIcmdWithABool* fScintillatorActivateCmd;
  G4UIcmdWithABool* fMCTruthActivateCmd;
};

#endif

