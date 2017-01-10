//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file electromagnetic/TestEm3/include/HEPDSWDetectorMessenger.hh
/// \brief Definition of the HEPDSWDetectorMessenger class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HEPDSWDetectorMessenger_h
#define HEPDSWDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class HEPDSWDetectorConstruction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HEPDSWDetectorMessenger: public G4UImessenger
{
public:
  HEPDSWDetectorMessenger(HEPDSWDetectorConstruction* );
  ~HEPDSWDetectorMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);

private:
  HEPDSWDetectorConstruction* Detector;

  G4UIdirectory*             fHepdDir;
  G4UIcommand*               fWorldSizeCmd;
  
  G4UIcmdWithABool* fCaloActivateCmd;
  G4UIcmdWithABool* fTrackerActivateCmd;
  //G4UIcmdWithABool* fScintillatorActivateCmd;
  G4UIcmdWithABool* fSatelliteActivateCmd;
  G4UIcmdWithABool* fHEPDBoxActivateCmd;

  G4UIcmdWithAString* fCaloConfigCmd;
  //G4UIcmdWithAString* fScintConfigCmd;
  G4UIcmdWithAString* fTrackerConfigCmd;
  G4UIcmdWithAString* fSatelliteConfigCmd;
  G4UIcmdWithAString* fHEPDBoxConfigCmd;

  //G4UIcmdWithAString* fCaloCaloMatConfigCmd;	      
  G4UIcommand*        fCaloCalo2MatConfigCmd;	      
  G4UIcmdWithAString* fCaloVetoMatConfigCmd;	      
  G4UIcmdWithAString* fCaloPoronMatConfigCmd;	      
  G4UIcmdWithAString* fCaloCarbonFiberMatConfigCmd;   
  G4UIcmdWithAString* fCaloHoneyCombMatConfigCmd;     
  G4UIcmdWithAnInteger* fCaloCrystalLayerConfigCmd;

  //G4UIcmdWithAString* fScintScintMatConfigCmd;	      
  //G4UIcmdWithAString* fScintPoronMatConfigCmd;	      
  //G4UIcmdWithAString* fScintCarbonFiberMatConfigCmd;  
		                                      
  G4UIcmdWithAString* fTrackerSiliconMatConfigCmd;    
  G4UIcmdWithAString* fTrackerKaptonMatConfigCmd;     
  G4UIcmdWithAString* fTrackerCarbonFiberMatConfigCmd;

  G4UIcommand*        fHEPDBoxMLBlanketMatConfigCmd;    
  G4UIcmdWithAString* fHEPDBoxBlanketMatConfigCmd;    
  G4UIcmdWithAString* fHEPDBoxWallMatConfigCmd;    
  G4UIcommand*        fHEPDBoxWallTwoMatConfigCmd;    

  G4UIcmdWithAString* fSatelliteBlanketMatConfigCmd;    
  G4UIcmdWithAString* fSatelliteWallMatConfigCmd;    



};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

