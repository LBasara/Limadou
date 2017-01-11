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
/// \file electromagnetic/TestEm5/include/HEPDSWPhysicsList.hh
/// \brief Definition of the HEPDSWPhysicsList class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HEPDSWPhysicsList_h
#define HEPDSWPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "StepMax.hh"

class G4VPhysicsConstructor;
class HEPDSWPhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HEPDSWPhysicsList: public G4VModularPhysicsList
{
public:
  HEPDSWPhysicsList();
 ~HEPDSWPhysicsList();

  virtual void ConstructParticle();
        
  void AddPhysicsList(const G4String& name);
    
  virtual void ConstructProcess();    
  void AddDecay();
  void AddStepMax();       
  void AddOpticalProcesses(); 
    
  virtual void SetCuts();
  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);
  void SetOpticalProcesses(bool b) {fUseOpticalProcesses = b;};
    
private:

  HEPDSWPhysicsListMessenger* fMessenger; 
  StepMax* stepMaxProcess;
  G4String fEmName;
  G4VPhysicsConstructor*  fEmHEPDSWPhysicsList;
    
  G4double fCutForGamma;
  G4double fCutForElectron;
  G4double fCutForPositron;    
  G4bool   fUseOpticalProcesses;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

