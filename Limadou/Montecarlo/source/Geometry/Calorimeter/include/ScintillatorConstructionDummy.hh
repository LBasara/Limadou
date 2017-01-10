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
/// \file electromagnetic/TestEm3/include/ScintillatorConstructionDummy.hh
/// \brief Definition of the ScintillatorConstructionDummy class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ScintillatorConstructionDummy_h
#define ScintillatorConstructionDummy_h 1

#include "globals.hh"
#include <vector>

class G4Box;
class G4SubtractionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class HEPDSWMaterial;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ScintillatorConstructionDummy
{
public:
  
  ScintillatorConstructionDummy();
  ~ScintillatorConstructionDummy();

  inline void SetScintillatorMaterial(G4String aMat){scintMaterial=aMat;}
  void Builder(G4VPhysicalVolume* motherVolume);

private:
  void ComputeObjectsPositioning();

  G4String scintMaterial;
  G4double fS1_X;
  G4double fS1_Y;
  G4double fS1_Z;

  G4double fS2_X;
  G4double fS2_Y;
  G4double fS2_Z;

  G4double fPhysiS1Scintillator_X;
  G4double fPhysiS1Scintillator_Y;
  G4double fPhysiS1Scintillator_Z;

  G4double fPhysiS2Scintillator_X;
  G4double fPhysiS2Scintillator_Y;
  G4double fPhysiS2Scintillator_Z;

  HEPDSWMaterial*        pMaterial;


  G4Box* fSolidS2Scint;
  G4Box* fSolidS1Scint;

  G4LogicalVolume* fLogicS2Scint;
  G4LogicalVolume* fLogicS1Scint;

  G4VPhysicalVolume* fPhysiS2Scint;
  G4VPhysicalVolume* fPhysiS1Scint;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

