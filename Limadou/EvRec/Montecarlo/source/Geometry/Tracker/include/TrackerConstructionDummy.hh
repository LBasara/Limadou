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
/// \file electromagnetic/TestEm3/include/TrackerConstructionDummy.hh
/// \brief Definition of the TrackerConstructionDummy class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef TrackerConstructionDummy_h
#define TrackerConstructionDummy_h 1

#include "globals.hh"
#include <vector>

class G4Box;
class G4SubtractionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class HEPDSWMaterial;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackerConstructionDummy
{
public:
  
   TrackerConstructionDummy();
  ~TrackerConstructionDummy();

  inline void SetSiliconMaterial(G4String aMat){siliconMaterial=aMat;}

  void Builder(G4VPhysicalVolume* motherVolume);
  
private:

  void ComputeObjectsPositioning();

  G4String siliconMaterial;

  G4int fSiSensNumber;

  G4int fLadderNumber;

  G4double fSiSens_X;
  G4double fSiSens_Y;
  G4double fSiSens_Z;

  G4double fActiveLadder_X;
  G4double fActiveLadder_Y;
  G4double fActiveLadder_Z;
  
  G4double fLadderPlane_X;
  G4double fLadderPlane_Y;
  G4double fLadderPlane_Z;

  G4double fLadder_X;
  G4double fLadder_Y;
  G4double fLadder_Z;

  G4double fL_X;
  G4double fL_Y;
  G4double fL_Z;

  G4double fPhysiL1_X;
  G4double fPhysiL1_Y;
  G4double fPhysiL1_Z;
  
  G4double fPhysiL2_X;
  G4double fPhysiL2_Y;
  G4double fPhysiL2_Z;

  HEPDSWMaterial* pMaterial;

  G4Box* fSolidSiSensor;
  G4Box* fSolidActiveLadder;
  G4Box* fSolidLadderPlane;
  G4Box* fSolidLadder;
  G4Box* fSolid;

  G4LogicalVolume* fLogicSiSensor;

  G4LogicalVolume* fLogicActiveLadderL1;
  G4LogicalVolume* fLogicLadderPlaneL1;
  G4LogicalVolume* fLogicLadderL1;
  G4LogicalVolume* fLogicL1;

  G4LogicalVolume* fLogicActiveLadderL2;
  G4LogicalVolume* fLogicLadderPlaneL2;
  G4LogicalVolume* fLogicLadderL2;
  G4LogicalVolume* fLogicL2;

  G4VPhysicalVolume* fPhysiL1;
  G4VPhysicalVolume* fPhysiLadderPlaneL1;
  G4VPhysicalVolume* fPhysiLadderL1;
  G4VPhysicalVolume* fPhysiActiveLadderL1;
  G4VPhysicalVolume* fPhysiSiSensorL1;

  G4VPhysicalVolume* fPhysiL2;
  G4VPhysicalVolume* fPhysiLadderPlaneL2;
  G4VPhysicalVolume* fPhysiLadderL2;
  G4VPhysicalVolume* fPhysiActiveLadderL2;
  G4VPhysicalVolume* fPhysiSiSensorL2;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

