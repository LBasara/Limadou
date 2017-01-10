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
/// \file electromagnetic/TestEm3/include/TrackerConstructionConfig3.hh
/// \brief Definition of the TrackerConstructionConfig3 class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef TrackerConstructionConfig3_h
#define TrackerConstructionConfig3_h 1

#include "globals.hh"
#include <vector>

class G4Box;
class G4SubtractionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class HEPDSWMaterial;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackerConstructionConfig3
{
public:
  
   TrackerConstructionConfig3();
  ~TrackerConstructionConfig3();

  inline void SetSiliconMaterial(G4String aMat){siliconMaterial=aMat;}
  inline void SetKaptonMaterial(G4String aMat){kaptonMaterial=aMat;}
  inline void SetCarbonFiberMaterial(G4String aMat){cfiberMaterial=aMat;}

  void Builder(G4VPhysicalVolume* motherVolume);
  
private:

  void ComputeObjectsPositioning();

  G4String cfiberMaterial;
  G4String kaptonMaterial;
  G4String siliconMaterial;

  G4int fSiSensNumber;

  G4int fLadderNumber;

  G4double fSiSens_X;
  G4double fSiSens_Y;
  G4double fSiSens_Z;

  G4double fActiveLadder_X;
  G4double fActiveLadder_Y;
  G4double fActiveLadder_Z;
  
  G4double fKaptonNLadder_X;
  G4double fKaptonNLadder_Y;
  G4double fKaptonNLadder_Z;

  G4double fKaptonPLadder_X;
  G4double fKaptonPLadder_Y;
  G4double fKaptonPLadder_Z;

  G4double fLadderPlane_X;
  G4double fLadderPlane_Y;
  G4double fLadderPlane_Z;

  G4double fLadder_X;
  G4double fLadder_Y;
  G4double fLadder_Z;

  G4double fSupport_X;
  G4double fSupport_Y;
  G4double fSupport_Z;

  G4double fSupportHole_X;
  G4double fSupportHole_Y;
  G4double fSupportHole_Z;

  G4double fL_X;
  G4double fL_Y;
  G4double fL_Z;

  G4double fSupportHoleStep_Y;
  G4double fSupportHole_Ymin;

  G4double fPhysiL1_X;
  G4double fPhysiL1_Y;
  G4double fPhysiL1_Z;
  
  G4double fPhysiL1SupportU_X;
  G4double fPhysiL1SupportU_Y;
  G4double fPhysiL1SupportU_Z;

  G4double fPhysiL1SupportD_X;
  G4double fPhysiL1SupportD_Y;
  G4double fPhysiL1SupportD_Z;

  G4double fPhysiLadderPlaneL1_X;
  G4double fPhysiLadderPlaneL1_Y;
  G4double fPhysiLadderPlaneL1_Z;

  G4double fPhysiKaptonNLadderL1_X;
  G4double fPhysiKaptonNLadderL1_Y;
  G4double fPhysiKaptonNLadderL1_Z;

  G4double fPhysiKaptonPLadderL1_X;
  G4double fPhysiKaptonPLadderL1_Y;
  G4double fPhysiKaptonPLadderL1_Z;

  G4double fPhysiActiveLadderL1_X;
  G4double fPhysiActiveLadderL1_Y;
  G4double fPhysiActiveLadderL1_Z;

  G4double fPhysiL2_X;
  G4double fPhysiL2_Y;
  G4double fPhysiL2_Z;

  G4double fPhysiL2SupportU_X;
  G4double fPhysiL2SupportU_Y;
  G4double fPhysiL2SupportU_Z;

  G4double fPhysiL2SupportD_X;
  G4double fPhysiL2SupportD_Y;
  G4double fPhysiL2SupportD_Z;

  G4double fPhysiLadderPlaneL2_X;
  G4double fPhysiLadderPlaneL2_Y;
  G4double fPhysiLadderPlaneL2_Z;

  G4double fPhysiKaptonNLadderL2_X;
  G4double fPhysiKaptonNLadderL2_Y;
  G4double fPhysiKaptonNLadderL2_Z;

  G4double fPhysiKaptonPLadderL2_X;
  G4double fPhysiKaptonPLadderL2_Y;
  G4double fPhysiKaptonPLadderL2_Z;

  G4double fPhysiActiveLadderL2_X;
  G4double fPhysiActiveLadderL2_Y;
  G4double fPhysiActiveLadderL2_Z;

  HEPDSWMaterial* pMaterial;

  G4Box* fSolidSiSensor;
  G4Box* fSolidActiveLadder;
  G4Box* fSolidKaptonNLadder;
  G4Box* fSolidKaptonPLadder;
  G4Box* fSolidLadderPlane;
  G4Box* fSolidLadder;
  G4Box* fSolidSupport;
  G4Box* fSolidSupportHole;
  std::vector<G4SubtractionSolid*> fSolidSupportT;
  G4SubtractionSolid* fSolidSupportFinal;
  G4Box* fSolid;

  G4LogicalVolume* fLogicSiSensor;

  G4LogicalVolume* fLogicActiveLadderL1;
  G4LogicalVolume* fLogicKaptonNLadderL1;
  G4LogicalVolume* fLogicKaptonPLadderL1;
  G4LogicalVolume* fLogicLadderPlaneL1;
  G4LogicalVolume* fLogicLadderL1;
  G4LogicalVolume* fLogicSupportL1;
  G4LogicalVolume* fLogicL1;

  G4LogicalVolume* fLogicActiveLadderL2;
  G4LogicalVolume* fLogicKaptonNLadderL2;
  G4LogicalVolume* fLogicKaptonPLadderL2;
  G4LogicalVolume* fLogicLadderPlaneL2;
  G4LogicalVolume* fLogicLadderL2;
  G4LogicalVolume* fLogicSupportL2;
  G4LogicalVolume* fLogicL2
;
  G4VPhysicalVolume* fPhysiL1;
  G4VPhysicalVolume* fPhysiL1SupportU;
  G4VPhysicalVolume* fPhysiL1SupportD;
  G4VPhysicalVolume* fPhysiLadderPlaneL1;
  G4VPhysicalVolume* fPhysiLadderL1;
  G4VPhysicalVolume* fPhysiKaptonNLadderL1;
  G4VPhysicalVolume* fPhysiKaptonPLadderL1;
  G4VPhysicalVolume* fPhysiActiveLadderL1;
  G4VPhysicalVolume* fPhysiSiSensorL1;

  G4VPhysicalVolume* fPhysiL2;
  G4VPhysicalVolume* fPhysiL2SupportU;
  G4VPhysicalVolume* fPhysiL2SupportD;
  G4VPhysicalVolume* fPhysiLadderPlaneL2;
  G4VPhysicalVolume* fPhysiLadderL2;
  G4VPhysicalVolume* fPhysiKaptonNLadderL2;
  G4VPhysicalVolume* fPhysiKaptonPLadderL2;
  G4VPhysicalVolume* fPhysiActiveLadderL2;
  G4VPhysicalVolume* fPhysiSiSensorL2;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

