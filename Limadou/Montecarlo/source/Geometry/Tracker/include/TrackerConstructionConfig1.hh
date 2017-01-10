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
/// \file electromagnetic/TestEm3/include/TrackerConstructionConfig1.hh
/// \brief Definition of the TrackerConstructionConfig1 class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef TrackerConstructionConfig1_h
#define TrackerConstructionConfig1_h 1

#include "globals.hh"
#include <vector>

class G4Box;
class G4SubtractionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class HEPDSWMaterial;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackerConstructionConfig1
{
public:
  
   TrackerConstructionConfig1();
  ~TrackerConstructionConfig1();

  inline void SetSiliconMaterial(G4String aMat){siliconMaterial=aMat;}
  inline void SetKaptonMaterial(G4String aMat){kaptonMaterial=aMat;}
  inline void SetCarbonFiberMaterial(G4String aMat){cfiberMaterial=aMat;}

  void Builder(G4VPhysicalVolume* motherVolume);
  
private:

  void ComputeObjectsPositioning();

  G4String cfiberMaterial;
  G4String kaptonMaterial;
  G4String siliconMaterial;

  G4int fSiSensNumberL1;
  G4int fSiSensNumberL2;

  G4int fLadderNumberL1;
  G4int fLadderNumberL2;

  G4double fSiSens_X;
  G4double fSiSens_Y;
  G4double fSiSens_Z;

  G4double fActiveLadderL1_X;
  G4double fActiveLadderL1_Y;
  G4double fActiveLadderL1_Z;
  
  G4double fKaptonNLadderL1_X;
  G4double fKaptonNLadderL1_Y;
  G4double fKaptonNLadderL1_Z;

  G4double fKaptonPLadderL1_X;
  G4double fKaptonPLadderL1_Y;
  G4double fKaptonPLadderL1_Z;

  G4double fLadderPlaneL1_X;
  G4double fLadderPlaneL1_Y;
  G4double fLadderPlaneL1_Z;

  G4double fLadderL1_X;
  G4double fLadderL1_Y;
  G4double fLadderL1_Z;

  G4double fL1Support_X;
  G4double fL1Support_Y;
  G4double fL1Support_Z;

  G4double fL1SupportHoleU_X;
  G4double fL1SupportHoleU_Y;
  G4double fL1SupportHoleU_Z;

  G4double fL1SupportHoleD_X;
  G4double fL1SupportHoleD_Y;
  G4double fL1SupportHoleD_Z;

  G4double fL1_X;
  G4double fL1_Y;
  G4double fL1_Z;

  G4double fActiveLadderL2_X;
  G4double fActiveLadderL2_Y;
  G4double fActiveLadderL2_Z;

  G4double fKaptonNLadderL2_X;
  G4double fKaptonNLadderL2_Y;
  G4double fKaptonNLadderL2_Z;

  G4double fKaptonPLadderL2_X;
  G4double fKaptonPLadderL2_Y;
  G4double fKaptonPLadderL2_Z;

  G4double fLadderPlaneL2_X;
  G4double fLadderPlaneL2_Y;
  G4double fLadderPlaneL2_Z;

  G4double fLadderL2_X;
  G4double fLadderL2_Y;
  G4double fLadderL2_Z;

  G4double fL2Support_X;
  G4double fL2Support_Y;
  G4double fL2Support_Z;

  G4double fL2SupportHoleU_X;
  G4double fL2SupportHoleU_Y;
  G4double fL2SupportHoleU_Z;

  G4double fL2SupportHoleD_X;
  G4double fL2SupportHoleD_Y;
  G4double fL2SupportHoleD_Z;

  G4double fL2_X;
  G4double fL2_Y;
  G4double fL2_Z;


  G4double fL1SupportHoleUStep_X;
  G4double fL1SupportHoleU_Xmin;

  G4double fL1SupportHoleDStep_X;
  G4double fL1SupportHoleD_Xmin;

  G4double fL1SupportHoleDStep_Y;
  G4double fL1SupportHoleD_Ymin;

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

  G4double fL2SupportHoleUStep_X;
  G4double fL2SupportHoleU_Xmin;

  G4double fL2SupportHoleDStep_X;
  G4double fL2SupportHoleD_Xmin;

  G4double fL2SupportHoleDStep_Y;
  G4double fL2SupportHoleD_Ymin;

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
  G4Box* fSolidActiveLadderL1;
  G4Box* fSolidKaptonNLadderL1;
  G4Box* fSolidKaptonPLadderL1;
  G4Box* fSolidLadderPlaneL1;
  G4Box* fSolidLadderL1;
  G4Box* fSolidL1Support;
  G4Box* fSolidL1SupportHoleU;
  G4Box* fSolidL1SupportHoleD;
  std::vector<G4SubtractionSolid*> fSolidL1SupportU;
  std::vector<G4SubtractionSolid*> fSolidL1SupportD;
  G4SubtractionSolid* fSolidL1SupportUFinal;
  G4SubtractionSolid* fSolidL1SupportDFinal;
  G4Box* fSolidL1;

  G4Box* fSolidActiveLadderL2;
  G4Box* fSolidKaptonNLadderL2;
  G4Box* fSolidKaptonPLadderL2;
  G4Box* fSolidLadderPlaneL2;
  G4Box* fSolidLadderL2;
  G4Box* fSolidL2Support;
  G4Box* fSolidL2SupportHoleU;
  G4Box* fSolidL2SupportHoleD;
  std::vector<G4SubtractionSolid*> fSolidL2SupportU;
  std::vector<G4SubtractionSolid*> fSolidL2SupportD;
  G4SubtractionSolid* fSolidL2SupportUFinal;
  G4SubtractionSolid* fSolidL2SupportDFinal;
  G4Box* fSolidL2;

  G4LogicalVolume* fLogicSiSensor;

  G4LogicalVolume* fLogicActiveLadderL1;
  G4LogicalVolume* fLogicKaptonNLadderL1;
  G4LogicalVolume* fLogicKaptonPLadderL1;
  G4LogicalVolume* fLogicLadderPlaneL1;
  G4LogicalVolume* fLogicLadderL1;
  G4LogicalVolume* fLogicL1SupportU;
  G4LogicalVolume* fLogicL1SupportD;
  G4LogicalVolume* fLogicL1;

  G4LogicalVolume* fLogicActiveLadderL2;
  G4LogicalVolume* fLogicKaptonNLadderL2;
  G4LogicalVolume* fLogicKaptonPLadderL2;
  G4LogicalVolume* fLogicLadderPlaneL2;
  G4LogicalVolume* fLogicLadderL2;
  G4LogicalVolume* fLogicL2SupportU;
  G4LogicalVolume* fLogicL2SupportD;
  G4LogicalVolume* fLogicL2;


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

