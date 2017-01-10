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
/// \file electromagnetic/TestEm3/include/HEPDBoxConstructionConfig4.hh
/// \brief Definition of the HEPDBoxConstructionConfig4 class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HEPDBoxConstructionConfig4_h
#define HEPDBoxConstructionConfig4_h 1

#include "globals.hh"
#include <vector>

class G4Box;
class G4SubtractionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class HEPDSWMaterial;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HEPDBoxConstructionConfig4
{
public:
  
  HEPDBoxConstructionConfig4();
  ~HEPDBoxConstructionConfig4();

  void SetBlanketMaterial(G4String aMat);
  void SetWallMaterial(G4String aMat1,G4String aMat2);

  void Builder(G4VPhysicalVolume* motherVolume);

private:
  void ComputeObjectsPositioning();

  G4String blanketMaterial;
  G4String wallMaterial;
  G4String wallHoneyCombMaterial;

  G4double fBlanket_X;
  G4double fBlanket_Y;
  G4double fBlanket_Z;

  G4double fWallExternal_X;
  G4double fWallExternal_Y;
  G4double fWallExternal_Z;

  G4double fWallHoneyComb_X;
  G4double fWallHoneyComb_Y;
  G4double fWallHoneyComb_Z;

  G4double fWallHole_X;
  G4double fWallHole_Y;
  G4double fWallHole_Z;

  G4double transWallHole_X;
  G4double transWallHole_Y;

  G4double fPhysiBlanket_Z;

  G4double fPhysiWallExternalIn_X;
  G4double fPhysiWallExternalIn_Y;
  G4double fPhysiWallExternalIn_Z;

  G4double fPhysiWallHoneyComb_X;
  G4double fPhysiWallHoneyComb_Y;
  G4double fPhysiWallHoneyComb_Z;

  G4double fPhysiWallExternalOut_X;
  G4double fPhysiWallExternalOut_Y;
  G4double fPhysiWallExternalOut_Z;

  HEPDSWMaterial*  pMaterial;

  G4Box* fSolidBlanket;
  G4SubtractionSolid* fSolidWallExternal;
  G4SubtractionSolid* fSolidWallHoneyComb;

  G4LogicalVolume* fLogicBlanket;
  G4LogicalVolume* fLogicWallExternal;
  G4LogicalVolume* fLogicWallHoneyComb;

  G4VPhysicalVolume* fPhysiBlanket;
  G4VPhysicalVolume* fPhysiWallExternalIn;
  G4VPhysicalVolume* fPhysiWallHoneyComb;
  G4VPhysicalVolume* fPhysiWallExternalOut;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

