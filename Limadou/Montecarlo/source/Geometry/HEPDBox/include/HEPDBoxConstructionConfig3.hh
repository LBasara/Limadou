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
/// \file electromagnetic/TestEm3/include/HEPDBoxConstructionConfig3.hh
/// \brief Definition of the HEPDBoxConstructionConfig3 class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HEPDBoxConstructionConfig3_h
#define HEPDBoxConstructionConfig3_h 1

#include "globals.hh"
#include <vector>

class G4Box;
class G4UnionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class HEPDSWMaterial;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HEPDBoxConstructionConfig3
{
public:
  
  HEPDBoxConstructionConfig3();
  ~HEPDBoxConstructionConfig3();

  void SetBlanketMaterial(G4String asuteMat,G4String asuteCoatMat,G4String asoftMat,G4String asupportMat);

  void Builder(G4VPhysicalVolume* motherVolume);

private:
  void ComputeObjectsPositioning();

  G4String blanketMaterial;
  G4String suteMaterial;
  G4String suteCoatMaterial;
  G4String softMaterial;
  G4String supportMaterial;

  G4double fBlanket_X;
  G4double fBlanket_Y;
  G4double fBlanket_Z;
             

  G4double fSuteCoat_X; 
  G4double fSuteCoat_Y; 
  G4double fSuteCoat_Z; 
             
  G4double fSuteThick_X;
  G4double fSuteThick_Y;
  G4double fSuteThick_Z;
             
  G4double fSuteThin_X; 
  G4double fSuteThin_Y; 
  G4double fSuteThin_Z; 
             
  G4double fSoftMat_X;  
  G4double fSoftMat_Y;  
  G4double fSoftMat_Z;  
             
  G4double fSupportO_X; 
  G4double fSupportO_Y; 
  G4double fSupportO_Z; 
             
  G4double fSupportV_X; 
  G4double fSupportV_Y; 
  G4double fSupportV_Z; 


  G4double fBlanketPos_Z;    
  G4double fSuteCoat1Pos_Z;  
  G4double fSuteThick1Pos_Z; 
  G4double fSoftMat1Pos_Z;  
  G4double fSuteCoat2Pos_Z;  
  G4double fSuteThin1Pos_Z;  
  G4double fSuteCoat3Pos_Z;  
  G4double fSoftMat2Pos_Z;   
  G4double fSuteCoat4Pos_Z;  
  G4double fSuteThin2Pos_Z;  
  G4double fSuteCoat5Pos_Z;  
  G4double fSoftMat3Pos_Z;   
  G4double fSuteCoat6Pos_Z;  
  G4double fSuteThin3Pos_Z;  
  G4double fSuteCoat7Pos_Z;  
  G4double fSoftMat4Pos_Z;   
  G4double fSuteCoat8Pos_Z;  
  G4double fSuteThick2Pos_Z; 
  G4double fSuteCoat9Pos_Z;  
  G4double fSupportPos_Z;    
             
  HEPDSWMaterial*  pMaterial;

  G4Box* fSolidBlanket;
  G4Box* fSolidSuteCoat;
  G4Box* fSolidSuteThin;
  G4Box* fSolidSuteThick;
  G4Box* fSolidSoftMat;
  G4Box* fSolidSupportO;
  G4Box* fSolidSupportV;
  G4UnionSolid* fSolidSupport;

  G4LogicalVolume* fLogicBlanket;
  G4LogicalVolume* fLogicSuteCoat;  
  G4LogicalVolume* fLogicSuteThin;
  G4LogicalVolume* fLogicSuteThick;
  G4LogicalVolume* fLogicSoftMat;
  G4LogicalVolume* fLogicSupport;

  G4VPhysicalVolume* fPhysiBlanket;
  G4VPhysicalVolume* fPhysiSuteCoat1;
  G4VPhysicalVolume* fPhysiSuteThick1;
  G4VPhysicalVolume* fPhysiSoftMat1;
  G4VPhysicalVolume* fPhysiSuteCoat2;
  G4VPhysicalVolume* fPhysiSuteThin1;
  G4VPhysicalVolume* fPhysiSuteCoat3;
  G4VPhysicalVolume* fPhysiSoftMat2;
  G4VPhysicalVolume* fPhysiSuteCoat4;
  G4VPhysicalVolume* fPhysiSuteThin2;
  G4VPhysicalVolume* fPhysiSuteCoat5;
  G4VPhysicalVolume* fPhysiSoftMat3;
  G4VPhysicalVolume* fPhysiSuteCoat6;
  G4VPhysicalVolume* fPhysiSuteThin3;
  G4VPhysicalVolume* fPhysiSuteCoat7;
  G4VPhysicalVolume* fPhysiSoftMat4;
  G4VPhysicalVolume* fPhysiSuteCoat8;
  G4VPhysicalVolume* fPhysiSuteThick2;
  G4VPhysicalVolume* fPhysiSuteCoat9;
  G4VPhysicalVolume* fPhysiSupport;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

