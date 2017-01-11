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
/// \file electromagnetic/TestEm3/include/HEPDBoxConstructionConfig1.hh
/// \brief Definition of the HEPDBoxConstructionConfig1 class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HEPDBoxConstructionConfig1_h
#define HEPDBoxConstructionConfig1_h 1

#include "globals.hh"
#include <vector>

class G4Box;
class G4SubtractionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class HEPDSWMaterial;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HEPDBoxConstructionConfig1
{
public:
  
  HEPDBoxConstructionConfig1();
  ~HEPDBoxConstructionConfig1();

  inline void SetBlanketMaterial(G4String aMat){blanketMaterial=aMat;}
  inline void SetWallMaterial(G4String aMat){wallMaterial=aMat;}

  void Builder(G4VPhysicalVolume* motherVolume);

private:
  void ComputeObjectsPositioning();

  G4String blanketMaterial;
  G4String wallMaterial;

  G4double fBlanket_X;
  G4double fBlanket_Y;
  G4double fBlanket_Z;
           
  G4double fWall_X;  
  G4double fWall_Y;  
  G4double fWall_Z;  

  G4double fWallHole_X;  
  G4double fWallHole_Y;  
  G4double fWallHole_Z;  
  
  G4double fBlanketPos_X;
  G4double fBlanketPos_Y;
  G4double fBlanketPos_Z;
           
  G4double fWallPos_X;
  G4double fWallPos_Y;
  G4double fWallPos_Z;


  HEPDSWMaterial*        pMaterial;

  G4Box* fSolidBlanket;
  G4Box* fSolidWallTemp;
  G4SubtractionSolid* fSolidWall;

  G4LogicalVolume* fLogicBlanket;
  G4LogicalVolume* fLogicWall;

  G4VPhysicalVolume* fPhysiBlanket;
  G4VPhysicalVolume* fPhysiWall;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

