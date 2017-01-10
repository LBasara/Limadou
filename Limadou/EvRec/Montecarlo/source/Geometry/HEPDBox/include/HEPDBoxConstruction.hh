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
/// \file electromagnetic/TestEm3/include/HEPDBoxConstruction.hh
/// \brief Definition of the HEPDBoxConstruction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HEPDBoxConstruction_h
#define HEPDBoxConstruction_h 1

#include "globals.hh"
#include "HEPDBoxConstructionConfig1.hh"
#include "HEPDBoxConstructionConfig2.hh"
#include "HEPDBoxConstructionConfig3.hh"
#include "HEPDBoxConstructionConfig4.hh"

class G4IntersectionSolid;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class HEPDMaterial;
class G4UniformMagField;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HEPDBoxConstruction 
{
public:
  
  HEPDBoxConstruction();
  ~HEPDBoxConstruction();
  
  void SetBlanketMaterial(G4String aMat);
  void SetBlanketMaterial(G4String aMat1,G4String aMat2,G4String aMat3,G4String aMat4);
  void SetWallMaterial(G4String aMat);
  void SetWallMaterial(G4String aMat1,G4String aMat2);

  void Builder(G4String config,G4VPhysicalVolume* motherVolume);
  
private:
  HEPDBoxConstructionConfig1* theConfig1;
  HEPDBoxConstructionConfig2* theConfig2;
  HEPDBoxConstructionConfig3* theConfig3;
  HEPDBoxConstructionConfig4* theConfig4;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif














				
				
				
