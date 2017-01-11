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
//          Filippo Ambroglini : filippo.ambroglini@pg.infn.it
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HEPDBoxConstructionConfig1.hh"


#include "G4NistManager.hh"
#include "HEPDSWMaterial.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HEPDBoxConstructionConfig1::HEPDBoxConstructionConfig1()
  :fSolidBlanket(0),fSolidWallTemp(0),fSolidWall(0),fLogicBlanket(0),fLogicWall(0),fPhysiBlanket(0),fPhysiWall(0)
{
  pMaterial     = new HEPDSWMaterial();
  fBlanket_X  = 256.980*mm;
  fBlanket_Y  = 277.135*mm;
  fBlanket_Z  = 0.1*mm;

  fWall_X = 350*mm;  
  fWall_Y = 350*mm;  
  fWall_Z = 6*mm;  

  fWallHole_X = 236.980*mm;
  fWallHole_Y = 257.135*mm;
  fWallHole_Z = 7*mm;


  ComputeObjectsPositioning();
  
  // materials
  blanketMaterial  ="mylar";
  wallMaterial     ="Aluminium";

  /*
  poronMaterial       ="Galactic";
  carbonFiberMaterial ="Galactic";
  */
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HEPDBoxConstructionConfig1::~HEPDBoxConstructionConfig1()
{
  if (pMaterial) delete pMaterial;
}
    
void HEPDBoxConstructionConfig1::ComputeObjectsPositioning(){

  fBlanketPos_X = 0*mm;
  fBlanketPos_Y = 0*mm;
  fBlanketPos_Z = 147.8*mm+74.94*mm+fBlanket_Z/2.;
  
  fWallPos_X = 0*mm;
  fWallPos_Y = 0*mm;
  fWallPos_Z = 147.8*mm+75.04*mm+fWall_Z/2.;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
void HEPDBoxConstructionConfig1::Builder(G4VPhysicalVolume* motherVolume)
{

  pMaterial->DefineMaterials();
  G4Material* blanketMat = pMaterial->GetMaterial(blanketMaterial);
  G4Material* wallMat = pMaterial->GetMaterial(wallMaterial);
  
  fSolidBlanket = new G4Box("Blanket",fBlanket_X/2.,fBlanket_Y/2.,fBlanket_Z/2);

  fSolidWallTemp = new G4Box("Wall",fWall_X/2.,fWall_Y/2.,fWall_Z/2.);

  fSolidWall = new G4SubtractionSolid("Wall",fSolidWallTemp,new G4Box("temp",fWallHole_X/2.,fWallHole_Y/2.,fWallHole_Z));

  fLogicBlanket = new G4LogicalVolume(fSolidBlanket,
				      blanketMat,
				      "Blanket");

  fLogicWall = new G4LogicalVolume(fSolidWall,
				      wallMat,
				      "Wall");
  

  fPhysiBlanket = new G4PVPlacement(0,
				    G4ThreeVector(fBlanketPos_X,fBlanketPos_Y,fBlanketPos_Z),
				    "Blanket",
				    fLogicBlanket,                
				    motherVolume,                
				    false,                       
				    0,true);       

  fPhysiWall = new G4PVPlacement(0,
				    G4ThreeVector(fWallPos_X,fWallPos_Y,fWallPos_Z),
				    "Wall",
				    fLogicWall,                
				    motherVolume,                
				    false,                       
				    0,true);       

  //Visualization Attribute

  G4VisAttributes* attYellow = new G4VisAttributes(G4Colour::Yellow());
  attYellow->SetVisibility(true);
  attYellow->SetForceAuxEdgeVisible(true);
  fLogicBlanket->SetVisAttributes(attYellow);

  G4VisAttributes* attGray = new G4VisAttributes(G4Colour::Gray());
  attGray->SetVisibility(true);
  attGray->SetForceAuxEdgeVisible(true);
  fLogicWall->SetVisAttributes(attGray);

}


