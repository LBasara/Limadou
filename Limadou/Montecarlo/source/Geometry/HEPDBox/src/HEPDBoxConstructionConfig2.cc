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

#include "HEPDBoxConstructionConfig2.hh"


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

HEPDBoxConstructionConfig2::HEPDBoxConstructionConfig2()
  :fSolidBlanket(0),fLogicBlanket(0),fPhysiBlanket(0)
{
  pMaterial     = new HEPDSWMaterial();
  fBlanket_X  = 277.135*mm;
  fBlanket_Y  = 256.980*mm;
  fBlanket_Z  = 0.1*mm;

  ComputeObjectsPositioning();
  
  // materials
  blanketMaterial  ="mylar";
  /*
  poronMaterial       ="Galactic";
  carbonFiberMaterial ="Galactic";
  */
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HEPDBoxConstructionConfig2::~HEPDBoxConstructionConfig2()
{
  if (pMaterial) delete pMaterial;
}
    
void HEPDBoxConstructionConfig2::ComputeObjectsPositioning(){

  fBlanketPos_X = 0*mm;
  fBlanketPos_Y = 0*mm;
  fBlanketPos_Z = 206.98*mm;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
void HEPDBoxConstructionConfig2::Builder(G4VPhysicalVolume* motherVolume)
{

  pMaterial->DefineMaterials();
  G4Material* blanketMat = pMaterial->GetMaterial(blanketMaterial);
  
  fSolidBlanket = new G4Box("Blanket",fBlanket_X/2.,fBlanket_Y/2.,fBlanket_Z/2);

  fLogicBlanket = new G4LogicalVolume(fSolidBlanket,
				      blanketMat,
				      "Blanket");


  fPhysiBlanket = new G4PVPlacement(0,
				    G4ThreeVector(fBlanketPos_X,fBlanketPos_Y,fBlanketPos_Z),
				    "Blanket",
				    fLogicBlanket,                
				    motherVolume,                
				    false,                       
				    0,true);       

  //Visualization Attribute

  G4VisAttributes* attYellow = new G4VisAttributes(G4Colour::Yellow());
  attYellow->SetVisibility(true);
  attYellow->SetForceAuxEdgeVisible(true);
  fLogicBlanket->SetVisAttributes(attYellow);

}


