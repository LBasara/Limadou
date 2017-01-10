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

#include "ScintillatorConstructionDummy.hh"


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
#include "CalorimeterSD.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintillatorConstructionDummy::ScintillatorConstructionDummy()
  :
  fSolidS2Scint(0),
  fSolidS1Scint(0),
  fLogicS2Scint(0),
  fLogicS1Scint(0),
  fPhysiS2Scint(0),
  fPhysiS1Scint(0)
{
  pMaterial     = new HEPDSWMaterial();

  fS1_X = 160*mm;
  fS1_Y = 160*mm;
  fS1_Z = 5*mm;

  fS2_X = 150*mm;
  fS2_Y = 150*mm;
  fS2_Z = 10*mm;

  ComputeObjectsPositioning();
  
  // materials
  scintMaterial       ="Scintillator";

}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ScintillatorConstructionDummy::~ScintillatorConstructionDummy()
{
  if (pMaterial) delete pMaterial;
}
    
void ScintillatorConstructionDummy::ComputeObjectsPositioning(){

  fPhysiS1Scintillator_X = 0;
  fPhysiS1Scintillator_Y = 0;
  fPhysiS1Scintillator_Z = 100*mm+fS2_Z+fS1_Z/2.;

  fPhysiS2Scintillator_X = 0;
  fPhysiS2Scintillator_Y = 0;
  fPhysiS2Scintillator_Z = 100*mm+fS2_Z/2.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
void ScintillatorConstructionDummy::Builder(G4VPhysicalVolume* motherVolume)
{

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String scintSDname = "/hepd/calorimeter";
  CalorimeterSD * scintSD = new CalorimeterSD(scintSDname);
  SDman->AddNewDetector(scintSD);

  pMaterial->DefineMaterials();
  G4Material* scintMat = pMaterial->GetMaterial(scintMaterial);
  
  fSolidS2Scint = new G4Box("S2Scintillator",fS2_X/2.,fS2_Y/2.,fS2_Z/2.);
  
  fSolidS1Scint = new G4Box("S1Scintillator",fS1_X/2.,fS1_Y/2.,fS1_Z/2.);

  fLogicS2Scint = new G4LogicalVolume(fSolidS2Scint,scintMat,"S2Scintillator");
  
  fLogicS1Scint = new G4LogicalVolume(fSolidS1Scint,scintMat,"S1Scintillator");
  
  fLogicS1Scint->SetSensitiveDetector(scintSD);
  fLogicS2Scint->SetSensitiveDetector(scintSD);

  fPhysiS2Scint = new G4PVPlacement(0,
			       G4ThreeVector(fPhysiS2Scintillator_X,fPhysiS2Scintillator_Y,fPhysiS2Scintillator_Z),
			       "S2Scintillator",
			       fLogicS2Scint,                
			       motherVolume,                
			       false,                       
			       0,true);       

  fPhysiS1Scint = new G4PVPlacement(0,
			       G4ThreeVector(fPhysiS1Scintillator_X,fPhysiS1Scintillator_Y,fPhysiS1Scintillator_Z),
			       "S1Scintillator",
			       fLogicS1Scint,                
			       motherVolume,                
			       false,                       
			       0,true);       
  
  //Visualization Attribute
  G4VisAttributes* attMagenta = new G4VisAttributes(G4Colour::Magenta());
  attMagenta->SetVisibility(true);
  attMagenta->SetForceAuxEdgeVisible(true);
  fLogicS1Scint->SetVisAttributes(attMagenta);
  fLogicS2Scint->SetVisAttributes(attMagenta);

}


