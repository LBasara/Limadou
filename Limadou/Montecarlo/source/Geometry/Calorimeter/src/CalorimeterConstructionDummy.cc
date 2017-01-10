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

#include "CalorimeterConstructionDummy.hh"

#include "G4NistManager.hh"
#include "HEPDSWMaterial.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4IntersectionSolid.hh"


#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "CalorimeterSD.hh"
#include "VetoSD.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterConstructionDummy::CalorimeterConstructionDummy()
  :fSolidCaloBox(0),fSolidScintBox(0),fSolidCrystalBox(0),fSolidLayer(0),fSolidVetoLatX(0),fSolidVetoLatY(0),fSolidVetoBottom(0),fSolidVetoUp(0),
   fLogicCaloBox(0),fLogicScintBox(0),fLogicCrystalBox(0),fLogicLayer(0),fLogicScintLayer(0),fLogicCrystalLayer(0),fLogicVetoLatX(0),fLogicVetoLatY(0),fLogicVetoBottom(0),fLogicVetoUp(0),
   fPhysiCaloBox(0),fPhysiScintBox(0),fPhysiCrystalBox(0),fPhysiLayer(0),fPhysiScintLayer(0),fPhysiCrystalLayer(0),fPhysiVetoLatXP(0),fPhysiVetoLatXM(0),fPhysiVetoLatYP(0),
  fPhysiVetoLatYM(0),fPhysiVetoBottom(0),fPhysiVetoUp(0)
{
  pMaterial           = new HEPDSWMaterial();

  fNbOfCaloLayers = 20;
  fNbOfCrystalLayers = 0;
  fNbOfScintLayers = fNbOfCaloLayers-fNbOfCrystalLayers;
  
  fCaloLayer_X    = 150*mm;
  fCaloLayer_Y    = 150*mm;
  fCaloLayer_Z    = 10*mm;

  fCalo_X = 150*mm;
  fCalo_Y = 150*mm;
  fCalo_Z = fCaloLayer_Z*fNbOfCaloLayers;

  fScint_X = 150*mm;
  fScint_Y = 150*mm;
  fScint_Z = fCaloLayer_Z*fNbOfScintLayers;

  fCrystal_X = 150*mm;
  fCrystal_Y = 150*mm;
  fCrystal_Z = fCaloLayer_Z*fNbOfCrystalLayers;

  fVetoUp_X = 160*mm;
  fVetoUp_Y = 160*mm;
  fVetoUp_Z = 10*mm;

  fVetoLatX_X = 160*mm;
  fVetoLatX_Y = 5*mm;
  fVetoLatX_Z = fCalo_Z;

  fVetoLatY_X = 5*mm;
  fVetoLatY_Y = 150*mm;
  fVetoLatY_Z = fCalo_Z;

  fVetoBottom_X = 160*mm;
  fVetoBottom_Y = 160*mm;
  fVetoBottom_Z = 5*mm;


  ComputeObjectsPositioning();
  
  // materials
  vetoMaterial      = "Scintillator";
  crystalMaterial   = "Scintillator";
  scintMaterial     = "Scintillator";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterConstructionDummy::~CalorimeterConstructionDummy()
{
  if (pMaterial) delete pMaterial;
}


void CalorimeterConstructionDummy::ComputeObjectsPositioning(){

  fPhysiVetoLatX_X = 0;
  fPhysiVetoLatX_Y = fCalo_Y/2.+fVetoLatX_Y/2.;
  fPhysiVetoLatX_Z = 0;

  fPhysiVetoLatY_X = fCalo_X/2.+fVetoLatY_X/2.;
  fPhysiVetoLatY_Y = 0;
  fPhysiVetoLatY_Z = 0;

  fPhysiVetoBottom_X = 0;
  fPhysiVetoBottom_Y = 0;
  fPhysiVetoBottom_Z = -fCalo_Z/2.-fVetoBottom_Z/2.;

  fPhysiVetoUp_X = 0;
  fPhysiVetoUp_Y = 0;
  fPhysiVetoUp_Z = fCalo_Z/2.+fVetoUp_Z/2.;

  fPhysiScintBox_X = 0;
  fPhysiScintBox_Y = 0;
  fPhysiScintBox_Z = fCalo_Z/2.-fScint_Z/2.;

  fPhysiCrystalBox_X = 0;
  fPhysiCrystalBox_Y = 0;
  fPhysiCrystalBox_Z = fCalo_Z/2.-fScint_Z-fCrystal_Z/2.;

}

void CalorimeterConstructionDummy::SetNumberOfCrystalLayer(G4int aVal){

  fNbOfCrystalLayers = aVal;
  fNbOfScintLayers = fNbOfCaloLayers-fNbOfCrystalLayers;

  fScint_X = 150*mm;
  fScint_Y = 150*mm;
  fScint_Z = fCaloLayer_Z*fNbOfScintLayers;

  fCrystal_X = 150*mm;
  fCrystal_Y = 150*mm;
  fCrystal_Z = fCaloLayer_Z*fNbOfCrystalLayers;

  ComputeObjectsPositioning();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterConstructionDummy::Builder(G4VPhysicalVolume* motherVolume)
{
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String caloSDname = "/hepd/calorimeter";
  CalorimeterSD * caloSD = new CalorimeterSD(caloSDname);
  SDman->AddNewDetector(caloSD);

  G4String vetoSDname = "/hepd/veto";
  VetoSD * vetoSD = new VetoSD(vetoSDname);
  SDman->AddNewDetector(vetoSD);

  pMaterial->DefineMaterials();
  G4Material* vacuum = pMaterial->GetMaterial("Galactic");
  G4Material* scintLayerMat = pMaterial->GetMaterial(scintMaterial);
  G4Material* crystalLayerMat = pMaterial->GetMaterial(crystalMaterial);
  G4Material* vetoLayerMat = pMaterial->GetMaterial(vetoMaterial);

  fSolidCaloBox = new G4Box("Calorimeter",fCalo_X/2.,fCalo_Y/2.,fCalo_Z/2.);
  
  if(fNbOfCrystalLayers>0){
    fSolidScintBox = new G4Box("CalorimeterScint",fScint_X/2.,fScint_Y/2.,fScint_Z/2.);
    fSolidCrystalBox = new G4Box("CalorimeterCrystal",fCrystal_X/2.,fCrystal_Y/2.,fCrystal_Z/2.);
  }

  fSolidLayer = new G4Box("ActiveLayer",fCaloLayer_X/2.,fCaloLayer_Y/2.,fCaloLayer_Z/2.);

  fSolidVetoLatX = new G4Box("Veto",fVetoLatX_X/2.,fVetoLatX_Y/2.,fVetoLatX_Z/2.);

  fSolidVetoLatY = new G4Box("Veto",fVetoLatY_X/2.,fVetoLatY_Y/2.,fVetoLatY_Z/2.);

  fSolidVetoBottom = new G4Box("Veto",fVetoBottom_X/2.,fVetoBottom_Y/2.,fVetoBottom_Z/2.);

  fSolidVetoUp = new G4SubtractionSolid("Veto",new G4Box("Veto",fVetoUp_X/2.,fVetoUp_Y/2.,fVetoUp_Z/2.),new G4Box("Veto",fVetoUp_X/2.-fVetoLatX_Y,fVetoUp_Y/2.-fVetoLatX_Y,fVetoUp_Z/2.));

  fLogicCaloBox = new G4LogicalVolume(fSolidCaloBox,vacuum,"Calorimeter");

  if(fNbOfCrystalLayers>0){
    fLogicScintBox = new G4LogicalVolume(fSolidScintBox,vacuum,"CalorimeterScint");
    fLogicCrystalBox = new G4LogicalVolume(fSolidCrystalBox,vacuum,"CalorimeterCrystall");
    fLogicScintLayer = new G4LogicalVolume(fSolidLayer,scintLayerMat,"ActiveLayerScint");
    fLogicCrystalLayer = new G4LogicalVolume(fSolidLayer,crystalLayerMat,"ActiveLayerCrystal"); 
  }else{
    fLogicLayer = new G4LogicalVolume(fSolidLayer,scintLayerMat,"ActiveLayer");
  }
  fLogicVetoLatX = new G4LogicalVolume(fSolidVetoLatX,vetoLayerMat,"Veto");

  fLogicVetoLatY = new G4LogicalVolume(fSolidVetoLatY,vetoLayerMat,"Veto");

  fLogicVetoBottom = new G4LogicalVolume(fSolidVetoBottom,vetoLayerMat,"Veto");

  fLogicVetoUp = new G4LogicalVolume(fSolidVetoUp,vetoLayerMat,"Veto");

  if(fNbOfCrystalLayers>0){
    fLogicScintLayer->SetSensitiveDetector(caloSD);
    fLogicCrystalLayer->SetSensitiveDetector(caloSD);
  }else{
    fLogicLayer->SetSensitiveDetector(caloSD);
  }
  fLogicVetoBottom->SetSensitiveDetector(vetoSD);
  fLogicVetoUp->SetSensitiveDetector(vetoSD);
  fLogicVetoLatY->SetSensitiveDetector(vetoSD);
  fLogicVetoLatX->SetSensitiveDetector(vetoSD);

  fPhysiCaloBox = new G4PVPlacement(0,                        
				    G4ThreeVector(0,0,0),   
				    "Calorimeter",               
				    fLogicCaloBox,                
				    motherVolume,                
				    false,                       
				    0,true);       
  if(fNbOfCrystalLayers>0){
    fPhysiScintBox = new G4PVPlacement(0,                        
				       G4ThreeVector(fPhysiScintBox_X,fPhysiScintBox_Y,fPhysiScintBox_Z),   
				       "CalorimeterScint",               
				       fLogicScintBox,                
				       fPhysiCaloBox,
				       false,                       
				       0,true); 

    fPhysiScintLayer = new G4PVReplica("ActiveLayerScint", 
				       fLogicScintLayer,            
				       fPhysiScintBox,            
				       kZAxis,                 

				       fNbOfScintLayers,        
				       fCaloLayer_Z);   
    
    fPhysiCrystalBox = new G4PVPlacement(0,                        
					 G4ThreeVector(fPhysiCrystalBox_X,fPhysiCrystalBox_Y,fPhysiCrystalBox_Z),     
					 "CalorimeterCrystall",               
					 fLogicCrystalBox,                
					 fPhysiCaloBox,
					 false,                       
					 0,true); 
    
    fPhysiCrystalLayer = new G4PVReplica("ActiveLayerCrystal", 
					 fLogicCrystalLayer,            
					 fPhysiCrystalBox,            
					 kZAxis,                 
					 fNbOfCrystalLayers,        
					 fCaloLayer_Z);   
    
  }else{
    fPhysiLayer = new G4PVReplica("ActiveLayer", 
				  fLogicLayer,            
				  fPhysiCaloBox,            
				  kZAxis,                 
				  fNbOfCaloLayers,        
				  fCaloLayer_Z);   
  }
  
  fPhysiVetoLatXP = new G4PVPlacement(0,
				      G4ThreeVector(fPhysiVetoLatX_X,fPhysiVetoLatX_Y,fPhysiVetoLatX_Z),     
				      "VetoLXP",              
				      fLogicVetoLatX,       
				      motherVolume,        
				      false,               
				      0,true); 

  fPhysiVetoLatXM = new G4PVPlacement(0,
				      G4ThreeVector(fPhysiVetoLatX_X,-fPhysiVetoLatX_Y,fPhysiVetoLatX_Z),     
				      "VetoLXM",              
				      fLogicVetoLatX,       
				      motherVolume,        
				      false,               
				      0,true); 
  
  fPhysiVetoLatYP = new G4PVPlacement(0,
				      G4ThreeVector(fPhysiVetoLatY_X,fPhysiVetoLatY_Y,fPhysiVetoLatY_Z),     
				      "VetoLYP",              
				      fLogicVetoLatY,       
				      motherVolume,        
				      false,               
				      0,true); 

  fPhysiVetoLatYM = new G4PVPlacement(0,
				      G4ThreeVector(-fPhysiVetoLatY_X,fPhysiVetoLatY_Y,fPhysiVetoLatY_Z),     
				      "VetoLYM",              
				      fLogicVetoLatY,       
				      motherVolume,        
				      false,               
				      0,true); 
  
  
  fPhysiVetoBottom = new G4PVPlacement(0,
				       G4ThreeVector(fPhysiVetoBottom_X,fPhysiVetoBottom_Y,fPhysiVetoBottom_Z),
				       "VetoB",
				       fLogicVetoBottom,
				       motherVolume,
				       false,
				       0,true);

  fPhysiVetoUp = new G4PVPlacement(0,
				   G4ThreeVector(fPhysiVetoUp_X,fPhysiVetoUp_Y,fPhysiVetoUp_Z),
				   "VetoU",
				   fLogicVetoUp,
				   motherVolume,
				   false,
				   0,true);


  //Visualization Attribute
  G4VisAttributes* attMagenta = new G4VisAttributes(G4Colour::Magenta());
  attMagenta->SetVisibility(true);
  attMagenta->SetForceAuxEdgeVisible(true);
  if(fNbOfCrystalLayers>0)
    fLogicScintLayer->SetVisAttributes(attMagenta);
  else
    fLogicLayer->SetVisAttributes(attMagenta);
  
  G4VisAttributes* attYellow = new G4VisAttributes(G4Colour::Yellow());
  attYellow->SetVisibility(true);
  attYellow->SetForceAuxEdgeVisible(true);
  if(fNbOfCrystalLayers>0)
    fLogicCrystalLayer->SetVisAttributes(attYellow);

  G4VisAttributes* attGreen = new G4VisAttributes(G4Colour::Green());
  attGreen->SetVisibility(true);
  attGreen->SetForceAuxEdgeVisible(true);
  fLogicVetoBottom->SetVisAttributes(attGreen);
  fLogicVetoUp->SetVisAttributes(attGreen);
  fLogicVetoLatY->SetVisAttributes(attGreen);
  fLogicVetoLatX->SetVisAttributes(attGreen);

  G4VisAttributes * attInvisible = new G4VisAttributes();
  attInvisible->SetVisibility(false);
  attInvisible->SetForceAuxEdgeVisible(false);
  fLogicCaloBox->SetVisAttributes(attInvisible);
  if(fNbOfCrystalLayers>0){
    fLogicScintBox->SetVisAttributes(attInvisible);
    fLogicCrystalBox->SetVisAttributes(attInvisible);
  } 
}


