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

#include "HEPDBoxConstructionConfig3.hh"


#include "G4NistManager.hh"
#include "HEPDSWMaterial.hh"
#include "G4Box.hh"
#include "G4UnionSolid.hh"
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

HEPDBoxConstructionConfig3::HEPDBoxConstructionConfig3()
  :fSolidBlanket(0),fLogicBlanket(0),fPhysiBlanket(0)
{
  pMaterial     = new HEPDSWMaterial();
  fBlanket_X  = 277.135*mm;
  fBlanket_Y  = 256.980*mm;
  fBlanket_Z  = 0.21*mm;

  fSuteCoat_X  = 277.135*mm;
  fSuteCoat_Y  = 256.980*mm;
  fSuteCoat_Z  = .004*mm;
              
  fSuteThick_X  = 277.135*mm;
  fSuteThick_Y  = 256.980*mm;
  fSuteThick_Z  = .018*mm;
              
  fSuteThin_X  = 277.135*mm;
  fSuteThin_Y  = 256.980*mm;
  fSuteThin_Z  = .006*mm;
              
  fSoftMat_X  = 277.135*mm;
  fSoftMat_Y  = 256.980*mm;
  fSoftMat_Z  = .004*mm;
              
  fSupportO_X = 30.0*mm;
  fSupportO_Y = 256.980*mm;
  fSupportO_Z = .1*mm;
              
  fSupportV_X = 277.135*mm;
  fSupportV_Y = 30.0*mm;
  fSupportV_Z = .1*mm;

  ComputeObjectsPositioning();
  
  // materials
  blanketMaterial   = "Galactic";
  suteMaterial      = "mylar";
  suteCoatMaterial  = "Aluminium";
  softMaterial      = "Dacron";
  supportMaterial   = "Kapton";

}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HEPDBoxConstructionConfig3::~HEPDBoxConstructionConfig3()
{
  if (pMaterial) delete pMaterial;
}

void HEPDBoxConstructionConfig3::SetBlanketMaterial(G4String asuteMat,G4String asuteCoatMat,G4String asoftMat,G4String asupportMat){
  asuteMat=suteMaterial; 
    asuteCoatMat=suteCoatMaterial; 
    asoftMat=softMaterial; 
    asupportMat=supportMaterial;
}    
void HEPDBoxConstructionConfig3::ComputeObjectsPositioning(){

  fBlanketPos_Z    = 207.035*mm;
  fSuteCoat1Pos_Z  = fBlanket_Z/2.-fSuteCoat_Z/2.;  
  fSuteThick1Pos_Z = fSuteCoat1Pos_Z-fSuteCoat_Z/2.-fSuteThick_Z/2.; 
  fSoftMat1Pos_Z   = fSuteThick1Pos_Z-fSuteThick_Z/2.-fSoftMat_Z/2.; 
  fSuteCoat2Pos_Z  = fSoftMat1Pos_Z-fSoftMat_Z/2.-fSuteCoat_Z/2.;
  fSuteThin1Pos_Z  = fSuteCoat2Pos_Z-fSuteCoat_Z/2.-fSuteThin_Z/2.;  
  fSuteCoat3Pos_Z  = fSuteThin1Pos_Z-fSuteThin_Z/2.-fSuteCoat_Z/2.;  
  fSoftMat2Pos_Z   = fSuteCoat3Pos_Z-fSuteCoat_Z/2.-fSoftMat_Z/2.;
  fSuteCoat4Pos_Z  = fSoftMat2Pos_Z-fSoftMat_Z/2.-fSuteCoat_Z/2.;  
  fSuteThin2Pos_Z  = fSuteCoat4Pos_Z-fSuteCoat_Z/2.-fSuteThin_Z/2.;  
  fSuteCoat5Pos_Z  = fSuteThin2Pos_Z-fSuteThin_Z/2.-fSuteCoat_Z/2.;  
  fSoftMat3Pos_Z   = fSuteCoat5Pos_Z-fSuteCoat_Z/2.-fSoftMat_Z/2.;   
  fSuteCoat6Pos_Z  = fSoftMat3Pos_Z-fSoftMat_Z/2.-fSuteCoat_Z/2.;  
  fSuteThin3Pos_Z  = fSuteCoat6Pos_Z-fSuteCoat_Z/2.-fSuteThin_Z/2.;  
  fSuteCoat7Pos_Z  = fSuteThin3Pos_Z-fSuteThin_Z/2.-fSuteCoat_Z/2.;  
  fSoftMat4Pos_Z   = fSuteCoat7Pos_Z-fSuteCoat_Z/2.-fSoftMat_Z/2.;   
  fSuteCoat8Pos_Z  = fSoftMat4Pos_Z-fSoftMat_Z/2.-fSuteCoat_Z/2.;  
  fSuteThick2Pos_Z = fSuteCoat8Pos_Z-fSuteCoat_Z/2.-fSuteThick_Z/2.; 
  fSuteCoat9Pos_Z  = fSuteThick2Pos_Z-fSuteThick_Z/2-fSuteCoat_Z/2.;  
  fSupportPos_Z    = fSuteCoat9Pos_Z-fSuteCoat_Z/2.-fSupportV_Z/2.;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
void HEPDBoxConstructionConfig3::Builder(G4VPhysicalVolume* motherVolume)
{

  pMaterial->DefineMaterials();
  G4Material* blanketMat = pMaterial->GetMaterial(blanketMaterial);
  G4Material* suteMat = pMaterial->GetMaterial(suteMaterial);
  G4Material* suteCoatMat = pMaterial->GetMaterial(suteCoatMaterial);
  G4Material* softMat = pMaterial->GetMaterial(softMaterial);
  G4Material* supportMat = pMaterial->GetMaterial(supportMaterial);
  
  fSolidBlanket = new G4Box("ThermalBlanket",fBlanket_X/2.,fBlanket_Y/2.,fBlanket_Z/2);

  fSolidSuteCoat = new G4Box("SuteCoating",fSuteCoat_X/2.,fSuteCoat_Y/2.,fSuteCoat_Z/2.);

  fSolidSuteThin = new G4Box("SuteThin",fSuteThin_X/2.,fSuteThin_Y/2.,fSuteThin_Z/2.);

  fSolidSuteThick = new G4Box("SuteThick",fSuteThick_X/2.,fSuteThick_Y/2.,fSuteThick_Z/2.);

  fSolidSoftMat = new G4Box("SoftMat",fSoftMat_X/2.,fSoftMat_Y/2.,fSoftMat_Z/2.);

  fSolidSupportO = new G4Box("SupportO",fSupportO_X/2.,fSupportO_Y/2.,fSupportO_Z/2.);

  fSolidSupportV = new G4Box("SupportV",fSupportV_X/2.,fSupportV_Y/2.,fSupportV_Z/2.);

  fSolidSupport = new G4UnionSolid("Support",fSolidSupportO,fSolidSupportV);
  

  fLogicBlanket = new G4LogicalVolume(fSolidBlanket,
				      blanketMat,
				      "ThermalBlanket");

  fLogicSuteCoat = new G4LogicalVolume(fSolidSuteCoat,
				      suteCoatMat,
				      "SuteCoating");

  fLogicSuteThin = new G4LogicalVolume(fSolidSuteThin,
				      suteMat,
				      "SuteThin");

  fLogicSuteThick = new G4LogicalVolume(fSolidSuteThick,
				      suteMat,
				      "SuteThick");

  fLogicSoftMat = new G4LogicalVolume(fSolidSoftMat,
				      softMat,
				      "SoftMat");

  fLogicSupport = new G4LogicalVolume(fSolidSupport,
				      supportMat,
				      "Support");

  fPhysiBlanket = new G4PVPlacement(0,
				    G4ThreeVector(0.,0.,fBlanketPos_Z),
				    "ThermalBlanket",
				    fLogicBlanket,                
				    motherVolume,                
				    false,                       
				    0,true);       

  fPhysiSuteCoat1 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSuteCoat1Pos_Z),
				      "SuteCoating",
				      fLogicSuteCoat,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       

  fPhysiSuteThick1 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSuteThick1Pos_Z),
				      "SuteThick",
				      fLogicSuteThick,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       

  fPhysiSoftMat1 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSoftMat1Pos_Z),
				      "SoftMat",
				      fLogicSoftMat,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       

  fPhysiSuteCoat2 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSuteCoat2Pos_Z),
				      "SuteCoating",
				      fLogicSuteCoat,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       

  fPhysiSuteThin1 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSuteThin1Pos_Z),
				      "SuteThin",
				      fLogicSuteThin,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       


  fPhysiSuteCoat3 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSuteCoat3Pos_Z),
				      "SuteCoating",
				      fLogicSuteCoat,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       

  fPhysiSoftMat2 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSoftMat2Pos_Z),
				      "SoftMat",
				      fLogicSoftMat,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       

  fPhysiSuteCoat4 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSuteCoat4Pos_Z),
				      "SuteCoating",
				      fLogicSuteCoat,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       

  fPhysiSuteThin2 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSuteThin2Pos_Z),
				      "SuteThin",
				      fLogicSuteThin,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       

  fPhysiSuteCoat5 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSuteCoat5Pos_Z),
				      "SuteCoating",
				      fLogicSuteCoat,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       

  fPhysiSoftMat3 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSoftMat3Pos_Z),
				      "SoftMat",
				      fLogicSoftMat,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       

  fPhysiSuteCoat6 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSuteCoat6Pos_Z),
				      "SuteCoating",
				      fLogicSuteCoat,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       

  fPhysiSuteThin3 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSuteThin3Pos_Z),
				      "SuteThin",
				      fLogicSuteThin,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       

  fPhysiSuteCoat7 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSuteCoat7Pos_Z),
				      "SuteCoating",
				      fLogicSuteCoat,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       

  fPhysiSoftMat4 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSoftMat4Pos_Z),
				      "SoftMat",
				      fLogicSoftMat,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       

  fPhysiSuteCoat8 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSuteCoat8Pos_Z),
				      "SuteCoating",
				      fLogicSuteCoat,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);   

  fPhysiSuteThick2 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSuteThick2Pos_Z),
				      "SuteThick",
				      fLogicSuteThick,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);       

  fPhysiSuteCoat9 = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSuteCoat9Pos_Z),
				      "SuteCoating",
				      fLogicSuteCoat,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);   
    
  fPhysiSupport = new G4PVPlacement(0,
				      G4ThreeVector(0,0,fSupportPos_Z),
				      "Support",
				      fLogicSupport,                
				      fPhysiBlanket,                
				      false,                       
				      0,true);   

  //Visualization Attribute
  G4VisAttributes * attInvisible = new G4VisAttributes();
  attInvisible->SetVisibility(false);
  attInvisible->SetForceAuxEdgeVisible(false);
  fLogicBlanket->SetVisAttributes(attInvisible);

  G4VisAttributes* attYellow = new G4VisAttributes(G4Colour::Yellow());
  attYellow->SetVisibility(true);
  attYellow->SetForceAuxEdgeVisible(true);
  fLogicSuteThin->SetVisAttributes(attYellow);
  fLogicSuteThick->SetVisAttributes(attYellow);

  G4VisAttributes* attCyan = new G4VisAttributes(G4Colour::Cyan());
  attCyan->SetVisibility(true);
  attCyan->SetForceAuxEdgeVisible(true);
  fLogicSoftMat->SetVisAttributes(attCyan);

  G4VisAttributes* attRed = new G4VisAttributes(G4Colour::Red());
  attRed->SetVisibility(true);
  attRed->SetForceAuxEdgeVisible(true);
  fLogicSupport->SetVisAttributes(attRed);

  G4VisAttributes* attGray = new G4VisAttributes(G4Colour::Gray());
  attGray->SetVisibility(true);
  attGray->SetForceAuxEdgeVisible(true);
  fLogicSuteCoat->SetVisAttributes(attGray);


}


