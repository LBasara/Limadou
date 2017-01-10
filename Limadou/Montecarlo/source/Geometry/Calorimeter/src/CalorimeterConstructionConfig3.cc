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

#include "CalorimeterConstructionConfig3.hh"

#include "G4NistManager.hh"
#include "HEPDSWMaterial.hh"
#include "G4Box.hh"
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

CalorimeterConstructionConfig3::CalorimeterConstructionConfig3()
  :fSolidCaloBox(0),fSolidScintBox(0),fSolidCrystalBox(0),fSolidLayer(0),fSolidActiveLayer(0),
   fSolidCFFront(0),fSolidCFLat(0),fSolidCFSuppO(0),fSolidCFSuppV(0),fSolidCFSuppLat(0),fSolidCFSuppFront(0),
   fSolidCFSuppStepOV(0),fSolidCFSuppStepLat1(0),fSolidCFSuppStepLat2(0),fSolidCFSuppStepFront1(0),fSolidCFSupp(0),
  fSolidPoronLat(0),fSolidPoronLatX(0),fSolidVetoLatX(0),fSolidPoronLatY(0),fSolidVetoLatY(0),fSolidPoronBottom(0),
   fSolidPoronBottomLatStick(0),fSolidVetoBottom(0),fSolidCFBottomLat(0),fSolidCFBottomFront(0),fSolidHoneyCombLatY(0),
  fSolidHoneyCombSkinLatY(0),fSolidHoneyCombLatX(0),fSolidHoneyCombSkinLatX(0),fSolidHoneyCombBottom(0),fSolidHoneyCombSkinBottom(0),
  fLogicCaloBox(0),fLogicScintBox(0),fLogicCrystalBox(0),fLogicLayer(0),fLogicScintLayer(0),fLogicCrystalLayer(0),fLogicActiveLayer(0),fLogicScintActiveLayer(0),fLogicCrystalActiveLayer(0),fLogicCFFront(0),fLogicCFLat(0),fLogicCFSupp(0),fLogicCFSuppTop(0),fLogicPoronLat(0),fLogicPoronLatX(0),	   
  fLogicVetoLatX(0),fLogicPoronLatY(0),fLogicVetoLatY(0),fLogicPoronBottom(0),fLogicPoronBottomLatStick(0),fLogicVetoBottom(0),fLogicCFBottomLat(0),	   
  fLogicCFBottomFront(0),fLogicHoneyCombLatY(0),fLogicHoneyCombSkinLatY(0),fLogicHoneyCombLatX(0),fLogicHoneyCombSkinLatX(0),fLogicHoneyCombBottom(0),	   
  fLogicHoneyCombSkinBottom(0),
  fPhysiCaloBox(0),fPhysiLayer(0),fPhysiCFSupp(0),fPhysiActiveLayer(0),fPhysiCFFrontP(0),fPhysiCFFrontM(0),fPhysiCFLatP(0),fPhysiCFLatM(0),fPhysiPoronLatP(0),fPhysiPoronLatM(0),
  fPhysiScintBox(0),fPhysiScintLayer(0),fPhysiScintCFSupp(0),fPhysiScintActiveLayer(0),fPhysiScintCFFrontP(0),fPhysiScintCFFrontM(0),fPhysiScintCFLatP(0),fPhysiScintCFLatM(0),fPhysiScintPoronLatP(0),fPhysiScintPoronLatM(0),
  fPhysiCrystalBox(0),fPhysiCrystalLayer(0),fPhysiCrystalCFSupp(0),fPhysiCrystalActiveLayer(0),fPhysiCrystalCFFrontP(0),fPhysiCrystalCFFrontM(0),fPhysiCrystalCFLatP(0),fPhysiCrystalCFLatM(0),fPhysiCrystalPoronLatP(0),fPhysiCrystalPoronLatM(0),
  fPhysiPoronLatXPInt(0),fPhysiPoronLatXMInt(0),fPhysiPoronLatXPExt(0),fPhysiPoronLatXMExt(0),fPhysiPoronLatYPInt(0),fPhysiPoronLatYMInt(0),
  fPhysiPoronLatYPExt(0),fPhysiPoronLatYMExt(0),fPhysiVetoLatXP(0),fPhysiVetoLatXM(0),fPhysiVetoLatYP(0),fPhysiVetoLatYM(0),fPhysiPoronBottomU(0),
  fPhysiPoronBottomD(0),fPhysiPoronBottomLatStickP(0),fPhysiPoronBottomLatStickM(0),fPhysiVetoBottom(0),fPhysiCFBottomLatP(0),fPhysiCFBottomLatM(0),
  fPhysiCFBottomFrontP(0),fPhysiCFBottomFrontM(0),fPhysiHoneyCombSkinLatYPInt(0),fPhysiHoneyCombLatYP(0),fPhysiHoneyCombSkinLatYPExt(0),fPhysiHoneyCombSkinLatYMInt(0),
  fPhysiHoneyCombLatYM(0),fPhysiHoneyCombSkinLatYMExt(0),fPhysiHoneyCombSkinLatXPInt(0),fPhysiHoneyCombLatXP(0),fPhysiHoneyCombSkinLatXPExt(0),
  fPhysiHoneyCombSkinLatXMInt(0),fPhysiHoneyCombLatXM(0),fPhysiHoneyCombSkinLatXMExt(0),fPhysiHoneyCombSkinBottomU(0),fPhysiHoneyCombBottom(0),fPhysiHoneyCombSkinBottomD(0),
  fPhysiCFSuppTop(0)
{
  pMaterial           = new HEPDSWMaterial();

  fNbOfCaloLayers    = 20;
  fNbOfCrystalLayers = 0;
  fNbOfScintLayers   = fNbOfCaloLayers-fNbOfCrystalLayers;

  fCaloLayer_X    = 177*mm;
  fCaloLayer_Y    = 170*mm;
  fCaloLayer_Z    = 14.78*mm;

  fCalo_X = 177*mm;
  fCalo_Y = 170*mm;
  fCalo_Z = fCaloLayer_Z*fNbOfCaloLayers;

  fScint_X = 177*mm;
  fScint_Y = 170*mm;
  fScint_Z = fCaloLayer_Z*fNbOfScintLayers;

  fCrystal_X = 177*mm;
  fCrystal_Y = 170*mm;
  fCrystal_Z = fCaloLayer_Z*fNbOfCrystalLayers;

  fCFLat_X = 10*mm;
  fCFLat_Y = 150*mm;
  fCFLat_Z = 10*mm;

  fCFFront_X = 30*mm;
  fCFFront_Y = 10*mm;
  fCFFront_Z = 10*mm;

  fCFSuppO_X = 157*mm;
  fCFSuppO_Y = 3*mm;
  fCFSuppO_Z = 4.78*mm;

  fCFSuppV_X = 3*mm;
  fCFSuppV_Y = 150*mm;
  fCFSuppV_Z = 4.78*mm;

  fCFSuppLat_X = 10*mm;
  fCFSuppLat_Y = 150*mm;
  fCFSuppLat_Z = 4.78*mm;
    
  fCFSuppFront_X = 30*mm;
  fCFSuppFront_Y = 10*mm;
  fCFSuppFront_Z = 4.78*mm;
    
  fPoronLat_X = 3.18*mm;
  fPoronLat_Y = 150*mm;
  fPoronLat_Z = 10*mm;

  fActiveLayer_X = 15*cm;
  fActiveLayer_Y = 15*cm;
  fActiveLayer_Z = 1*cm;

  fPoronLatX_X = 199.72*mm;
  fPoronLatX_Y = 3.18*mm;
  fPoronLatX_Z = 322.24*mm;

  fVetoLatX_X = 199.72*mm;
  fVetoLatX_Y = 5*mm;
  fVetoLatX_Z = 322.24*mm;

  fPoronLatY_X = 3.18*mm;
  fPoronLatY_Y = 170*mm;
  fPoronLatY_Z = 322.24*mm;

  fVetoLatY_X = 5*mm;
  fVetoLatY_Y = 170*mm;
  fVetoLatY_Z = 322.24*mm;

  fPoronBottom_X = 157*mm;
  fPoronBottom_Y = 150*mm;
  fPoronBottom_Z = 3.18*mm;

  fPoronBottomLatStick_X = 3.18*mm;
  fPoronBottomLatStick_Y = 150*mm;
  fPoronBottomLatStick_Z = 5*mm;

  fVetoBottom_X = 150*mm;
  fVetoBottom_Y = 150*mm;
  fVetoBottom_Z = 5*mm;

  fCFBottomLat_X = 10*mm;
  fCFBottomLat_Y = 150*mm;
  fCFBottomLat_Z = 11.36*mm;

  fCFBottomFront_X = 50*mm;
  fCFBottomFront_Y = 10*mm;
  fCFBottomFront_Z = 11.36*mm;

  fHoneyCombLatY_X = 10*mm;
  fHoneyCombLatY_Y = 192.72*mm;
  fHoneyCombLatY_Z = 322.24*mm;

  fHoneyCombSkinLatY_X = 1*mm;
  fHoneyCombSkinLatY_Y = 192.72*mm;
  fHoneyCombSkinLatY_Z = 322.24*mm;

  fHoneyCombLatX_X = 199.72*mm;
  fHoneyCombLatX_Y = 10*mm;
  fHoneyCombLatX_Z = 322.24*mm;

  fHoneyCombSkinLatX_X = 199.72*mm;
  fHoneyCombSkinLatX_Y = 1*mm;
  fHoneyCombSkinLatX_Z = 322.24*mm;

  fHoneyCombBottom_X = 300*mm;
  fHoneyCombBottom_Y = 300*mm;
  fHoneyCombBottom_Z = 10*mm;

  fHoneyCombSkinBottom_X = 300*mm;
  fHoneyCombSkinBottom_Y = 300*mm;
  fHoneyCombSkinBottom_Z = 1*mm;

  ComputeObjectsPositioning();
  
  // materials
  vetoMaterial      = "Scintillator";
  caloMaterial      = "Scintillator";
  scintMaterial     = "Scintillator";
  crystalMaterial   = "BGO";

  poronMaterial     = "Poron";
  cfiberMaterial    = "CarbonFiber";
  honeycombMaterial = "nextel";
  /*
    poronMaterial     = "Galactic";
    cfiberMaterial    = "Galactic";
    honeycombMaterial = "Galactic";
  */  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterConstructionConfig3::~CalorimeterConstructionConfig3()
{
  if (pMaterial) delete pMaterial;
}

void CalorimeterConstructionConfig3::SetNumberOfCrystalLayer(G4int aVal){
  
  fNbOfCrystalLayers = aVal;
  fNbOfScintLayers = fNbOfCaloLayers-fNbOfCrystalLayers;
  
  fScint_X = 177*mm;
  fScint_Y = 170*mm;
  fScint_Z = fCaloLayer_Z*fNbOfScintLayers;

  fCrystal_X = 177*mm;
  fCrystal_Y = 170*mm;
  fCrystal_Z = fCaloLayer_Z*fNbOfCrystalLayers;

  ComputeObjectsPositioning();

}

void CalorimeterConstructionConfig3::ComputeObjectsPositioning(){

  fPhysiScintBox_X = 0;
  fPhysiScintBox_Y = 0;
  fPhysiScintBox_Z = fCalo_Z/2.-fScint_Z/2.;

  fPhysiCrystalBox_X = 0;
  fPhysiCrystalBox_Y = 0;
  fPhysiCrystalBox_Z = fCalo_Z/2.-fScint_Z-fCrystal_Z/2.;

  fSolidCFSuppStepLat_X = fCaloLayer_X/2.-fCFLat_X/2.;
  fSolidCFSuppStepLat_Y = 0;
  fSolidCFSuppStepLat_Z = 0;

  fSolidCFSuppStepFront_X = 0;
  fSolidCFSuppStepFront_Y = fCaloLayer_Y/2.-fCFFront_Y/2.;
  fSolidCFSuppStepFront_Z = 0;

  fPhysiCFSupp_X = 0;
  fPhysiCFSupp_Y = 0;
  fPhysiCFSupp_Z = -fCaloLayer_Z/2.+fCFSuppO_Z/2.;

  fPhysiActiveLayer_X = 0;
  fPhysiActiveLayer_Y = 0;
  fPhysiActiveLayer_Z = fCaloLayer_Z/2.-fActiveLayer_Z/2.;

  fPhysiCFFront_X = 0;
  fPhysiCFFront_Y = fCaloLayer_Y/2.-fCFFront_Y/2.;
  fPhysiCFFront_Z = fCaloLayer_Z/2.-fCFFront_Z/2.;

  fPhysiCFLat_X = fCaloLayer_X/2.-fCFLat_X/2.;
  fPhysiCFLat_Y = 0;
  fPhysiCFLat_Z = fCaloLayer_Z/2.-fCFLat_Z/2.;

  fPhysiPoronLat_X = fCaloLayer_X/2.-fCFLat_X-fPoronLat_X/2.;
  fPhysiPoronLat_Y = 0;
  fPhysiPoronLat_Z = fCaloLayer_Z/2.-fPoronLat_Z/2.;

  fPhysiPoronLatXInt_X = 0;
  fPhysiPoronLatXInt_Y = fCalo_Y/2.+fPoronLatX_Y/2.;
  fPhysiPoronLatXInt_Z = 1.96*mm;

  fPhysiPoronLatXExt_X = 0;
  fPhysiPoronLatXExt_Y = fCalo_Y/2.+fPoronLatX_Y+fVetoLatX_Y+fPoronLatX_Y/2.;
  fPhysiPoronLatXExt_Z = 1.96*mm;

  fPhysiPoronLatYInt_X = fCalo_X/2.+fPoronLatY_X/2.;
  fPhysiPoronLatYInt_Y = 0;
  fPhysiPoronLatYInt_Z = 1.96*mm;

  fPhysiPoronLatYExt_X = fCalo_X/2.+fPoronLatY_X+fVetoLatY_X+fPoronLatY_X/2.;
  fPhysiPoronLatYExt_Y = 0;
  fPhysiPoronLatYExt_Z = 1.96*mm;

  fPhysiVetoLatX_X = 0;
  fPhysiVetoLatX_Y = fCalo_Y/2.+fPoronLatX_Y+fVetoLatX_Y/2.;
  fPhysiVetoLatX_Z = 1.96*mm;

  fPhysiVetoLatY_X = fCalo_X/2.+fPoronLatY_X+fVetoLatY_X/2.;
  fPhysiVetoLatY_Y = 0;
  fPhysiVetoLatY_Z = 1.96*mm;

  fPhysiPoronBottomU_X = 0;
  fPhysiPoronBottomU_Y = 0;
  fPhysiPoronBottomU_Z = -fCalo_Z/2.-fPoronBottom_Z/2.;

  fPhysiPoronBottomD_X = 0;
  fPhysiPoronBottomD_Y = 0;
  fPhysiPoronBottomD_Z = -fCalo_Z/2.-fPoronBottom_Z-fVetoBottom_Z-fPoronBottom_Z/2.;

  fPhysiPoronBottomLatStick_X = fCalo_X/2-fCFBottomLat_X-fPoronBottomLatStick_X/2.;
  fPhysiPoronBottomLatStick_Y = 0;
  fPhysiPoronBottomLatStick_Z = -fCalo_Z/2.-fPoronBottom_Z-fPoronBottomLatStick_Z/2.;

  fPhysiVetoBottom_X = 0;
  fPhysiVetoBottom_Y = 0;
  fPhysiVetoBottom_Z = -fCalo_Z/2.-fPoronBottom_Z-fVetoBottom_Z/2.;

  fPhysiCFBottomLat_X = fCalo_X/2-fCFBottomLat_X/2.;
  fPhysiCFBottomLat_Y = 0;
  fPhysiCFBottomLat_Z = -fCalo_Z/2.-fCFBottomLat_Z/2.;

  fPhysiCFBottomFront_X = 0;
  fPhysiCFBottomFront_Y = fCalo_Y/2-fCFBottomFront_Y/2.;
  fPhysiCFBottomFront_Z = -fCalo_Z/2.-fCFBottomFront_Z/2.;

  fPhysiHoneyCombSkinLatYInt_X = fPhysiPoronLatYExt_X+fPoronLatY_X/2.+fHoneyCombSkinLatY_X/2.;
  fPhysiHoneyCombSkinLatYInt_Y = 0;
  fPhysiHoneyCombSkinLatYInt_Z = 1.96*mm;

  fPhysiHoneyCombSkinLatYExt_X = fPhysiHoneyCombSkinLatYInt_X+fHoneyCombSkinLatY_X+fHoneyCombLatY_X;
  fPhysiHoneyCombSkinLatYExt_Y = 0;
  fPhysiHoneyCombSkinLatYExt_Z = 1.96*mm;

  fPhysiHoneyCombLatY_X = fPhysiHoneyCombSkinLatYInt_X+fHoneyCombSkinLatY_X/2.+fHoneyCombLatY_X/2.;
  fPhysiHoneyCombLatY_Y = 0;
  fPhysiHoneyCombLatY_Z = 1.96*mm;

  fPhysiHoneyCombSkinLatXInt_X = 0;
  fPhysiHoneyCombSkinLatXInt_Y = fPhysiPoronLatXExt_Y+fPoronLatX_Y/2.+fHoneyCombSkinLatX_Y/2.;
  fPhysiHoneyCombSkinLatXInt_Z = 1.96*mm;

  fPhysiHoneyCombSkinLatXExt_X = 0;
  fPhysiHoneyCombSkinLatXExt_Y = fPhysiHoneyCombSkinLatXInt_Y+fHoneyCombSkinLatX_Y+fHoneyCombLatX_Y;
  fPhysiHoneyCombSkinLatXExt_Z = 1.96*mm;

  fPhysiHoneyCombLatX_X = 0;
  fPhysiHoneyCombLatX_Y = fPhysiHoneyCombSkinLatXInt_Y+fHoneyCombSkinLatX_Y/2.+fHoneyCombLatX_Y/2.;
  fPhysiHoneyCombLatX_Z = 1.96*mm;

  fPhysiHoneyCombSkinBottomU_X = 0;
  fPhysiHoneyCombSkinBottomU_Y = 0;
  fPhysiHoneyCombSkinBottomU_Z = -fCalo_Z/2.-fCFBottomFront_Z-fHoneyCombSkinBottom_Z/2.;

  fPhysiHoneyCombSkinBottomD_X = 0;
  fPhysiHoneyCombSkinBottomD_Y = 0;
  fPhysiHoneyCombSkinBottomD_Z = -fCalo_Z/2.-fCFBottomFront_Z-fHoneyCombSkinBottom_Z-fHoneyCombBottom_Z-fHoneyCombSkinBottom_Z/2.;

  fPhysiHoneyCombBottom_X = 0;
  fPhysiHoneyCombBottom_Y = 0;
  fPhysiHoneyCombBottom_Z = -fCalo_Z/2.-fCFBottomFront_Z-fHoneyCombSkinBottom_Z-fHoneyCombBottom_Z/2.;

  fPhysiCFSuppTop_X = 0;
  fPhysiCFSuppTop_Y = 0;
  fPhysiCFSuppTop_Z = fCalo_Z/2.+fCFSuppO_Z/2.;

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterConstructionConfig3::Builder(G4VPhysicalVolume* motherVolume)
{
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String caloSDname = "/hepd/calorimeter";
  CalorimeterSD * caloSD = new CalorimeterSD(caloSDname);
  SDman->AddNewDetector(caloSD);

  G4String vetoSDname = "/hepd/veto";
  VetoSD * vetoSD = new VetoSD(vetoSDname);
  SDman->AddNewDetector(vetoSD);

  pMaterial->DefineMaterials();
  G4Material* vacuum          = pMaterial->GetMaterial("Galactic");
  G4Material* caloLayerMat    = pMaterial->GetMaterial(caloMaterial);
  G4Material* scintLayerMat   = pMaterial->GetMaterial(scintMaterial);
  G4Material* crystalLayerMat = pMaterial->GetMaterial(crystalMaterial);
  G4Material* vetoLayerMat    = pMaterial->GetMaterial(vetoMaterial);
  G4Material* cfMat           = pMaterial->GetMaterial(cfiberMaterial);
  G4Material* porMat          = pMaterial->GetMaterial(poronMaterial);
  G4Material* hcMat           = pMaterial->GetMaterial(honeycombMaterial);


  fSolidCaloBox = new G4Box("Calorimeter",fCalo_X/2.,fCalo_Y/2.,fCalo_Z/2.);

  if(fNbOfCrystalLayers>0){
    fSolidScintBox = new G4Box("CalorimeterScint",fScint_X/2.,fScint_Y/2.,fScint_Z/2.);
    fSolidCrystalBox = new G4Box("CalorimeterCrystal",fCrystal_X/2.,fCrystal_Y/2.,fCrystal_Z/2.);
  }

  fSolidLayer = new G4Box("Layer",fCaloLayer_X/2.,fCaloLayer_Y/2.,fCaloLayer_Z/2.);

  fSolidActiveLayer = new G4Box("ActiveLayer",fActiveLayer_X/2.,fActiveLayer_Y/2.,fActiveLayer_Z/2.);  

  fSolidCFFront = new  G4Box("CFSupport",fCFFront_X/2.,fCFFront_Y/2.,fCFFront_Z/2.); 

  fSolidCFLat = new  G4Box("CFSupport",fCFLat_X/2.,fCFLat_Y/2.,fCFLat_Z/2.); 

  fSolidCFSuppO = new G4Box("CFSupport",fCFSuppO_X/2.,fCFSuppO_Y/2.,fCFSuppO_Z/2.);

  fSolidCFSuppV = new G4Box("CFSupport",fCFSuppV_X/2.,fCFSuppV_Y/2.,fCFSuppV_Z/2.);

  fSolidCFSuppLat = new G4Box("CFSupport",fCFSuppLat_X/2.,fCFSuppLat_Y/2.,fCFSuppLat_Z/2.);

  fSolidCFSuppFront = new G4Box("CFSupport",fCFSuppFront_X/2.,fCFSuppFront_Y/2.,fCFSuppFront_Z/2.);

  fSolidCFSuppStepOV = new G4UnionSolid("CFSupport",fSolidCFSuppO,fSolidCFSuppV);

  fSolidCFSuppStepLat1 = new G4UnionSolid("CFSupport",fSolidCFSuppStepOV,fSolidCFSuppLat,0,
					  G4ThreeVector(fSolidCFSuppStepLat_X,fSolidCFSuppStepLat_Y,fSolidCFSuppStepLat_Z));

  fSolidCFSuppStepLat2 = new G4UnionSolid("CFSupport",fSolidCFSuppStepLat1,fSolidCFSuppLat,0,
					  G4ThreeVector(-fSolidCFSuppStepLat_X,fSolidCFSuppStepLat_Y,fSolidCFSuppStepLat_Z));

  fSolidCFSuppStepFront1 = new G4UnionSolid("CFSupport",fSolidCFSuppStepLat2,fSolidCFSuppFront,0,
					    G4ThreeVector(fSolidCFSuppStepFront_X,fSolidCFSuppStepFront_Y,fSolidCFSuppStepFront_Z));

  fSolidCFSupp = new G4UnionSolid("CFSupport",fSolidCFSuppStepFront1,fSolidCFSuppFront,0,
				  G4ThreeVector(fSolidCFSuppStepFront_X,-fSolidCFSuppStepFront_Y,fSolidCFSuppStepFront_Z));

  fSolidPoronLat = new G4Box("Poron",fPoronLat_X/2.,fPoronLat_Y/2.,fPoronLat_Z/2.);

  fSolidPoronLatX = new G4Box("Poron",fPoronLatX_X/2.,fPoronLatX_Y/2.,fPoronLatX_Z/2.);

  fSolidVetoLatX = new G4Box("Veto",fVetoLatX_X/2.,fVetoLatX_Y/2.,fVetoLatX_Z/2.);
  
  fSolidPoronLatY = new G4Box("Poron",fPoronLatY_X/2.,fPoronLatY_Y/2.,fPoronLatY_Z/2.);

  fSolidVetoLatY = new G4Box("Veto",fVetoLatY_X/2.,fVetoLatY_Y/2.,fVetoLatY_Z/2.);
  
  fSolidPoronBottom = new G4Box("Poron",fPoronBottom_X/2.,fPoronBottom_Y/2.,fPoronBottom_Z/2.);

  fSolidPoronBottomLatStick = new G4Box("Poron",fPoronBottomLatStick_X/2.,fPoronBottomLatStick_Y/2.,fPoronBottomLatStick_Z/2.);

  fSolidVetoBottom = new G4Box("Veto",fVetoBottom_X/2.,fVetoBottom_Y/2.,fVetoBottom_Z/2.);

  fSolidCFBottomLat = new G4Box("CFSupport",fCFBottomLat_X/2.,fCFBottomLat_Y/2.,fCFBottomLat_Z/2.);
  
  fSolidCFBottomFront = new G4Box("CFSupport",fCFBottomFront_X/2.,fCFBottomFront_Y/2.,fCFBottomFront_Z/2.);

  fSolidHoneyCombLatY = new G4Box("HoneyCombCore",fHoneyCombLatY_X/2.,fHoneyCombLatY_Y/2.,fHoneyCombLatY_Z/2.);

  fSolidHoneyCombSkinLatY = new G4Box("HoneyCombSkin",fHoneyCombSkinLatY_X/2.,fHoneyCombSkinLatY_Y/2.,fHoneyCombSkinLatY_Z/2.);

  fSolidHoneyCombLatX = new G4Box("HoneyCombCore",fHoneyCombLatX_X/2.,fHoneyCombLatX_Y/2.,fHoneyCombLatX_Z/2.);

  fSolidHoneyCombSkinLatX = new G4Box("HoneyCombSkin",fHoneyCombSkinLatX_X/2.,fHoneyCombSkinLatX_Y/2.,fHoneyCombSkinLatX_Z/2.);

  fSolidHoneyCombBottom = new G4Box("HoneyCombCore",fHoneyCombBottom_X/2.,fHoneyCombBottom_Y/2.,fHoneyCombBottom_Z/2.);

  fSolidHoneyCombSkinBottom = new G4Box("HoneyCombSkin",fHoneyCombSkinBottom_X/2.,fHoneyCombSkinBottom_Y/2.,fHoneyCombSkinBottom_Z/2.);


  fLogicCaloBox = new G4LogicalVolume(fSolidCaloBox,vacuum,"Calorimeter");

  if(fNbOfCrystalLayers>0){
    fLogicScintBox = new G4LogicalVolume(fSolidScintBox,vacuum,"CalorimeterScint");
    fLogicCrystalBox = new G4LogicalVolume(fSolidCrystalBox,vacuum,"CalorimeterCrystall");
    fLogicScintLayer = new G4LogicalVolume(fSolidLayer,vacuum,"LayerScint");
    fLogicCrystalLayer = new G4LogicalVolume(fSolidLayer,vacuum,"LayerCrystal"); 
    fLogicScintActiveLayer = new G4LogicalVolume(fSolidActiveLayer,scintLayerMat,"ActiveLayerScint");
    fLogicCrystalActiveLayer = new G4LogicalVolume(fSolidActiveLayer,crystalLayerMat,"ActiveLayerCrystal"); 
  }else{
    fLogicLayer = new G4LogicalVolume(fSolidLayer,vacuum,"Layer");
    fLogicActiveLayer = new G4LogicalVolume(fSolidActiveLayer,caloLayerMat,"ActiveLayer");  
  }
  
  fLogicCFFront = new G4LogicalVolume(fSolidCFFront,cfMat,"CFSupport");

  fLogicCFLat = new G4LogicalVolume(fSolidCFLat,cfMat,"CFSupport");

  fLogicCFSupp = new G4LogicalVolume(fSolidCFSupp,cfMat,"CFSupport");

  fLogicCFSuppTop = new G4LogicalVolume(fSolidCFSupp,cfMat,"CFSupport");

  fLogicPoronLat = new G4LogicalVolume(fSolidPoronLat,porMat,"Poron");

  fLogicPoronLatX = new G4LogicalVolume(fSolidPoronLatX,porMat,"Poron");

  fLogicVetoLatX = new G4LogicalVolume(fSolidVetoLatX,vetoLayerMat,"Veto");

  fLogicPoronLatY = new G4LogicalVolume(fSolidPoronLatY,porMat,"Poron");

  fLogicVetoLatY = new G4LogicalVolume(fSolidVetoLatY,vetoLayerMat,"Veto");

  fLogicPoronBottom = new G4LogicalVolume(fSolidPoronBottom,porMat,"Poron");

  fLogicPoronBottomLatStick = new G4LogicalVolume(fSolidPoronBottomLatStick,porMat,"Poron");

  fLogicVetoBottom = new G4LogicalVolume(fSolidVetoBottom,vetoLayerMat,"Veto");

  fLogicCFBottomLat = new G4LogicalVolume(fSolidCFBottomLat,cfMat,"CFSupport"); 

  fLogicCFBottomFront = new G4LogicalVolume(fSolidCFBottomFront,cfMat,"CFSupport"); 

  fLogicHoneyCombLatY = new G4LogicalVolume(fSolidHoneyCombLatY,hcMat,"HoneyCombCore"); 

  fLogicHoneyCombSkinLatY = new G4LogicalVolume(fSolidHoneyCombSkinLatY,cfMat,"HoneyCombSkin"); 

  fLogicHoneyCombLatX = new G4LogicalVolume(fSolidHoneyCombLatX,hcMat,"HoneyCombCore"); 

  fLogicHoneyCombSkinLatX = new G4LogicalVolume(fSolidHoneyCombSkinLatX,cfMat,"HoneyCombSkin"); 

  fLogicHoneyCombBottom = new G4LogicalVolume(fSolidHoneyCombBottom,hcMat,"HoneyCombCore"); 

  fLogicHoneyCombSkinBottom = new G4LogicalVolume(fSolidHoneyCombSkinBottom,cfMat,"HoneyCombSkin"); 

  if(fNbOfCrystalLayers>0){
    fLogicScintActiveLayer->SetSensitiveDetector(caloSD);
    fLogicCrystalActiveLayer->SetSensitiveDetector(caloSD);
  }else{
    fLogicActiveLayer->SetSensitiveDetector(caloSD);
  }
  fLogicVetoBottom->SetSensitiveDetector(vetoSD);
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
      
    fPhysiScintLayer = new G4PVReplica("LayerScint", 
				       fLogicScintLayer,            
				       fPhysiScintBox,            
				       kZAxis,                 
				       fNbOfScintLayers,        
				       fCaloLayer_Z);   
    
    fPhysiScintCFSupp = new G4PVPlacement(0,                   
					  G4ThreeVector(fPhysiCFSupp_X,fPhysiCFSupp_Y,fPhysiCFSupp_Z),     
					  "CFSupport",              
					  fLogicCFSupp,       
					  fPhysiScintLayer,        
					  false,               
					  0,true); 

    fPhysiScintActiveLayer = new G4PVPlacement(0,                   
					       G4ThreeVector(fPhysiActiveLayer_X,fPhysiActiveLayer_Y,fPhysiActiveLayer_Z),     
					       "ActiveLayerScint",              
					       fLogicScintActiveLayer,       
					       fPhysiScintLayer,        
					       false,               
					       0,true); 

    fPhysiScintCFFrontP = new G4PVPlacement(0,                   
					    G4ThreeVector(fPhysiCFFront_X,fPhysiCFFront_Y,fPhysiCFFront_Z),     
					    "CFSupportFP",              
					    fLogicCFFront,       
					    fPhysiScintLayer,        
					    false,               
					    0,true); 

    fPhysiScintCFFrontM = new G4PVPlacement(0,                   
					    G4ThreeVector(fPhysiCFFront_X,-fPhysiCFFront_Y,fPhysiCFFront_Z),     
					    "CFSupportFM",              
					    fLogicCFFront,       
					    fPhysiScintLayer,        
					    false,               
					    0,true); 

    fPhysiScintCFLatP = new G4PVPlacement(0,                   
					  G4ThreeVector(fPhysiCFLat_X,fPhysiCFLat_Y,fPhysiCFLat_Z),     
					  "CFSupportLP",              
					  fLogicCFLat,       
					  fPhysiScintLayer,        
					  false,               
					  0,true); 

    fPhysiScintCFLatM = new G4PVPlacement(0,                   
					  G4ThreeVector(-fPhysiCFLat_X,fPhysiCFLat_Y,fPhysiCFLat_Z),     
					  "CFSupportLM",              
					  fLogicCFLat,       
					  fPhysiScintLayer,        
					  false,               
					  0,true); 

    fPhysiScintPoronLatP = new G4PVPlacement(0,                   
					     G4ThreeVector(fPhysiPoronLat_X,fPhysiPoronLat_Y,fPhysiPoronLat_Z),     
					     "PoronLP",              
					     fLogicPoronLat,       
					     fPhysiScintLayer,        
					     false,               
					     0,true); 

    fPhysiScintPoronLatM = new G4PVPlacement(0,                   
					     G4ThreeVector(-fPhysiPoronLat_X,fPhysiPoronLat_Y,fPhysiPoronLat_Z),     
					     "PoronLM",              
					     fLogicPoronLat,       
					     fPhysiScintLayer,        
					     false,               
					     0,true); 

    fPhysiCrystalBox = new G4PVPlacement(0,                        
				      G4ThreeVector(fPhysiCrystalBox_X,fPhysiCrystalBox_Y,fPhysiCrystalBox_Z),   
				      "CalorimeterCrystal",               
				      fLogicCrystalBox,                
				      fPhysiCaloBox,                
				      false,                       
				      0,true);       
  

    fPhysiCrystalLayer = new G4PVReplica("LayerCrystal", 
					 fLogicCrystalLayer,            
					 fPhysiCrystalBox,            
					 kZAxis,                 
					 fNbOfCrystalLayers,        
					 fCaloLayer_Z);   
    
    fPhysiCrystalCFSupp = new G4PVPlacement(0,                   
					    G4ThreeVector(fPhysiCFSupp_X,fPhysiCFSupp_Y,fPhysiCFSupp_Z),     
					    "CFSupport",              
					    fLogicCFSupp,       
					    fPhysiCrystalLayer,        
					    false,               
					    0,true); 

    fPhysiCrystalActiveLayer = new G4PVPlacement(0,                   
						 G4ThreeVector(fPhysiActiveLayer_X,fPhysiActiveLayer_Y,fPhysiActiveLayer_Z),     
						 "ActiveLayerCrystal",              
						 fLogicCrystalActiveLayer,       
						 fPhysiCrystalLayer,        
						 false,               
						 0,true); 

    fPhysiCrystalCFFrontP = new G4PVPlacement(0,                   
					      G4ThreeVector(fPhysiCFFront_X,fPhysiCFFront_Y,fPhysiCFFront_Z),     
					      "CFSupportFP",              
					      fLogicCFFront,       
					      fPhysiCrystalLayer,        
					      false,               
					      0,true); 

    fPhysiCrystalCFFrontM = new G4PVPlacement(0,                   
					      G4ThreeVector(fPhysiCFFront_X,-fPhysiCFFront_Y,fPhysiCFFront_Z),     
					      "CFSupportFM",              
					      fLogicCFFront,       
					      fPhysiCrystalLayer,        
					      false,               
					      0,true); 

    fPhysiCrystalCFLatP = new G4PVPlacement(0,                   
					    G4ThreeVector(fPhysiCFLat_X,fPhysiCFLat_Y,fPhysiCFLat_Z),     
					    "CFSupportLP",              
					    fLogicCFLat,       
					    fPhysiCrystalLayer,        
					    false,               
					    0,true); 

    fPhysiCrystalCFLatM = new G4PVPlacement(0,                   
					    G4ThreeVector(-fPhysiCFLat_X,fPhysiCFLat_Y,fPhysiCFLat_Z),     
					    "CFSupportLM",              
					    fLogicCFLat,       
					    fPhysiCrystalLayer,        
					    false,               
					    0,true); 

    fPhysiCrystalPoronLatP = new G4PVPlacement(0,                   
					       G4ThreeVector(fPhysiPoronLat_X,fPhysiPoronLat_Y,fPhysiPoronLat_Z),     
					       "PoronLP",              
					       fLogicPoronLat,       
					       fPhysiCrystalLayer,        
					       false,               
					       0,true); 

    fPhysiCrystalPoronLatM = new G4PVPlacement(0,                   
					       G4ThreeVector(-fPhysiPoronLat_X,fPhysiPoronLat_Y,fPhysiPoronLat_Z),     
					       "PoronLM",              
					       fLogicPoronLat,       
					       fPhysiCrystalLayer,        
					       false,               
					       0,true); 


  }else{
    fPhysiLayer = new G4PVReplica("Layer", 
				  fLogicLayer,            
				  fPhysiCaloBox,            
				  kZAxis,                 
				  fNbOfCaloLayers,        
				  fCaloLayer_Z);   
    
    fPhysiCFSupp = new G4PVPlacement(0,                   
				     G4ThreeVector(fPhysiCFSupp_X,fPhysiCFSupp_Y,fPhysiCFSupp_Z),     
				     "CFSupport",              
				     fLogicCFSupp,       
				     fPhysiLayer,        
				     false,               
				     0,true); 
  
    fPhysiActiveLayer = new G4PVPlacement(0,                   
					  G4ThreeVector(fPhysiActiveLayer_X,fPhysiActiveLayer_Y,fPhysiActiveLayer_Z),     
					  "ActiveLayer",              
					  fLogicActiveLayer,       
					  fPhysiLayer,        
					  false,               
					  0,true); 

    fPhysiCFFrontP = new G4PVPlacement(0,                   
				       G4ThreeVector(fPhysiCFFront_X,fPhysiCFFront_Y,fPhysiCFFront_Z),     
				       "CFSupportFP",              
				       fLogicCFFront,       
				       fPhysiLayer,        
				       false,               
				       0,true); 

    fPhysiCFFrontM = new G4PVPlacement(0,                   
				       G4ThreeVector(fPhysiCFFront_X,-fPhysiCFFront_Y,fPhysiCFFront_Z),     
				       "CFSupportFM",              
				       fLogicCFFront,       
				       fPhysiLayer,        
				       false,               
				       0,true); 

    fPhysiCFLatP = new G4PVPlacement(0,                   
				     G4ThreeVector(fPhysiCFLat_X,fPhysiCFLat_Y,fPhysiCFLat_Z),     
				     "CFSupportLP",              
				     fLogicCFLat,       
				     fPhysiLayer,        
				     false,               
				     0,true); 

    fPhysiCFLatM = new G4PVPlacement(0,                   
				     G4ThreeVector(-fPhysiCFLat_X,fPhysiCFLat_Y,fPhysiCFLat_Z),     
				     "CFSupportLM",              
				     fLogicCFLat,       
				     fPhysiLayer,        
				     false,               
				     0,true); 

    fPhysiPoronLatP = new G4PVPlacement(0,                   
					G4ThreeVector(fPhysiPoronLat_X,fPhysiPoronLat_Y,fPhysiPoronLat_Z),     
					"PoronLP",              
					fLogicPoronLat,       
					fPhysiLayer,        
					false,               
					0,true); 

    fPhysiPoronLatM = new G4PVPlacement(0,                   
					G4ThreeVector(-fPhysiPoronLat_X,fPhysiPoronLat_Y,fPhysiPoronLat_Z),     
					"PoronLM",              
					fLogicPoronLat,       
					fPhysiLayer,        
					false,               
					0,true); 

  }

  fPhysiPoronLatXPInt = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiPoronLatXInt_X,fPhysiPoronLatXInt_Y,fPhysiPoronLatXInt_Z),     
					  "PoronLXPI",              
					  fLogicPoronLatX,       
					  motherVolume,        
					  false,               
					  0,true); 
  
  fPhysiPoronLatXMInt = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiPoronLatXInt_X,-fPhysiPoronLatXInt_Y,fPhysiPoronLatXInt_Z),     
					  "PoronLXMI",              
					  fLogicPoronLatX,       
					  motherVolume,        
					  false,               
					  0,true); 

  fPhysiPoronLatXPExt = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiPoronLatXExt_X,fPhysiPoronLatXExt_Y,fPhysiPoronLatXExt_Z),     
					  "PoronLXPE",              
					  fLogicPoronLatX,       
					  motherVolume,        
					  false,               
					  0,true); 
  
  fPhysiPoronLatXMExt = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiPoronLatXExt_X,-fPhysiPoronLatXExt_Y,fPhysiPoronLatXExt_Z),     
					  "PoronLXME",              
					  fLogicPoronLatX,       
					  motherVolume,        
					  false,               
					  0,true); 
  

  fPhysiPoronLatYPInt = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiPoronLatYInt_X,fPhysiPoronLatYInt_Y,fPhysiPoronLatYInt_Z),     
					  "PoronLYPI",              
					  fLogicPoronLatY,       
					  motherVolume,        
					  false,               
					  0,true); 
  
  fPhysiPoronLatYMInt = new G4PVPlacement(0,
					  G4ThreeVector(-fPhysiPoronLatYInt_X,fPhysiPoronLatYInt_Y,fPhysiPoronLatYInt_Z),     
					  "PoronLYMI",              
					  fLogicPoronLatY,       
					  motherVolume,        
					  false,               
					  0,true); 

  fPhysiPoronLatYPExt = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiPoronLatYExt_X,fPhysiPoronLatYExt_Y,fPhysiPoronLatYExt_Z),     
					  "PoronLYPE",              
					  fLogicPoronLatY,       
					  motherVolume,        
					  false,               
					  0,true); 
  
  fPhysiPoronLatYMExt = new G4PVPlacement(0,
					  G4ThreeVector(-fPhysiPoronLatYExt_X,fPhysiPoronLatYExt_Y,fPhysiPoronLatYExt_Z),     
					  "PoronLYME",              
					  fLogicPoronLatY,       
					  motherVolume,        
					  false,               
					  0,true); 
  
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
  
  fPhysiPoronBottomU = new G4PVPlacement(0,
					 G4ThreeVector(fPhysiPoronBottomU_X,fPhysiPoronBottomU_Y,fPhysiPoronBottomU_Z),     
					 "PoronBU",              
					 fLogicPoronBottom,       
					 motherVolume,        
					 false,               
					 0,true); 
  
  fPhysiPoronBottomD = new G4PVPlacement(0,
					 G4ThreeVector(fPhysiPoronBottomD_X,fPhysiPoronBottomD_Y,fPhysiPoronBottomD_Z),     
					 "PoronBD",              
					 fLogicPoronBottom,       
					 motherVolume,        
					 false,               
					 0,true); 
  
  fPhysiPoronBottomLatStickP = new G4PVPlacement(0,
						 G4ThreeVector(fPhysiPoronBottomLatStick_X,fPhysiPoronBottomLatStick_Y,fPhysiPoronBottomLatStick_Z),
						 "PoronBLP",
						 fLogicPoronBottomLatStick,
						 motherVolume,        
						 false,               
						 0,true); 

  fPhysiPoronBottomLatStickM = new G4PVPlacement(0,
						 G4ThreeVector(-fPhysiPoronBottomLatStick_X,fPhysiPoronBottomLatStick_Y,fPhysiPoronBottomLatStick_Z),
						 "PoronBLM",
						 fLogicPoronBottomLatStick,
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

  fPhysiCFBottomLatP = new G4PVPlacement(0,
					 G4ThreeVector(fPhysiCFBottomLat_X,fPhysiCFBottomLat_Y,fPhysiCFBottomLat_Z),
					 "CFSupportBLP",
					 fLogicCFBottomLat,
					 motherVolume,
					 false,
					 0,true);

  fPhysiCFBottomLatM = new G4PVPlacement(0,
					 G4ThreeVector(-fPhysiCFBottomLat_X,fPhysiCFBottomLat_Y,fPhysiCFBottomLat_Z),
					 "CFSupportBLM",
					 fLogicCFBottomLat,
					 motherVolume,
					 false,
					 0,true);

  fPhysiCFBottomFrontP = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiCFBottomFront_X,fPhysiCFBottomFront_Y,fPhysiCFBottomFront_Z),
					   "CFSupportBFP",
					   fLogicCFBottomFront,
					   motherVolume,
					   false,
					   0,true);

  fPhysiCFBottomFrontM = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiCFBottomFront_X,-fPhysiCFBottomFront_Y,fPhysiCFBottomFront_Z),
					   "CFSupportBFM",
					   fLogicCFBottomFront,
					   motherVolume,
					   false,
					   0,true);

  fPhysiHoneyCombSkinLatYPInt = new G4PVPlacement(0,
						  G4ThreeVector(fPhysiHoneyCombSkinLatYInt_X,fPhysiHoneyCombSkinLatYInt_Y,fPhysiHoneyCombSkinLatYInt_Z),
						  "HoneyCombSkinLYPI",
						  fLogicHoneyCombSkinLatY,
						  motherVolume,
						  false,
						  0,true);
  
  fPhysiHoneyCombLatYP = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiHoneyCombLatY_X,fPhysiHoneyCombLatY_Y,fPhysiHoneyCombLatY_Z),
					   "HoneyCombCoreLYP",
					   fLogicHoneyCombLatY,
					   motherVolume,
					   false,
					   0,true);

  fPhysiHoneyCombSkinLatYPExt = new G4PVPlacement(0,
						  G4ThreeVector(fPhysiHoneyCombSkinLatYExt_X,fPhysiHoneyCombSkinLatYExt_Y,fPhysiHoneyCombSkinLatYExt_Z),
						  "HoneyCombSkinLYPE",
						  fLogicHoneyCombSkinLatY,
						  motherVolume,
						  false,
						  0,true);
  
  fPhysiHoneyCombSkinLatYMInt = new G4PVPlacement(0,
						  G4ThreeVector(-fPhysiHoneyCombSkinLatYInt_X,fPhysiHoneyCombSkinLatYInt_Y,fPhysiHoneyCombSkinLatYInt_Z),
						  "HoneyCombSkinLYMI",
						  fLogicHoneyCombSkinLatY,
						  motherVolume,
						  false,
						  0,true);

  fPhysiHoneyCombLatYM = new G4PVPlacement(0,
					   G4ThreeVector(-fPhysiHoneyCombLatY_X,fPhysiHoneyCombLatY_Y,fPhysiHoneyCombLatY_Z),
					   "HoneyCombCoreLYM",
					   fLogicHoneyCombLatY,
					   motherVolume,
					   false,
					   0,true);
  
  fPhysiHoneyCombSkinLatYMExt = new G4PVPlacement(0,
						  G4ThreeVector(-fPhysiHoneyCombSkinLatYExt_X,fPhysiHoneyCombSkinLatYExt_Y,fPhysiHoneyCombSkinLatYExt_Z),
						  "HoneyCombSkinLYME",
						  fLogicHoneyCombSkinLatY,
						  motherVolume,
						  false,
						  0,true);

  fPhysiHoneyCombSkinLatXPInt = new G4PVPlacement(0,
						  G4ThreeVector(fPhysiHoneyCombSkinLatXInt_X,fPhysiHoneyCombSkinLatXInt_Y,fPhysiHoneyCombSkinLatXInt_Z),
						  "HoneyCombSkinLXPI",
						  fLogicHoneyCombSkinLatX,
						  motherVolume,
						  false,
						  0,true);
  
  fPhysiHoneyCombLatXP = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiHoneyCombLatX_X,fPhysiHoneyCombLatX_Y,fPhysiHoneyCombLatX_Z),
					   "HoneyCombCoreLXP",
					   fLogicHoneyCombLatX,
					   motherVolume,
					   false,
					   0,true);

  fPhysiHoneyCombSkinLatXPExt = new G4PVPlacement(0,
						  G4ThreeVector(fPhysiHoneyCombSkinLatXExt_X,fPhysiHoneyCombSkinLatXExt_Y,fPhysiHoneyCombSkinLatXExt_Z),
						  "HoneyCombSkinLXPE",
						  fLogicHoneyCombSkinLatX,
						  motherVolume,
						  false,
						  0,true);
  
  fPhysiHoneyCombSkinLatXMInt = new G4PVPlacement(0,
						  G4ThreeVector(fPhysiHoneyCombSkinLatXInt_X,-fPhysiHoneyCombSkinLatXInt_Y,fPhysiHoneyCombSkinLatXInt_Z),
						  "HoneyCombSkinLXMI",
						  fLogicHoneyCombSkinLatX,
						  motherVolume,
						  false,
						  0,true);

  fPhysiHoneyCombLatXM = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiHoneyCombLatX_X,-fPhysiHoneyCombLatX_Y,fPhysiHoneyCombLatX_Z),
					   "HoneyCombCoreLXM",
					   fLogicHoneyCombLatX,
					   motherVolume,
					   false,
					   0,true);
  
  fPhysiHoneyCombSkinLatXMExt = new G4PVPlacement(0,
						  G4ThreeVector(fPhysiHoneyCombSkinLatXExt_X,-fPhysiHoneyCombSkinLatXExt_Y,fPhysiHoneyCombSkinLatXExt_Z),
						  "HoneyCombSkinLXME",
						  fLogicHoneyCombSkinLatX,
						  motherVolume,
						  false,
						  0,true);

  fPhysiHoneyCombSkinBottomU = new G4PVPlacement(0,
						 G4ThreeVector(fPhysiHoneyCombSkinBottomU_X,fPhysiHoneyCombSkinBottomU_Y,fPhysiHoneyCombSkinBottomU_Z),
						 "HoneyCombSkinBU",
						 fLogicHoneyCombSkinBottom,
						 motherVolume,
						 false,
						 0,true);
  
  fPhysiHoneyCombBottom = new G4PVPlacement(0,
					    G4ThreeVector(fPhysiHoneyCombBottom_X,fPhysiHoneyCombBottom_Y,fPhysiHoneyCombBottom_Z),	
					    "HoneyCombCoreB",
					    fLogicHoneyCombBottom,
					    motherVolume,
					    false,
					    0,true);
  
  fPhysiHoneyCombSkinBottomD = new G4PVPlacement(0,
						 G4ThreeVector(fPhysiHoneyCombSkinBottomD_X,fPhysiHoneyCombSkinBottomD_Y,fPhysiHoneyCombSkinBottomD_Z),
						 "HoneyCombSkinBD",
						 fLogicHoneyCombSkinBottom,
						 motherVolume,
						 false,

						 0,true);


  fPhysiCFSuppTop = new G4PVPlacement(0,
				      G4ThreeVector(fPhysiCFSuppTop_X,fPhysiCFSuppTop_Y,fPhysiCFSuppTop_Z),
				      "CFSupportT",
				      fLogicCFSuppTop,
				      motherVolume,
				      false,
				      0,true);

  //Visualization Attribute
  G4VisAttributes* attCyan = new G4VisAttributes(G4Colour::Cyan());
  attCyan->SetVisibility(true);
  attCyan->SetForceAuxEdgeVisible(true);
  fLogicHoneyCombLatY->SetVisAttributes(attCyan);
  fLogicHoneyCombLatX->SetVisAttributes(attCyan);
  fLogicHoneyCombBottom->SetVisAttributes(attCyan);

  G4VisAttributes* attMagenta = new G4VisAttributes(G4Colour::Magenta());
  attMagenta->SetVisibility(true);
  attMagenta->SetForceAuxEdgeVisible(true);
  if(fNbOfCrystalLayers>0){
    fLogicScintActiveLayer->SetVisAttributes(attMagenta);
  }else{
    fLogicActiveLayer->SetVisAttributes(attMagenta);
  }
  G4VisAttributes* attYellow = new G4VisAttributes(G4Colour::Yellow());
  attYellow->SetVisibility(true);
  attYellow->SetForceAuxEdgeVisible(true);
  fLogicPoronLat->SetVisAttributes(attYellow);
  fLogicPoronLatX->SetVisAttributes(attYellow);
  fLogicPoronLatY->SetVisAttributes(attYellow);
  fLogicPoronBottom->SetVisAttributes(attYellow);
  fLogicPoronBottomLatStick->SetVisAttributes(attYellow);

  G4VisAttributes* attBrown = new G4VisAttributes(G4Colour::Brown());
  attBrown->SetVisibility(true);
  attBrown->SetForceAuxEdgeVisible(true);
  fLogicCFLat->SetVisAttributes(attBrown);
  fLogicCFFront->SetVisAttributes(attBrown);
  fLogicCFSupp->SetVisAttributes(attBrown);
  fLogicCFSuppTop->SetVisAttributes(attBrown);
  fLogicCFBottomLat->SetVisAttributes(attBrown);
  fLogicCFBottomFront->SetVisAttributes(attBrown);
  fLogicHoneyCombSkinLatY->SetVisAttributes(attBrown);
  fLogicHoneyCombSkinLatX->SetVisAttributes(attBrown);
  fLogicHoneyCombSkinBottom->SetVisAttributes(attBrown);


  G4VisAttributes* attGreen = new G4VisAttributes(G4Colour::Green());
  attGreen->SetVisibility(true);
  attGreen->SetForceAuxEdgeVisible(true);
  fLogicVetoBottom->SetVisAttributes(attGreen);
  fLogicVetoLatY->SetVisAttributes(attGreen);
  fLogicVetoLatX->SetVisAttributes(attGreen);

  G4VisAttributes* attRed = new G4VisAttributes(G4Colour::Red());
  attRed->SetVisibility(true);
  attRed->SetForceAuxEdgeVisible(true);
  if(fNbOfCrystalLayers>0){
    fLogicCrystalActiveLayer->SetVisAttributes(attRed);
  }


  G4VisAttributes * attInvisible = new G4VisAttributes();
  attInvisible->SetVisibility(false);
  attInvisible->SetForceAuxEdgeVisible(false);
  fLogicCaloBox->SetVisAttributes(attInvisible); 
  if(fNbOfCrystalLayers>0){
    fLogicScintBox->SetVisAttributes(attInvisible); 
    fLogicCrystalBox->SetVisAttributes(attInvisible); 
    fLogicScintLayer->SetVisAttributes(attInvisible); 
    fLogicCrystalLayer->SetVisAttributes(attInvisible); 
  }else{
    fLogicLayer->SetVisAttributes(attInvisible); 
  }
}


