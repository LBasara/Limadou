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

#include "CalorimeterConstructionConfig1.hh"

#include "G4NistManager.hh"
#include "HEPDSWMaterial.hh"
#include "G4Box.hh"
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

CalorimeterConstructionConfig1::CalorimeterConstructionConfig1()
  :fSolidCaloBox(0),fSolidLayer(0),fSolidActiveLayer(0),
   fSolidCFPlane(0),fSolidCFBS(0),fSolidCFSS(0),fSolidCFYPStick(0),fSolidCFYSBStick(0),fSolidCFYBBStick(0),
   fSolidPoronLayer(0),fSolidPoronStick(0),fSolidVetoBottom(0),fSolidCFStickLat(0),fSolidCFStickFront(0),
  fSolidHCBMSkin(0),fSolidHCBLSkin(0),fSolidHCBMCore(0),fSolidHCBLCore(0),fSolidVetoLatY(0),fSolidPoronVetoLatY(0),
  fSolidVetoLatX(0),fSolidPoronVetoLatX(0),fSolidHCSkinLatX(0),fSolidHCCoreLatX(0),fSolidHCSkinLatY(0),fSolidHCCoreLatY(0),
  fLogicCaloBox(0),fLogicLayer(0),fLogicActiveLayer(0),
  fLogicCFPlane(0),fLogicCFBS(0),fLogicCFSS(0),fLogicCFYPStick(0),fLogicCFYSBStick(0),fLogicCFYBBStick(0),
  fLogicPoronLayer(0),fLogicPoronStick(0),fLogicVetoBottom(0),fLogicCFStickLat(0),fLogicCFStickFront(0),
  fLogicHCBMSkin(0),fLogicHCBLSkin(0),fLogicHCBMCore(0),fLogicHCBLCore(0),fLogicVetoLatY(0),fLogicPoronVetoLatY(0),
  fLogicVetoLatX(0),fLogicPoronVetoLatX(0),fLogicHCSkinLatX(0),fLogicHCCoreLatX(0),fLogicHCSkinLatY(0),fLogicHCCoreLatY(0),
  fPhysiCaloBox(0),fPhysiLayer(0),fPhysiActiveLayer(0),
  fPhysiCFPlaneU(0),fPhysiCFPlaneD(0),fPhysiPoronLayer(0),fPhysiPoronStickR(0),fPhysiPoronStickL(0),
  fPhysiCFBSR(0),fPhysiCFBSL(0),fPhysiCFSSR(0),fPhysiCFSSL(0),fPhysiCFYPStickUP(0),fPhysiCFYPStickUM(0),
  fPhysiCFYBBStickP(0),fPhysiCFYBBStickM(0),fPhysiCFYPStickDP(0),
  fPhysiCFYPStickDM(0),fPhysiCFYSBStickP(0),fPhysiCFYSBStickM(0),fPhysiVetoBottom(0),fPhysiPoronVetoBottom(0),
  fPhysiPoronVetoLatR(0),fPhysiPoronVetoLatL(0),fPhysiCFStickLatR(0),fPhysiCFStickLatL(0),fPhysiCFStickFrontP(0),
  fPhysiCFStickFrontM(0),fPhysiHCMSkinU(0),fPhysiHCMCore(0),fPhysiHCMSkinD(0),fPhysiHCLSkinUR(0),
  fPhysiHCLCoreR(0),fPhysiHCLSkinDR(0),fPhysiHCLSkinUL(0),fPhysiHCLCoreL(0),fPhysiHCLSkinDL(0),
  fPhysiPoronVetoLatYPInt(0),fPhysiVetoLatYP(0),fPhysiPoronVetoLatYPExt(0),fPhysiPoronVetoLatYMInt(0),
  fPhysiVetoLatYM(0),fPhysiPoronVetoLatYMExt(0),fPhysiPoronVetoLatXPInt(0),fPhysiVetoLatXP(0),
  fPhysiPoronVetoLatXPExt(0),fPhysiPoronVetoLatXMInt(0),fPhysiVetoLatXM(0),fPhysiPoronVetoLatXMExt(0),
  fPhysiHCSkinLatXPInt(0),fPhysiHCCoreLatXP(0),fPhysiHCSkinLatXPExt(0),
  fPhysiHCSkinLatXMInt(0),fPhysiHCCoreLatXM(0),fPhysiHCSkinLatXMExt(0),
  fPhysiHCSkinLatYPInt(0),fPhysiHCCoreLatYP(0),fPhysiHCSkinLatYPExt(0),
  fPhysiHCSkinLatYMInt(0),fPhysiHCCoreLatYM(0),fPhysiHCSkinLatYMExt(0)
{
  pMaterial           = new HEPDSWMaterial();

  fNbOfCaloLayers = 20;
  fCaloLayer_X    = 22*cm;
  fCaloLayer_Y    = 22.7*cm;
  fCaloLayer_Z    = 15.18*mm;

  fCalo_X = 22*cm;
  fCalo_Y = 22.7*cm;
  fCalo_Z = fCaloLayer_Z*fNbOfCaloLayers;

  fCFPlane_X = 20*cm;
  fCFPlane_Y = 22.7*cm;
  fCFPlane_Z = 1*mm;

  fCFBS_X = 20*cm;
  fCFBS_Y = 1*cm;
  fCFBS_Z = 1*cm;

  fCFSS_X = 20*cm;
  fCFSS_Y = 1*cm;
  fCFSS_Z = 3.18*mm;

  fCFYPStick_X = 1*cm;
  fCFYPStick_Y = 10*cm;
  fCFYPStcik_Z = 1*mm;

  fCFYSBStick_X = 1*cm;
  fCFYSBStick_Y = 10*cm;
  fCFYSBStick_Z = 3.18*mm;

  fCFYBBStick_X = 1*cm;
  fCFYBBStick_Y = 10*cm;
  fCFYBBStick_Z = 1*cm;

  fPoronL_X = 20*cm;
  fPoronL_Y = 20.7*cm;
  fPoronL_Z = 3.18*mm;
  
  fPoronS_X = 20*cm;
  fPoronS_Y = 3.18*mm;
  fPoronS_Z = 1*cm;

  fActiveLayer_X = 20*cm;
  fActiveLayer_Y = 20*cm;
  fActiveLayer_Z = 1*cm;

  fVetoBottom_X = 20*cm;
  fVetoBottom_Y = 20*cm;
  fVetoBottom_Z = 1*cm;

  fCFStickLat_X = 20*cm;
  fCFStickLat_Y = 10*mm;
  fCFStickLat_Z = 13.18*mm;

  fCFStickFront_X = 10*mm;
  fCFStickFront_Y = 100*mm;
  fCFStickFront_Z = 13.18*mm;

  fHCBMSkin_X = 252.72*mm;
  fHCBMSkin_Y = 336*mm;
  fHCBMSkin_Z = 1*mm;

  fHCBLSkin_X = 43.64*mm;
  fHCBLSkin_Y = 220*mm;
  fHCBLSkin_Z = 1*mm;

  fHCBMCore_X = 252.72*mm;
  fHCBMCore_Y = 336*mm;
  fHCBMCore_Z = 10*mm;

  fHCBLCore_X = 43.64*mm;
  fHCBLCore_Y = 220*mm;
  fHCBLCore_Z = 10*mm;

  fVetoLatY_X = 220*mm;
  fVetoLatY_Y = 10*mm;
  fVetoLatY_Z = 316.78*mm;
  
  fPoronVetoLatY_X = 220*mm;
  fPoronVetoLatY_Y = 3.18*mm;
  fPoronVetoLatY_Z = 316.78*mm;
  
  fVetoLatX_X = 10*mm;
  fVetoLatX_Y = 259.72*mm;
  fVetoLatX_Z = 316.78*mm;
  
  fPoronVetoLatX_X = 3.18*mm;
  fPoronVetoLatX_Y = 259.72*mm;
  fPoronVetoLatX_Z = 316.78*mm;
  
  fHCSkinLatX_X = 1*mm;
  fHCSkinLatX_Y = 259.72*mm;
  fHCSkinLatX_Z = 316.78*mm;

  fHCCoreLatX_X = 10*mm;
  fHCCoreLatX_Y = 259.72*mm;
  fHCCoreLatX_Z = 316.78*mm;

  fHCSkinLatY_X = 252.72*mm;
  fHCSkinLatY_Y = 1*mm;
  fHCSkinLatY_Z = 316.78*mm;

  fHCCoreLatY_X = 252.72*mm;
  fHCCoreLatY_Y = 10*mm;
  fHCCoreLatY_Z = 316.78*mm;


  ComputeObjectsPositioning();


  // materials
  vetoMaterial      = "Scintillator";
  caloMaterial      = "Scintillator";

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

CalorimeterConstructionConfig1::~CalorimeterConstructionConfig1()
{
  if (pMaterial) delete pMaterial;
}


void CalorimeterConstructionConfig1::ComputeObjectsPositioning(){

  fPhysiCFPlaneU_X = 0;
  fPhysiCFPlaneU_Y = 0;
  fPhysiCFPlaneU_Z = fCaloLayer_Z/2.-fCFPlane_Z/2.;

  fPhysiActiveLayer_X = 0;
  fPhysiActiveLayer_Y = 0;
  fPhysiActiveLayer_Z = fCaloLayer_Z/2.-fCFPlane_Z-fActiveLayer_Z/2.;

  fPhysiCFPlaneD_X = 0;
  fPhysiCFPlaneD_Y = 0;
  fPhysiCFPlaneD_Z = fCaloLayer_Z/2.-fCFPlane_Z-fActiveLayer_Z-fCFPlane_Z/2.;

  fPhysiPoronLayer_X = 0;
  fPhysiPoronLayer_Y = 0;
  fPhysiPoronLayer_Z = fCaloLayer_Z/2.-fCFPlane_Z-fActiveLayer_Z-fCFPlane_Z-fPoronL_Z/2.;

  fPhysiPoronStick_X = 0;
  fPhysiPoronStick_Y = fActiveLayer_Y/2.+fPoronS_Y/2.;
  fPhysiPoronStick_Z = fCaloLayer_Z/2.-fCFPlane_Z-fActiveLayer_Z/2.;

  fPhysiCFBS_X = 0;
  fPhysiCFBS_Y = fPoronL_Y/2.+fCFBS_Y/2.;
  fPhysiCFBS_Z = fCaloLayer_Z/2.-fCFPlane_Z-fActiveLayer_Z/2.;

  fPhysiCFSS_X = 0;
  fPhysiCFSS_Y = fPoronL_Y/2.+fCFSS_Y/2.;
  fPhysiCFSS_Z = fCaloLayer_Z/2.-fCFPlane_Z-fActiveLayer_Z-fCFPlane_Z-fPoronL_Z/2.;

  fPhysiCFYPStickU_X = fCFPlane_X/2.+fCFYPStick_X/2.;
  fPhysiCFYPStickU_Y = 0;
  fPhysiCFYPStickU_Z = fCaloLayer_Z/2.-fCFPlane_Z/2.;

  fPhysiCFYBBStick_X = fCFPlane_X/2.+fCFYBBStick_X/2.;
  fPhysiCFYBBStick_Y = 0;
  fPhysiCFYBBStick_Z = fCaloLayer_Z/2.-fCFPlane_Z-fActiveLayer_Z/2.;

  fPhysiCFYPStickD_X = fCFPlane_X/2.+fCFYPStick_X/2.;
  fPhysiCFYPStickD_Y = 0;
  fPhysiCFYPStickD_Z = fCaloLayer_Z/2.-fCFPlane_Z-fActiveLayer_Z-fCFPlane_Z/2.;

  fPhysiCFYSBStick_X = fCFPlane_X/2.+fCFYSBStick_X/2.;
  fPhysiCFYSBStick_Y = 0;
  fPhysiCFYSBStick_Z = fCaloLayer_Z/2.-fCFPlane_Z-fActiveLayer_Z-fCFPlane_Z-fPoronL_Z/2.;

  fPhysiVetoBottom_X = 0;
  fPhysiVetoBottom_Y = 0;
  fPhysiVetoBottom_Z = -fCalo_Z/2.-fVetoBottom_Z/2.;

  fPhysiPoronVetoBottom_X = 0;
  fPhysiPoronVetoBottom_Y = 0;
  fPhysiPoronVetoBottom_Z = -fCalo_Z/2.-fVetoBottom_Z-fPoronL_Z/2.;

  fPhysiPoronVetoLat_X = 0;
  fPhysiPoronVetoLat_Y = fVetoBottom_Y/2.+fPoronVetoLatY_Y/2.;
  fPhysiPoronVetoLat_Z = -fCalo_Z/2.-fVetoBottom_Z/2.;

  fPhysiCFStickLat_X = 0;
  fPhysiCFStickLat_Y = fPoronL_Y/2.+fCFStickLat_Y/2.;
  fPhysiCFStickLat_Z = -fCalo_Z/2.-fVetoBottom_Z/2.-fPoronL_Z/2.;

  fPhysiCFStickFront_X = fPoronL_X/2.+fCFStickFront_X/2.;
  fPhysiCFStickFront_Y = 0;
  fPhysiCFStickFront_Z = -fCalo_Z/2.-fVetoBottom_Z/2.-fPoronL_Z/2.;

  fPhysiHCMSkinU_X = 0;
  fPhysiHCMSkinU_Y = 0;
  fPhysiHCMSkinU_Z = -fCalo_Z/2.-fVetoBottom_Z-fPoronL_Z-fHCBMSkin_Z/2.;

  fPhysiHCMCore_X = 0;
  fPhysiHCMCore_Y = 0;
  fPhysiHCMCore_Z = -fCalo_Z/2.-fVetoBottom_Z-fPoronL_Z-fHCBMSkin_Z-fHCBMCore_Z/2.;

  fPhysiHCMSkinD_X = 0;
  fPhysiHCMSkinD_Y = 0;
  fPhysiHCMSkinD_Z = -fCalo_Z/2.-fVetoBottom_Z-fPoronL_Z-fHCBMSkin_Z-fHCBMCore_Z-fHCBMSkin_Z/2.;

  fPhysiHCLSkinU_X = fHCBMSkin_X/2.+fHCBLSkin_X/2.;
  fPhysiHCLSkinU_Y = 0;
  fPhysiHCLSkinU_Z = -fCalo_Z/2.-fVetoBottom_Z-fPoronL_Z-fHCBMSkin_Z/2.;

  fPhysiHCLCore_X = fHCBMSkin_X/2.+fHCBLSkin_X/2.;
  fPhysiHCLCore_Y = 0;
  fPhysiHCLCore_Z = -fCalo_Z/2.-fVetoBottom_Z-fPoronL_Z-fHCBMSkin_Z-fHCBMCore_Z/2.;

  fPhysiHCLSkinD_X = fHCBMSkin_X/2.+fHCBLSkin_X/2.;
  fPhysiHCLSkinD_Y = 0;
  fPhysiHCLSkinD_Z = -fCalo_Z/2.-fVetoBottom_Z-fPoronL_Z-fHCBMSkin_Z-fHCBMCore_Z-fHCBMSkin_Z/2.;

  fPhysiPoronVetoLatYInt_X = 0;
  fPhysiPoronVetoLatYInt_Y = fCalo_Y/2.+fPoronVetoLatY_Y/2.;
  fPhysiPoronVetoLatYInt_Z = (fCalo_Z-fPoronVetoLatY_Z)/2.;

  fPhysiVetoLatY_X = 0;
  fPhysiVetoLatY_Y = fCalo_Y/2.+fPoronVetoLatY_Y+fVetoLatY_Y/2.;;
  fPhysiVetoLatY_Z = (fCalo_Z-fVetoLatY_Z)/2.;

  fPhysiPoronVetoLatYExt_X = 0;
  fPhysiPoronVetoLatYExt_Y = fCalo_Y/2.+fPoronVetoLatY_Y+fVetoLatY_Y+fPoronVetoLatY_Y/2.;
  fPhysiPoronVetoLatYExt_Z = (fCalo_Z-fPoronVetoLatY_Z)/2.;

  fPhysiPoronVetoLatXInt_X = fCalo_X/2.+fPoronVetoLatX_X/2.;
  fPhysiPoronVetoLatXInt_Y = 0;
  fPhysiPoronVetoLatXInt_Z = (fCalo_Z-fPoronVetoLatX_Z)/2.;

  fPhysiVetoLatX_X = fCalo_X/2.+fPoronVetoLatX_X+fVetoLatX_X/2.;
  fPhysiVetoLatX_Y = 0;
  fPhysiVetoLatX_Z = (fCalo_Z-fVetoLatX_Z)/2.;

  fPhysiPoronVetoLatXExt_X = fCalo_X/2.+fPoronVetoLatX_X+fVetoLatX_X+fPoronVetoLatX_X/2.;
  fPhysiPoronVetoLatXExt_Y = 0;
  fPhysiPoronVetoLatXExt_Z = (fCalo_Z-fPoronVetoLatX_Z)/2.;

  fPhysiHCSkinLatXInt_X = fCalo_X/2.+fPoronVetoLatX_X+fVetoLatX_X+fPoronVetoLatX_X+fHCSkinLatX_X/2.;
  fPhysiHCSkinLatXInt_Y = 0;
  fPhysiHCSkinLatXInt_Z = (fCalo_Z-fHCSkinLatX_Z)/2.;

  fPhysiHCCoreLatX_X =fCalo_X/2.+fPoronVetoLatX_X+fVetoLatX_X+fPoronVetoLatX_X+fHCSkinLatX_X+fHCCoreLatX_X/2.;
  fPhysiHCCoreLatX_Y =0;
  fPhysiHCCoreLatX_Z =(fCalo_Z-fHCCoreLatX_Z)/2.;;

  fPhysiHCSkinLatXExt_X = fCalo_X/2.+fPoronVetoLatX_X+fVetoLatX_X+fPoronVetoLatX_X+fHCSkinLatX_X+fHCCoreLatX_X+fHCSkinLatX_X/2.;
  fPhysiHCSkinLatXExt_Y = 0;
  fPhysiHCSkinLatXExt_Z = (fCalo_Z-fHCSkinLatX_Z)/2.;

  fPhysiHCSkinLatYInt_X = 0;
  fPhysiHCSkinLatYInt_Y = fCalo_Y/2.+fPoronVetoLatY_Y+fVetoLatY_Y+fPoronVetoLatY_Y+fHCSkinLatY_Y/2.;
  fPhysiHCSkinLatYInt_Z = (fCalo_Z-fHCSkinLatY_Z)/2.;

  fPhysiHCCoreLatY_X = 0;
  fPhysiHCCoreLatY_Y = fCalo_Y/2.+fPoronVetoLatY_Y+fVetoLatY_Y+fPoronVetoLatY_Y+fHCSkinLatY_Y+fHCCoreLatY_Y/2.;
  fPhysiHCCoreLatY_Z = (fCalo_Z-fHCCoreLatY_Z)/2.;;

  fPhysiHCSkinLatYExt_X = 0;
  fPhysiHCSkinLatYExt_Y = fCalo_Y/2.+fPoronVetoLatY_Y+fVetoLatY_Y+fPoronVetoLatY_Y+fHCSkinLatY_Y+fHCCoreLatY_Y+fHCSkinLatY_Y/2.;
  fPhysiHCSkinLatYExt_Z = (fCalo_Z-fHCSkinLatY_Z)/2.;

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterConstructionConfig1::Builder(G4VPhysicalVolume* motherVolume)
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
  G4Material* caloLayerMat = pMaterial->GetMaterial(caloMaterial);
  G4Material* vetoLayerMat = pMaterial->GetMaterial(vetoMaterial);
  G4Material* cfMat        = pMaterial->GetMaterial(cfiberMaterial);
  G4Material* porMat       = pMaterial->GetMaterial(poronMaterial);
  G4Material* hcMat        = pMaterial->GetMaterial(honeycombMaterial);

  fSolidCaloBox = new G4Box("Calorimeter",fCalo_X/2.,fCalo_Y/2.,fCalo_Z/2.);

  fSolidLayer = new G4Box("Layer",fCaloLayer_X/2.,fCaloLayer_Y/2.,fCaloLayer_Z/2.);

  fSolidActiveLayer = new G4Box("ActiveLayer",fActiveLayer_X/2.,fActiveLayer_Y/2.,fActiveLayer_Z/2.);  

  fSolidCFPlane = new  G4Box("CFSupport",fCFPlane_X/2.,fCFPlane_Y/2.,fCFPlane_Z/2.); 

  fSolidCFBS = new G4Box("CFSupport",fCFBS_X/2.,fCFBS_Y/2.,fCFBS_Z/2.);

  fSolidCFSS = new G4Box("CFSupport",fCFSS_X/2.,fCFSS_Y/2.,fCFSS_Z/2.);

  fSolidCFYPStick = new G4Box("CFSupport",fCFYPStick_X/2.,fCFYPStick_Y/2.,fCFYPStcik_Z/2.);

  fSolidCFYSBStick = new G4Box("CFSupport",fCFYSBStick_X/2.,fCFYSBStick_Y/2.,fCFYSBStick_Z/2.);

  fSolidCFYBBStick = new G4Box("CFSupport",fCFYBBStick_X/2.,fCFYBBStick_Y/2.,fCFYBBStick_Z/2.);

  fSolidPoronLayer = new  G4Box("Poron",fPoronL_X/2.,fPoronL_Y/2.,fPoronL_Z/2.); 

  fSolidPoronStick = new  G4Box("Poron",fPoronS_X/2.,fPoronS_Y/2.,fPoronS_Z/2.); 
  
  fSolidVetoBottom = new G4Box("Veto",fVetoBottom_X/2.,fVetoBottom_Y/2.,fVetoBottom_Z/2.);

  fSolidCFStickLat = new G4Box("CFSupport",fCFStickLat_X/2.,fCFStickLat_Y/2.,fCFStickLat_Z/2.);

  fSolidCFStickFront = new G4Box("CFSupport",fCFStickFront_X/2.,fCFStickFront_Y/2.,fCFStickFront_Z/2.);

  fSolidHCBMSkin = new G4Box("HoneyCombSkin",fHCBMSkin_X/2.,fHCBMSkin_Y/2.,fHCBMSkin_Z/2.);

  fSolidHCBLSkin = new G4Box("HoneyCombSkin",fHCBLSkin_X/2.,fHCBLSkin_Y/2.,fHCBLSkin_Z/2.);

  fSolidHCBMCore = new G4Box("HoneyCombCore",fHCBMCore_X/2.,fHCBMCore_Y/2.,fHCBMCore_Z/2.);

  fSolidHCBLCore = new G4Box("HoneyCombCore",fHCBLCore_X/2.,fHCBLCore_Y/2.,fHCBLCore_Z/2.);
  
  fSolidVetoLatY = new G4Box("Veto",fVetoLatY_X/2.,fVetoLatY_Y/2.,fVetoLatY_Z/2.);

  fSolidPoronVetoLatY = new G4Box("Veto",fPoronVetoLatY_X/2.,fPoronVetoLatY_Y/2.,fPoronVetoLatY_Z/2.);

  fSolidVetoLatX = new G4Box("Veto",fVetoLatX_X/2.,fVetoLatX_Y/2.,fVetoLatX_Z/2.);

  fSolidPoronVetoLatX = new G4Box("Veto",fPoronVetoLatX_X/2.,fPoronVetoLatX_Y/2.,fPoronVetoLatX_Z/2.);

  fSolidHCSkinLatX = new G4Box("HoneyCombSkin",fHCSkinLatX_X/2.,fHCSkinLatX_Y/2.,fHCSkinLatX_Z/2.);
		   
  fSolidHCCoreLatX = new G4Box("HoneyCombCore",fHCCoreLatX_X/2.,fHCCoreLatX_Y/2.,fHCCoreLatX_Z/2.);
		   
  fSolidHCSkinLatY = new G4Box("HoneyCombSkin",fHCSkinLatY_X/2.,fHCSkinLatY_Y/2.,fHCSkinLatY_Z/2.);
		   
  fSolidHCCoreLatY = new G4Box("HoneyCombCore",fHCCoreLatY_X/2.,fHCCoreLatY_Y/2.,fHCCoreLatY_Z/2.);


  fLogicCaloBox = new G4LogicalVolume(fSolidCaloBox,vacuum,"Calorimeter");
  
  fLogicLayer = new G4LogicalVolume(fSolidLayer,vacuum,"Layer");

  fLogicActiveLayer = new G4LogicalVolume(fSolidActiveLayer,caloLayerMat,"ActiveLayer");  
  
  fLogicCFPlane = new G4LogicalVolume(fSolidCFPlane,cfMat,"CFSupport");

  fLogicCFBS = new G4LogicalVolume(fSolidCFBS,cfMat,"CFSupport");

  fLogicCFSS = new G4LogicalVolume(fSolidCFSS,cfMat,"CFSupport");

  fLogicCFYPStick = new G4LogicalVolume(fSolidCFYPStick,cfMat,"CFSupport");

  fLogicCFYSBStick = new G4LogicalVolume(fSolidCFYSBStick,cfMat,"CFSupport");

  fLogicCFYBBStick = new G4LogicalVolume(fSolidCFYBBStick,cfMat,"CFSupport");

  fLogicPoronLayer = new G4LogicalVolume(fSolidPoronLayer,porMat,"Poron");
  
  fLogicPoronStick = new G4LogicalVolume(fSolidPoronStick,porMat,"Poron");
  
  fLogicVetoBottom = new G4LogicalVolume(fSolidVetoBottom,vetoLayerMat,"Veto");

  fLogicCFStickLat = new G4LogicalVolume(fSolidCFStickLat,cfMat,"CFSupport");
  
  fLogicCFStickFront = new G4LogicalVolume(fSolidCFStickFront,cfMat,"CFSupport");

  fLogicHCBMSkin = new G4LogicalVolume(fSolidHCBMSkin,cfMat,"HoneyCombSkin");
  
  fLogicHCBLSkin = new G4LogicalVolume(fSolidHCBLSkin,cfMat,"HoneyCombSkin");

  fLogicHCBMCore = new G4LogicalVolume(fSolidHCBMCore,hcMat,"HoneyCombCore");
  
  fLogicHCBLCore = new G4LogicalVolume(fSolidHCBLCore,hcMat,"HoneyCombCore");
  
  fLogicVetoLatY = new G4LogicalVolume(fSolidVetoLatY,vetoLayerMat,"Veto");

  fLogicPoronVetoLatY = new G4LogicalVolume(fSolidPoronVetoLatY,porMat,"Poron");

  fLogicVetoLatX = new G4LogicalVolume(fSolidVetoLatX,vetoLayerMat,"Veto");

  fLogicPoronVetoLatX = new G4LogicalVolume(fSolidPoronVetoLatX,porMat,"Poron");

  fLogicHCSkinLatX = new G4LogicalVolume(fSolidHCSkinLatX,cfMat,"HoneyCombSkin");
  					 		   
  fLogicHCCoreLatX = new G4LogicalVolume(fSolidHCCoreLatX,hcMat,"HoneyCombCore");
  					 		   
  fLogicHCSkinLatY = new G4LogicalVolume(fSolidHCSkinLatY,cfMat,"HoneyCombSkin");
  					 		   
  fLogicHCCoreLatY = new G4LogicalVolume(fSolidHCCoreLatY,hcMat,"HoneyCombCore");

  fLogicActiveLayer->SetSensitiveDetector(caloSD);
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
  
  fPhysiLayer = new G4PVReplica("Layer", 
				fLogicLayer,            
				fPhysiCaloBox,            
				kZAxis,                 
				fNbOfCaloLayers,        
				fCaloLayer_Z);   
  
  fPhysiCFPlaneU = new G4PVPlacement(0,                   
				     G4ThreeVector(fPhysiCFPlaneU_X,fPhysiCFPlaneU_Y,fPhysiCFPlaneU_Z),     
				     "CFSupport",              
				     fLogicCFPlane,       
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

  fPhysiCFPlaneD = new G4PVPlacement(0,                   
				     G4ThreeVector(fPhysiCFPlaneD_X,fPhysiCFPlaneD_Y,fPhysiCFPlaneD_Z),     
				     "CFSupport",              
				     fLogicCFPlane,       
				     fPhysiLayer,        
				     false,               
				     0,true); 

  fPhysiPoronLayer = new G4PVPlacement(0,                   
				       G4ThreeVector(fPhysiPoronLayer_X,fPhysiPoronLayer_Y,fPhysiPoronLayer_Z),     
				       "Poron",              
				       fLogicPoronLayer,       
				       fPhysiLayer,        
				       false,               
				       0,true); 

  fPhysiPoronStickR = new G4PVPlacement(0,                   
					G4ThreeVector(fPhysiPoronStick_X,fPhysiPoronStick_Y,fPhysiPoronStick_Z),
					"Poron",              
					fLogicPoronStick,       
					fPhysiLayer,        
					false,               
					0,true); 

  fPhysiPoronStickL = new G4PVPlacement(0,                   
					G4ThreeVector(fPhysiPoronStick_X,-fPhysiPoronStick_Y,fPhysiPoronStick_Z),
					"Poron",              
					fLogicPoronStick,       
					fPhysiLayer,        
					false,               
					0,true); 

  fPhysiCFBSR = new G4PVPlacement(0,                   
				  G4ThreeVector(fPhysiCFBS_X,fPhysiCFBS_Y,fPhysiCFBS_Z),
				  "CFSupport",              
				  fLogicCFBS,       
				  fPhysiLayer,        
				  false,               
				  0,true); 

  fPhysiCFBSL = new G4PVPlacement(0,                   
				  G4ThreeVector(fPhysiCFBS_X,-fPhysiCFBS_Y,fPhysiCFBS_Z),
				  "CFSupport",              
				  fLogicCFBS,       
				  fPhysiLayer,        
				  false,               
				  0,true); 

  fPhysiCFSSR = new G4PVPlacement(0,                   
				  G4ThreeVector(fPhysiCFSS_X,fPhysiCFSS_Y,fPhysiCFSS_Z),
				  "CFSupport",              
				  fLogicCFSS,       
				  fPhysiLayer,        
				  false,               
				  0,true); 

  fPhysiCFSSL = new G4PVPlacement(0,                   
				  G4ThreeVector(fPhysiCFSS_X,-fPhysiCFSS_Y,fPhysiCFSS_Z),
				  "CFSupport",              
				  fLogicCFSS,       
				  fPhysiLayer,        
				  false,               
				  0,true); 

  fPhysiCFYPStickUP = new G4PVPlacement(0,                   
					G4ThreeVector(fPhysiCFYPStickU_X,fPhysiCFYPStickU_Y,fPhysiCFYPStickU_Z),
					"CFSupport",              
					fLogicCFYPStick,       
					fPhysiLayer,        
					false,               
					0,true); 

  fPhysiCFYPStickUM = new G4PVPlacement(0,                   
					G4ThreeVector(-fPhysiCFYPStickU_X,fPhysiCFYPStickU_Y,fPhysiCFYPStickU_Z),     
					"CFSupport",              
					fLogicCFYPStick,       
					fPhysiLayer,        
					false,               
					0,true); 

  fPhysiCFYBBStickP = new G4PVPlacement(0,                   
					G4ThreeVector(fPhysiCFYBBStick_X,fPhysiCFYBBStick_Y,fPhysiCFYBBStick_Z),
					"CFSupport",              
					fLogicCFYBBStick,       
					fPhysiLayer,        
					false,               
					0,true); 

  fPhysiCFYBBStickM = new G4PVPlacement(0,                   
					G4ThreeVector(-fPhysiCFYBBStick_X,fPhysiCFYBBStick_Y,fPhysiCFYBBStick_Z),
					"CFSupport",              
					fLogicCFYBBStick,       
					fPhysiLayer,        
					false,               
					0,true); 

  fPhysiCFYPStickDP = new G4PVPlacement(0,                   
					G4ThreeVector(fPhysiCFYPStickD_X,fPhysiCFYPStickD_Y,fPhysiCFYPStickD_Z),
					"CFSupport",              
					fLogicCFYPStick,       
					fPhysiLayer,        
					false,               
					0,true); 

  fPhysiCFYPStickDM = new G4PVPlacement(0,                   
					G4ThreeVector(-fPhysiCFYPStickD_X,fPhysiCFYPStickD_Y,fPhysiCFYPStickD_Z),
					"CFSupport",              
					fLogicCFYPStick,       
					fPhysiLayer,        
					false,               
					0,true); 

  fPhysiCFYSBStickP = new G4PVPlacement(0,                   
					G4ThreeVector(fPhysiCFYSBStick_X,fPhysiCFYSBStick_Y,fPhysiCFYSBStick_Z),
					"CFSupport",              
					fLogicCFYSBStick,       
					fPhysiLayer,        
					false,               
					0,true); 

  fPhysiCFYSBStickM = new G4PVPlacement(0,                   
					G4ThreeVector(-fPhysiCFYSBStick_X,fPhysiCFYSBStick_Y,fPhysiCFYSBStick_Z),     
					"CFSupport",              
					fLogicCFYSBStick,       
					fPhysiLayer,        
					false,               
					0,true); 

  fPhysiVetoBottom = new G4PVPlacement(0,
				       G4ThreeVector(fPhysiVetoBottom_X,fPhysiVetoBottom_Y,fPhysiVetoBottom_Z),
				       "VetoB",
				       fLogicVetoBottom,
				       motherVolume,
				       false,
				       0,true);

  fPhysiPoronVetoBottom = new G4PVPlacement(0,
					    G4ThreeVector(fPhysiPoronVetoBottom_X,fPhysiPoronVetoBottom_Y,fPhysiPoronVetoBottom_Z),
					    "Poron",
					    fLogicPoronLayer,
					    motherVolume,
					    false,
					    0,true);

  fPhysiPoronVetoLatR = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiPoronVetoLat_X,fPhysiPoronVetoLat_Y,fPhysiPoronVetoLat_Z),
					  "Poron",
					  fLogicPoronStick,
					  motherVolume,
					  false,
					  0,true);

  fPhysiPoronVetoLatL = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiPoronVetoLat_X,-fPhysiPoronVetoLat_Y,fPhysiPoronVetoLat_Z),
					  "Poron",
					  fLogicPoronStick,
					  motherVolume,
					  false,
					  0,true);

  fPhysiCFStickLatR = new G4PVPlacement(0,
					G4ThreeVector(fPhysiCFStickLat_X,fPhysiCFStickLat_Y,fPhysiCFStickLat_Z),
					"CFSupport",
					fLogicCFStickLat,
					motherVolume,
					false,
					0,true);

  fPhysiCFStickLatL = new G4PVPlacement(0,
					G4ThreeVector(fPhysiCFStickLat_X,-fPhysiCFStickLat_Y,fPhysiCFStickLat_Z),
					"CFSupport",
					fLogicCFStickLat,
					motherVolume,
					false,
					0,true);
  
  fPhysiCFStickFrontP = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiCFStickFront_X,fPhysiCFStickFront_Y,fPhysiCFStickFront_Z),
					  "CFSupport",
					  fLogicCFStickFront,
					  motherVolume,
					  false,
					  0,true);

  fPhysiCFStickFrontM = new G4PVPlacement(0,
					  G4ThreeVector(-fPhysiCFStickFront_X,fPhysiCFStickFront_Y,fPhysiCFStickFront_Z),
					  "CFSupport",
					  fLogicCFStickFront,
					  motherVolume,
					  false,
					  0,true);
  
  fPhysiHCMSkinU = new G4PVPlacement(0,
				     G4ThreeVector(fPhysiHCMSkinU_X,fPhysiHCMSkinU_Y,fPhysiHCMSkinU_Z),
				     "HoneyCombSkin",
				     fLogicHCBMSkin,
				     motherVolume,
				     false,
				     0,true);

  fPhysiHCMCore = new G4PVPlacement(0,
				    G4ThreeVector(fPhysiHCMCore_X,fPhysiHCMCore_Y,fPhysiHCMCore_Z),
				    "HoneyCombCore",
				    fLogicHCBMCore,
				    motherVolume,
				    false,
				    0,true);

  fPhysiHCMSkinD = new G4PVPlacement(0,
				     G4ThreeVector(fPhysiHCMSkinD_X,fPhysiHCMSkinD_Y,fPhysiHCMSkinD_Z),
				     "HoneyCombSkin",
				     fLogicHCBMSkin,
				     motherVolume,
				     false,
				     0,true);

  fPhysiHCLSkinUR = new G4PVPlacement(0,
				      G4ThreeVector(fPhysiHCLSkinU_X,fPhysiHCLSkinU_Y,fPhysiHCLSkinU_Z),
				      "HoneyCombSkin",
				      fLogicHCBLSkin,
				      motherVolume,
				      false,
				      0,true);

  fPhysiHCLCoreR = new G4PVPlacement(0,
				     G4ThreeVector(fPhysiHCLCore_X,fPhysiHCLCore_Y,fPhysiHCLCore_Z),
				     "HoneyCombCore",
				     fLogicHCBLCore,
				     motherVolume,
				     false,
				     0,true);

  fPhysiHCLSkinDR = new G4PVPlacement(0,
				      G4ThreeVector(fPhysiHCLSkinD_X,fPhysiHCLSkinD_Y,fPhysiHCLSkinD_Z),
				      "HoneyCombSkin",
				      fLogicHCBLSkin,
				      motherVolume,
				      false,
				      0,true);

  fPhysiHCLSkinUL = new G4PVPlacement(0,
				      G4ThreeVector(-fPhysiHCLSkinU_X,fPhysiHCLSkinU_Y,fPhysiHCLSkinU_Z),
				      "HoneyCombSkin",
				      fLogicHCBLSkin,
				      motherVolume,
				      false,
				      0,true);

  fPhysiHCLCoreL = new G4PVPlacement(0,
				     G4ThreeVector(-fPhysiHCLCore_X,fPhysiHCLCore_Y,fPhysiHCLCore_Z),
				     "HoneyCombCore",
				     fLogicHCBLCore,
				     motherVolume,
				     false,
				     0,true);

  fPhysiHCLSkinDL = new G4PVPlacement(0,
				      G4ThreeVector(-fPhysiHCLSkinD_X,fPhysiHCLSkinD_Y,fPhysiHCLSkinD_Z),
				      "HoneyCombSkin",
				      fLogicHCBLSkin,
				      motherVolume,
				      false,
				      0,true);

  fPhysiPoronVetoLatYPInt = new G4PVPlacement(0,
					      G4ThreeVector(fPhysiPoronVetoLatYInt_X,fPhysiPoronVetoLatYInt_Y,fPhysiPoronVetoLatYInt_Z),
					      "Poron",
					      fLogicPoronVetoLatY,
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
  
  fPhysiPoronVetoLatYPExt = new G4PVPlacement(0,
					      G4ThreeVector(fPhysiPoronVetoLatYExt_X,fPhysiPoronVetoLatYExt_Y,fPhysiPoronVetoLatYExt_Z),
					      "Poron",
					      fLogicPoronVetoLatY,
					      motherVolume,
					      false,
					      0,true);

  fPhysiPoronVetoLatYMInt = new G4PVPlacement(0,
					      G4ThreeVector(fPhysiPoronVetoLatYInt_X,-fPhysiPoronVetoLatYInt_Y,fPhysiPoronVetoLatYInt_Z),
					      "Poron",
					      fLogicPoronVetoLatY,
					      motherVolume,
					      false,
					      0,true);
  fPhysiVetoLatYM = new G4PVPlacement(0,
				      G4ThreeVector(fPhysiVetoLatY_X,-fPhysiVetoLatY_Y,fPhysiVetoLatY_Z),
				      "VetoLYM",
				      fLogicVetoLatY,
				      motherVolume,
				      false,
				      0,true);
  
  fPhysiPoronVetoLatYMExt = new G4PVPlacement(0,
					      G4ThreeVector(fPhysiPoronVetoLatYExt_X,-fPhysiPoronVetoLatYExt_Y,fPhysiPoronVetoLatYExt_Z),
					      "Poron",
					      fLogicPoronVetoLatY,
					      motherVolume,
					      false,
					      0,true);

  fPhysiPoronVetoLatXPInt = new G4PVPlacement(0,
					      G4ThreeVector(fPhysiPoronVetoLatXInt_X,fPhysiPoronVetoLatXInt_Y,fPhysiPoronVetoLatXInt_Z),
					      "Poron",
					      fLogicPoronVetoLatX,
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
  
  fPhysiPoronVetoLatXPExt = new G4PVPlacement(0,
					      G4ThreeVector(fPhysiPoronVetoLatXExt_X,fPhysiPoronVetoLatXExt_Y,fPhysiPoronVetoLatXExt_Z),
					      "Poron",
					      fLogicPoronVetoLatX,
					      motherVolume,
					      false,
					      0,true);

  fPhysiPoronVetoLatXMInt = new G4PVPlacement(0,
					      G4ThreeVector(-fPhysiPoronVetoLatXInt_X,fPhysiPoronVetoLatXInt_Y,fPhysiPoronVetoLatXInt_Z),
					      "Poron",
					      fLogicPoronVetoLatX,
					      motherVolume,
					      false,
					      0,true);
  fPhysiVetoLatXM = new G4PVPlacement(0,
				      G4ThreeVector(-fPhysiVetoLatX_X,fPhysiVetoLatX_Y,fPhysiVetoLatX_Z),
				      "VetoLXM",
				      fLogicVetoLatX,
				      motherVolume,
				      false,
				      0,true);
  
  fPhysiPoronVetoLatXMExt = new G4PVPlacement(0,
					      G4ThreeVector(-fPhysiPoronVetoLatXExt_X,fPhysiPoronVetoLatXExt_Y,fPhysiPoronVetoLatXExt_Z),
					      "Poron",
					      fLogicPoronVetoLatX,
					      motherVolume,
					      false,
					      0,true);
					      
  fPhysiHCSkinLatXPInt = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiHCSkinLatXInt_X,fPhysiHCSkinLatXInt_Y,fPhysiHCSkinLatXInt_Z),
					   "HoneyCombSkin",
					   fLogicHCSkinLatX,
					   motherVolume,
					   false,
					   0,true);
  fPhysiHCCoreLatXP = new G4PVPlacement(0,
					G4ThreeVector(fPhysiHCCoreLatX_X,fPhysiHCCoreLatX_Y,fPhysiHCCoreLatX_Z),
					"HoneyCombCore",
					fLogicHCCoreLatX,
					motherVolume,
					false,
					0,true);
  				
  fPhysiHCSkinLatXPExt = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiHCSkinLatXExt_X,fPhysiHCSkinLatXExt_Y,fPhysiHCSkinLatXExt_Z),
					   "HoneyCombSkin",
					   fLogicHCSkinLatX,
					   motherVolume,
					   false,
					   0,true);
				
  fPhysiHCSkinLatXMInt = new G4PVPlacement(0,
					   G4ThreeVector(-fPhysiHCSkinLatXInt_X,fPhysiHCSkinLatXInt_Y,fPhysiHCSkinLatXInt_Z),
					   "HoneyCombSkin",
					   fLogicHCSkinLatX,
					   motherVolume,
					   false,
					   0,true);
  fPhysiHCCoreLatXM = new G4PVPlacement(0,
					G4ThreeVector(-fPhysiHCCoreLatX_X,fPhysiHCCoreLatX_Y,fPhysiHCCoreLatX_Z),
					"HoneyCombCore",
					fLogicHCCoreLatX,
					motherVolume,
					false,
					0,true);
  				
  fPhysiHCSkinLatXMExt = new G4PVPlacement(0,
					   G4ThreeVector(-fPhysiHCSkinLatXExt_X,fPhysiHCSkinLatXExt_Y,fPhysiHCSkinLatXExt_Z),
					   "HoneyCombSkin",
					   fLogicHCSkinLatX,
					   motherVolume,
					   false,
					   0,true);
				
				
  fPhysiHCSkinLatYPInt = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiHCSkinLatYInt_X,fPhysiHCSkinLatYInt_Y,fPhysiHCSkinLatYInt_Z),
					   "HoneyCombSkin",
					   fLogicHCSkinLatY,
					   motherVolume,
					   false,
					   0,true);
  fPhysiHCCoreLatYP = new G4PVPlacement(0,
					G4ThreeVector(fPhysiHCCoreLatY_X,fPhysiHCCoreLatY_Y,fPhysiHCCoreLatY_Z),
					"HoneyCombCore",
					fLogicHCCoreLatY,
					motherVolume,
					false,
					0,true);
  				
  fPhysiHCSkinLatYPExt = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiHCSkinLatYExt_X,fPhysiHCSkinLatYExt_Y,fPhysiHCSkinLatYExt_Z),
					   "HoneyCombSkin",
					   fLogicHCSkinLatY,
					   motherVolume,
					   false,
					   0,true);
				
  fPhysiHCSkinLatYMInt = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiHCSkinLatYInt_X,-fPhysiHCSkinLatYInt_Y,fPhysiHCSkinLatYInt_Z),
					   "HoneyCombSkin",
					   fLogicHCSkinLatY,
					   motherVolume,
					   false,
					   0,true);
  fPhysiHCCoreLatYM = new G4PVPlacement(0,
					G4ThreeVector(fPhysiHCCoreLatY_X,-fPhysiHCCoreLatY_Y,fPhysiHCCoreLatY_Z),
					"HoneyCombCore",
					fLogicHCCoreLatY,
					motherVolume,
					false,
					0,true);
  				
  fPhysiHCSkinLatYMExt = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiHCSkinLatYExt_X,-fPhysiHCSkinLatYExt_Y,fPhysiHCSkinLatYExt_Z),
					   "HoneyCombSkin",
					   fLogicHCSkinLatY,
					   motherVolume,
					   false,
					   0,true);
  
  //Visualization Attribute
  G4VisAttributes* attCyan = new G4VisAttributes(G4Colour::Cyan());
  attCyan->SetVisibility(true);
  attCyan->SetForceAuxEdgeVisible(true);
  fLogicHCBMCore->SetVisAttributes(attCyan);
  fLogicHCBLCore->SetVisAttributes(attCyan);
  fLogicHCCoreLatX->SetVisAttributes(attCyan);
  fLogicHCCoreLatY->SetVisAttributes(attCyan);

  G4VisAttributes* attMagenta = new G4VisAttributes(G4Colour::Magenta());
  attMagenta->SetVisibility(true);
  attMagenta->SetForceAuxEdgeVisible(true);
  fLogicActiveLayer->SetVisAttributes(attMagenta);

  G4VisAttributes* attYellow = new G4VisAttributes(G4Colour::Yellow());
  attYellow->SetVisibility(true);
  attYellow->SetForceAuxEdgeVisible(true);
  fLogicPoronLayer->SetVisAttributes(attYellow);
  fLogicPoronStick->SetVisAttributes(attYellow);
  fLogicPoronVetoLatY->SetVisAttributes(attYellow);
  fLogicPoronVetoLatX->SetVisAttributes(attYellow);

  G4VisAttributes* attBrown = new G4VisAttributes(G4Colour::Brown());
  attBrown->SetVisibility(true);
  attBrown->SetForceAuxEdgeVisible(true);
  fLogicCFPlane->SetVisAttributes(attBrown);
  fLogicCFBS->SetVisAttributes(attBrown);
  fLogicCFSS->SetVisAttributes(attBrown);
  fLogicCFYSBStick->SetVisAttributes(attBrown);
  fLogicCFYPStick->SetVisAttributes(attBrown);
  fLogicCFYBBStick->SetVisAttributes(attBrown);
  fLogicCFStickLat->SetVisAttributes(attBrown);
  fLogicCFStickFront->SetVisAttributes(attBrown);
  fLogicHCBMSkin->SetVisAttributes(attBrown);
  fLogicHCBLSkin->SetVisAttributes(attBrown);
  fLogicHCSkinLatX->SetVisAttributes(attBrown);
  fLogicHCSkinLatY->SetVisAttributes(attBrown);

  G4VisAttributes* attGreen = new G4VisAttributes(G4Colour::Green());
  attGreen->SetVisibility(true);
  attGreen->SetForceAuxEdgeVisible(true);
  fLogicVetoBottom->SetVisAttributes(attGreen);
  fLogicVetoLatY->SetVisAttributes(attGreen);
  fLogicVetoLatX->SetVisAttributes(attGreen);

  G4VisAttributes * attInvisible = new G4VisAttributes();
  attInvisible->SetVisibility(false);
  attInvisible->SetForceAuxEdgeVisible(false);
  fLogicCaloBox->SetVisAttributes(attInvisible); 
  fLogicLayer->SetVisAttributes(attInvisible); 
}


