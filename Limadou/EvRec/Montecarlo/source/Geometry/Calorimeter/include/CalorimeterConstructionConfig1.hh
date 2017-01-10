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
/// \file electromagnetic/TestEm3/include/CalorimeterConstructionConfig1.hh
/// \brief Definition of the CalorimeterConstructionConfig1 class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CalorimeterConstructionConfig1_h
#define CalorimeterConstructionConfig1_h 1

#include "globals.hh"

class G4IntersectionSolid;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class HEPDSWMaterial;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CalorimeterConstructionConfig1 
{
public:
  
   CalorimeterConstructionConfig1();
  ~CalorimeterConstructionConfig1();

  inline void SetVetoMaterial(G4String aMat){vetoMaterial=aMat;}
  inline void SetCaloMaterial(G4String aMat){caloMaterial=aMat;}
  inline void SetPoronMaterial(G4String aMat){poronMaterial=aMat;}
  inline void SetCarbonFiberMaterial(G4String aMat){cfiberMaterial=aMat;}
  inline void SetHoneyCombMaterial(G4String aMat){honeycombMaterial=aMat;}

  void Builder(G4VPhysicalVolume* motherVolume);
  
private:
  void ComputeObjectsPositioning();

  G4int              fNbOfCaloLayers;

  G4String           vetoMaterial;
  G4String           caloMaterial;
  G4String           poronMaterial;
  G4String           cfiberMaterial;
  G4String           honeycombMaterial;

  G4double  fCaloLayer_X;
  G4double  fCaloLayer_Y;
  G4double  fCaloLayer_Z;

  G4double  fCalo_X;
  G4double  fCalo_Y;
  G4double  fCalo_Z;

  G4double  fCFPlane_X;
  G4double  fCFPlane_Y;
  G4double  fCFPlane_Z;

  G4double  fCFBS_X;
  G4double  fCFBS_Y;
  G4double  fCFBS_Z;

  G4double  fCFSS_X;
  G4double  fCFSS_Y;
  G4double  fCFSS_Z;

  G4double  fCFYPStick_X;
  G4double  fCFYPStick_Y;
  G4double  fCFYPStcik_Z;

  G4double  fCFYSBStick_X;
  G4double  fCFYSBStick_Y;
  G4double  fCFYSBStick_Z;

  G4double  fCFYBBStick_X;
  G4double  fCFYBBStick_Y;
  G4double  fCFYBBStick_Z;

  G4double  fPoronL_X;
  G4double  fPoronL_Y;
  G4double  fPoronL_Z;
  
  G4double  fPoronS_X;
  G4double  fPoronS_Y;
  G4double  fPoronS_Z;

  G4double  fActiveLayer_X;
  G4double  fActiveLayer_Y;
  G4double  fActiveLayer_Z;

  G4double  fVetoBottom_X;
  G4double  fVetoBottom_Y;
  G4double  fVetoBottom_Z;

  G4double  fCFStickLat_X;
  G4double  fCFStickLat_Y;
  G4double  fCFStickLat_Z;

  G4double  fCFStickFront_X;
  G4double  fCFStickFront_Y;
  G4double  fCFStickFront_Z;

  G4double  fHCBMSkin_X;
  G4double  fHCBMSkin_Y;
  G4double  fHCBMSkin_Z;

  G4double  fHCBLSkin_X;
  G4double  fHCBLSkin_Y;
  G4double  fHCBLSkin_Z;

  G4double  fHCBMCore_X;
  G4double  fHCBMCore_Y;
  G4double  fHCBMCore_Z;

  G4double  fHCBLCore_X;
  G4double  fHCBLCore_Y;
  G4double  fHCBLCore_Z;

  G4double fVetoLatY_X;
  G4double fVetoLatY_Y;
  G4double fVetoLatY_Z;
  
  G4double fPoronVetoLatY_X;
  G4double fPoronVetoLatY_Y;
  G4double fPoronVetoLatY_Z;
  
  G4double fVetoLatX_X;
  G4double fVetoLatX_Y;
  G4double fVetoLatX_Z;
  
  G4double fPoronVetoLatX_X;
  G4double fPoronVetoLatX_Y;
  G4double fPoronVetoLatX_Z;

  G4double fHCSkinLatX_X;
  G4double fHCSkinLatX_Y;
  G4double fHCSkinLatX_Z;

  G4double fHCCoreLatX_X;
  G4double fHCCoreLatX_Y;
  G4double fHCCoreLatX_Z;

  G4double fHCSkinLatY_X;
  G4double fHCSkinLatY_Y;
  G4double fHCSkinLatY_Z;

  G4double fHCCoreLatY_X;
  G4double fHCCoreLatY_Y;
  G4double fHCCoreLatY_Z;
  

  G4double fPhysiCFPlaneU_X;
  G4double fPhysiCFPlaneU_Y;
  G4double fPhysiCFPlaneU_Z;

  G4double fPhysiActiveLayer_X;
  G4double fPhysiActiveLayer_Y;
  G4double fPhysiActiveLayer_Z;

  G4double fPhysiCFPlaneD_X;
  G4double fPhysiCFPlaneD_Y;
  G4double fPhysiCFPlaneD_Z;

  G4double fPhysiPoronLayer_X;
  G4double fPhysiPoronLayer_Y;
  G4double fPhysiPoronLayer_Z;

  G4double fPhysiPoronStick_X;
  G4double fPhysiPoronStick_Y;
  G4double fPhysiPoronStick_Z;

  G4double fPhysiCFBS_X;
  G4double fPhysiCFBS_Y;
  G4double fPhysiCFBS_Z;

  G4double fPhysiCFSS_X;
  G4double fPhysiCFSS_Y;
  G4double fPhysiCFSS_Z;

  G4double fPhysiCFYPStickU_X;
  G4double fPhysiCFYPStickU_Y;
  G4double fPhysiCFYPStickU_Z;

  G4double fPhysiCFYBBStick_X;
  G4double fPhysiCFYBBStick_Y;
  G4double fPhysiCFYBBStick_Z;

  G4double fPhysiCFYPStickD_X;
  G4double fPhysiCFYPStickD_Y;
  G4double fPhysiCFYPStickD_Z;

  G4double fPhysiCFYSBStick_X;
  G4double fPhysiCFYSBStick_Y;
  G4double fPhysiCFYSBStick_Z;

  G4double fPhysiVetoBottom_X;
  G4double fPhysiVetoBottom_Y;
  G4double fPhysiVetoBottom_Z;

  G4double fPhysiPoronVetoBottom_X;
  G4double fPhysiPoronVetoBottom_Y;
  G4double fPhysiPoronVetoBottom_Z;

  G4double fPhysiPoronVetoLat_X;
  G4double fPhysiPoronVetoLat_Y;
  G4double fPhysiPoronVetoLat_Z;

  G4double fPhysiCFStickLat_X;
  G4double fPhysiCFStickLat_Y;
  G4double fPhysiCFStickLat_Z;

  G4double fPhysiCFStickFront_X;
  G4double fPhysiCFStickFront_Y;
  G4double fPhysiCFStickFront_Z;

  G4double fPhysiHCMSkinU_X;
  G4double fPhysiHCMSkinU_Y;
  G4double fPhysiHCMSkinU_Z;

  G4double fPhysiHCMCore_X;
  G4double fPhysiHCMCore_Y;
  G4double fPhysiHCMCore_Z;

  G4double fPhysiHCMSkinD_X;
  G4double fPhysiHCMSkinD_Y;
  G4double fPhysiHCMSkinD_Z;

  G4double fPhysiHCLSkinU_X;
  G4double fPhysiHCLSkinU_Y;
  G4double fPhysiHCLSkinU_Z;

  G4double fPhysiHCLCore_X;
  G4double fPhysiHCLCore_Y;
  G4double fPhysiHCLCore_Z;

  G4double fPhysiHCLSkinD_X;
  G4double fPhysiHCLSkinD_Y;
  G4double fPhysiHCLSkinD_Z;

  G4double fPhysiPoronVetoLatYInt_X;
  G4double fPhysiPoronVetoLatYInt_Y;
  G4double fPhysiPoronVetoLatYInt_Z;

  G4double fPhysiVetoLatY_X;
  G4double fPhysiVetoLatY_Y;
  G4double fPhysiVetoLatY_Z;

  G4double fPhysiPoronVetoLatYExt_X;
  G4double fPhysiPoronVetoLatYExt_Y;
  G4double fPhysiPoronVetoLatYExt_Z;

  G4double fPhysiPoronVetoLatXInt_X;
  G4double fPhysiPoronVetoLatXInt_Y;
  G4double fPhysiPoronVetoLatXInt_Z;

  G4double fPhysiVetoLatX_X;
  G4double fPhysiVetoLatX_Y;
  G4double fPhysiVetoLatX_Z;

  G4double fPhysiPoronVetoLatXExt_X;
  G4double fPhysiPoronVetoLatXExt_Y;
  G4double fPhysiPoronVetoLatXExt_Z;

  G4double fPhysiHCSkinLatXInt_X;
  G4double fPhysiHCSkinLatXInt_Y;
  G4double fPhysiHCSkinLatXInt_Z;

  G4double fPhysiHCCoreLatX_X;
  G4double fPhysiHCCoreLatX_Y;
  G4double fPhysiHCCoreLatX_Z;

  G4double fPhysiHCSkinLatXExt_X;
  G4double fPhysiHCSkinLatXExt_Y;
  G4double fPhysiHCSkinLatXExt_Z;

  G4double fPhysiHCSkinLatYInt_X;
  G4double fPhysiHCSkinLatYInt_Y;
  G4double fPhysiHCSkinLatYInt_Z;

  G4double fPhysiHCCoreLatY_X;
  G4double fPhysiHCCoreLatY_Y;
  G4double fPhysiHCCoreLatY_Z;

  G4double fPhysiHCSkinLatYExt_X;
  G4double fPhysiHCSkinLatYExt_Y;
  G4double fPhysiHCSkinLatYExt_Z;

  HEPDSWMaterial*      pMaterial;

  G4Box* fSolidCaloBox;
  G4Box* fSolidLayer;
  G4Box* fSolidActiveLayer; 
  G4Box* fSolidCFPlane;
  G4Box* fSolidCFBS;
  G4Box* fSolidCFSS;
  G4Box* fSolidCFYPStick;
  G4Box* fSolidCFYSBStick;
  G4Box* fSolidCFYBBStick;
  G4Box* fSolidPoronLayer;
  G4Box* fSolidPoronStick;
  G4Box* fSolidVetoBottom;
  G4Box* fSolidCFStickLat;
  G4Box* fSolidCFStickFront;
  G4Box* fSolidHCBMSkin;
  G4Box* fSolidHCBLSkin;
  G4Box* fSolidHCBMCore;
  G4Box* fSolidHCBLCore;
  G4Box* fSolidVetoLatY;
  G4Box* fSolidPoronVetoLatY;
  G4Box* fSolidVetoLatX;
  G4Box* fSolidPoronVetoLatX;
  G4Box* fSolidHCSkinLatX;
  G4Box* fSolidHCCoreLatX;
  G4Box* fSolidHCSkinLatY;
  G4Box* fSolidHCCoreLatY;

  G4LogicalVolume* fLogicCaloBox;
  G4LogicalVolume* fLogicLayer;
  G4LogicalVolume* fLogicActiveLayer;
  G4LogicalVolume* fLogicCFPlane;
  G4LogicalVolume* fLogicCFBS;
  G4LogicalVolume* fLogicCFSS;
  G4LogicalVolume* fLogicCFYPStick;
  G4LogicalVolume* fLogicCFYSBStick;
  G4LogicalVolume* fLogicCFYBBStick;
  G4LogicalVolume* fLogicPoronLayer;
  G4LogicalVolume* fLogicPoronStick;
  G4LogicalVolume* fLogicVetoBottom;
  G4LogicalVolume* fLogicCFStickLat;
  G4LogicalVolume* fLogicCFStickFront;
  G4LogicalVolume* fLogicHCBMSkin;
  G4LogicalVolume* fLogicHCBLSkin;
  G4LogicalVolume* fLogicHCBMCore;
  G4LogicalVolume* fLogicHCBLCore;
  G4LogicalVolume* fLogicVetoLatY;
  G4LogicalVolume* fLogicPoronVetoLatY;
  G4LogicalVolume* fLogicVetoLatX;
  G4LogicalVolume* fLogicPoronVetoLatX;
  G4LogicalVolume* fLogicHCSkinLatX;
  G4LogicalVolume* fLogicHCCoreLatX;
  G4LogicalVolume* fLogicHCSkinLatY;
  G4LogicalVolume* fLogicHCCoreLatY;

  G4VPhysicalVolume* fPhysiCaloBox;
  G4VPhysicalVolume* fPhysiLayer;
  G4VPhysicalVolume* fPhysiActiveLayer;
  G4VPhysicalVolume* fPhysiCFPlaneU;
  G4VPhysicalVolume* fPhysiCFPlaneD;
  G4VPhysicalVolume* fPhysiPoronLayer;
  G4VPhysicalVolume* fPhysiPoronStickR;
  G4VPhysicalVolume* fPhysiPoronStickL;
  G4VPhysicalVolume* fPhysiCFBSR;
  G4VPhysicalVolume* fPhysiCFBSL;
  G4VPhysicalVolume* fPhysiCFSSR;
  G4VPhysicalVolume* fPhysiCFSSL;
  G4VPhysicalVolume* fPhysiCFYPStickUP;
  G4VPhysicalVolume* fPhysiCFYPStickUM;
  G4VPhysicalVolume* fPhysiCFYBBStickP;
  G4VPhysicalVolume* fPhysiCFYBBStickM;
  G4VPhysicalVolume* fPhysiCFYPStickDP;
  G4VPhysicalVolume* fPhysiCFYPStickDM;
  G4VPhysicalVolume* fPhysiCFYSBStickP;
  G4VPhysicalVolume* fPhysiCFYSBStickM;
  G4VPhysicalVolume* fPhysiVetoBottom;
  G4VPhysicalVolume* fPhysiPoronVetoBottom;
  G4VPhysicalVolume* fPhysiPoronVetoLatR;
  G4VPhysicalVolume* fPhysiPoronVetoLatL;
  G4VPhysicalVolume* fPhysiCFStickLatR;
  G4VPhysicalVolume* fPhysiCFStickLatL;
  G4VPhysicalVolume* fPhysiCFStickFrontP;
  G4VPhysicalVolume* fPhysiCFStickFrontM;
  G4VPhysicalVolume* fPhysiHCMSkinU;
  G4VPhysicalVolume* fPhysiHCMCore;
  G4VPhysicalVolume* fPhysiHCMSkinD;
  G4VPhysicalVolume* fPhysiHCLSkinUR;
  G4VPhysicalVolume* fPhysiHCLCoreR;
  G4VPhysicalVolume* fPhysiHCLSkinDR;
  G4VPhysicalVolume* fPhysiHCLSkinUL;
  G4VPhysicalVolume* fPhysiHCLCoreL;
  G4VPhysicalVolume* fPhysiHCLSkinDL;
  G4VPhysicalVolume* fPhysiPoronVetoLatYPInt;
  G4VPhysicalVolume* fPhysiVetoLatYP;
  G4VPhysicalVolume* fPhysiPoronVetoLatYPExt;
  G4VPhysicalVolume* fPhysiPoronVetoLatYMInt;
  G4VPhysicalVolume* fPhysiVetoLatYM;
  G4VPhysicalVolume* fPhysiPoronVetoLatYMExt;
  G4VPhysicalVolume* fPhysiPoronVetoLatXPInt;
  G4VPhysicalVolume* fPhysiVetoLatXP;
  G4VPhysicalVolume* fPhysiPoronVetoLatXPExt;
  G4VPhysicalVolume* fPhysiPoronVetoLatXMInt;
  G4VPhysicalVolume* fPhysiVetoLatXM;
  G4VPhysicalVolume* fPhysiPoronVetoLatXMExt;
  G4VPhysicalVolume* fPhysiHCSkinLatXPInt;
  G4VPhysicalVolume* fPhysiHCCoreLatXP;   
  G4VPhysicalVolume* fPhysiHCSkinLatXPExt;
  G4VPhysicalVolume* fPhysiHCSkinLatXMInt;
  G4VPhysicalVolume* fPhysiHCCoreLatXM;   
  G4VPhysicalVolume* fPhysiHCSkinLatXMExt;
  G4VPhysicalVolume* fPhysiHCSkinLatYPInt;
  G4VPhysicalVolume* fPhysiHCCoreLatYP;   
  G4VPhysicalVolume* fPhysiHCSkinLatYPExt;
  G4VPhysicalVolume* fPhysiHCSkinLatYMInt;
  G4VPhysicalVolume* fPhysiHCCoreLatYM;   
  G4VPhysicalVolume* fPhysiHCSkinLatYMExt;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif














				
				
				
