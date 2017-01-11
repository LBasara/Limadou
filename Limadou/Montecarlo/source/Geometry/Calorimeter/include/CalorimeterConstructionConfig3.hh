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
/// \file electromagnetic/TestEm3/include/CalorimeterConstructionConfig3.hh
/// \brief Definition of the CalorimeterConstructionConfig3 class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CalorimeterConstructionConfig3_h
#define CalorimeterConstructionConfig3_h 1

#include "globals.hh"

class G4IntersectionSolid;
class G4Box;
class G4UnionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class HEPDSWMaterial;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CalorimeterConstructionConfig3 
{
public:
  
   CalorimeterConstructionConfig3();
  ~CalorimeterConstructionConfig3();

  inline void SetVetoMaterial(G4String aMat){vetoMaterial=aMat;}
  inline void SetCaloMaterial(G4String aMat){caloMaterial=aMat;}
  inline void SetCaloMaterial(G4String aMat1,G4String aMat2){scintMaterial=aMat1;crystalMaterial=aMat2;}
  inline void SetPoronMaterial(G4String aMat){poronMaterial=aMat;}
  inline void SetCarbonFiberMaterial(G4String aMat){cfiberMaterial=aMat;}
  inline void SetHoneyCombMaterial(G4String aMat){honeycombMaterial=aMat;}

  void SetNumberOfCrystalLayer(G4int aVal);
  void Builder(G4VPhysicalVolume* motherVolume);
  
private:
  void ComputeObjectsPositioning();

  G4int              fNbOfCaloLayers;
  G4int              fNbOfCrystalLayers;
  G4int              fNbOfScintLayers;

  G4String           vetoMaterial;
  G4String           caloMaterial;
  G4String           crystalMaterial;
  G4String           scintMaterial;
  G4String           poronMaterial;
  G4String           cfiberMaterial;
  G4String           honeycombMaterial;

  G4double  fCaloLayer_X;
  G4double  fCaloLayer_Y;
  G4double  fCaloLayer_Z;

  G4double  fCalo_X;
  G4double  fCalo_Y;
  G4double  fCalo_Z;

  G4double  fScint_X;
  G4double  fScint_Y;
  G4double  fScint_Z;

  G4double  fCrystal_X;
  G4double  fCrystal_Y;
  G4double  fCrystal_Z;

  G4double fCFLat_X;
  G4double fCFLat_Y;
  G4double fCFLat_Z;

  G4double fCFFront_X;
  G4double fCFFront_Y;
  G4double fCFFront_Z;

  G4double fCFSuppO_X;
  G4double fCFSuppO_Y;
  G4double fCFSuppO_Z;

  G4double fCFSuppV_X;
  G4double fCFSuppV_Y;
  G4double fCFSuppV_Z;

  G4double fCFSuppLat_X;
  G4double fCFSuppLat_Y;
  G4double fCFSuppLat_Z;
    
  G4double fCFSuppFront_X;
  G4double fCFSuppFront_Y;
  G4double fCFSuppFront_Z;
    
  G4double fPoronLat_X;
  G4double fPoronLat_Y;
  G4double fPoronLat_Z;

  G4double  fActiveLayer_X;
  G4double  fActiveLayer_Y;
  G4double  fActiveLayer_Z;

  G4double fPoronLatX_X;
  G4double fPoronLatX_Y;
  G4double fPoronLatX_Z;

  G4double fVetoLatX_X;
  G4double fVetoLatX_Y;
  G4double fVetoLatX_Z;

  G4double fPoronLatY_X;
  G4double fPoronLatY_Y;
  G4double fPoronLatY_Z;

  G4double fVetoLatY_X;
  G4double fVetoLatY_Y;
  G4double fVetoLatY_Z;

  G4double fPoronBottom_X;
  G4double fPoronBottom_Y;
  G4double fPoronBottom_Z;

  G4double fPoronBottomLatStick_X;
  G4double fPoronBottomLatStick_Y;
  G4double fPoronBottomLatStick_Z;

  G4double fVetoBottom_X;
  G4double fVetoBottom_Y;
  G4double fVetoBottom_Z;

  G4double fCFBottomLat_X;
  G4double fCFBottomLat_Y;
  G4double fCFBottomLat_Z;

  G4double fCFBottomFront_X;
  G4double fCFBottomFront_Y;
  G4double fCFBottomFront_Z;

  G4double fHoneyCombLatY_X;
  G4double fHoneyCombLatY_Y;
  G4double fHoneyCombLatY_Z;

  G4double fHoneyCombSkinLatY_X;
  G4double fHoneyCombSkinLatY_Y;
  G4double fHoneyCombSkinLatY_Z;

  G4double fHoneyCombLatX_X;
  G4double fHoneyCombLatX_Y;
  G4double fHoneyCombLatX_Z;

  G4double fHoneyCombSkinLatX_X;
  G4double fHoneyCombSkinLatX_Y;
  G4double fHoneyCombSkinLatX_Z;

  G4double fHoneyCombBottom_X;
  G4double fHoneyCombBottom_Y;
  G4double fHoneyCombBottom_Z;

  G4double fHoneyCombSkinBottom_X;
  G4double fHoneyCombSkinBottom_Y;
  G4double fHoneyCombSkinBottom_Z;

  G4double fSolidCFSuppStepLat_X;
  G4double fSolidCFSuppStepLat_Y;
  G4double fSolidCFSuppStepLat_Z;

  G4double fSolidCFSuppStepFront_X;
  G4double fSolidCFSuppStepFront_Y;
  G4double fSolidCFSuppStepFront_Z;

  G4double  fPhysiScintBox_X;
  G4double  fPhysiScintBox_Y;
  G4double  fPhysiScintBox_Z;

  G4double  fPhysiCrystalBox_X;
  G4double  fPhysiCrystalBox_Y;
  G4double  fPhysiCrystalBox_Z;

  G4double fPhysiCFSupp_X;
  G4double fPhysiCFSupp_Y;
  G4double fPhysiCFSupp_Z;

  G4double fPhysiActiveLayer_X;
  G4double fPhysiActiveLayer_Y;
  G4double fPhysiActiveLayer_Z;

  G4double fPhysiCFFront_X;
  G4double fPhysiCFFront_Y;
  G4double fPhysiCFFront_Z;

  G4double fPhysiCFLat_X;
  G4double fPhysiCFLat_Y;
  G4double fPhysiCFLat_Z;

  G4double fPhysiPoronLat_X;
  G4double fPhysiPoronLat_Y;
  G4double fPhysiPoronLat_Z;

  G4double fPhysiPoronLatXInt_X;
  G4double fPhysiPoronLatXInt_Y;
  G4double fPhysiPoronLatXInt_Z;

  G4double fPhysiPoronLatXExt_X;
  G4double fPhysiPoronLatXExt_Y;
  G4double fPhysiPoronLatXExt_Z;

  G4double fPhysiPoronLatYInt_X;
  G4double fPhysiPoronLatYInt_Y;
  G4double fPhysiPoronLatYInt_Z;

  G4double fPhysiPoronLatYExt_X;
  G4double fPhysiPoronLatYExt_Y;
  G4double fPhysiPoronLatYExt_Z;

  G4double fPhysiVetoLatX_X;
  G4double fPhysiVetoLatX_Y;
  G4double fPhysiVetoLatX_Z;

  G4double fPhysiVetoLatY_X;
  G4double fPhysiVetoLatY_Y;
  G4double fPhysiVetoLatY_Z;

  G4double fPhysiPoronBottomU_X;
  G4double fPhysiPoronBottomU_Y;
  G4double fPhysiPoronBottomU_Z;

  G4double fPhysiPoronBottomD_X;
  G4double fPhysiPoronBottomD_Y;
  G4double fPhysiPoronBottomD_Z;

  G4double fPhysiPoronBottomLatStick_X;
  G4double fPhysiPoronBottomLatStick_Y;
  G4double fPhysiPoronBottomLatStick_Z;

  G4double fPhysiVetoBottom_X;
  G4double fPhysiVetoBottom_Y;
  G4double fPhysiVetoBottom_Z;

  G4double fPhysiCFBottomLat_X;
  G4double fPhysiCFBottomLat_Y;
  G4double fPhysiCFBottomLat_Z;

  G4double fPhysiCFBottomFront_X;
  G4double fPhysiCFBottomFront_Y;
  G4double fPhysiCFBottomFront_Z;

  G4double fPhysiHoneyCombSkinLatYInt_X;
  G4double fPhysiHoneyCombSkinLatYInt_Y;
  G4double fPhysiHoneyCombSkinLatYInt_Z;

  G4double fPhysiHoneyCombSkinLatYExt_X;
  G4double fPhysiHoneyCombSkinLatYExt_Y;
  G4double fPhysiHoneyCombSkinLatYExt_Z;

  G4double fPhysiHoneyCombLatY_X;
  G4double fPhysiHoneyCombLatY_Y;
  G4double fPhysiHoneyCombLatY_Z;

  G4double fPhysiHoneyCombSkinLatXInt_X;
  G4double fPhysiHoneyCombSkinLatXInt_Y;
  G4double fPhysiHoneyCombSkinLatXInt_Z;

  G4double fPhysiHoneyCombSkinLatXExt_X;
  G4double fPhysiHoneyCombSkinLatXExt_Y;
  G4double fPhysiHoneyCombSkinLatXExt_Z;

  G4double fPhysiHoneyCombLatX_X;
  G4double fPhysiHoneyCombLatX_Y;
  G4double fPhysiHoneyCombLatX_Z;

  G4double fPhysiHoneyCombSkinBottomU_X;
  G4double fPhysiHoneyCombSkinBottomU_Y;
  G4double fPhysiHoneyCombSkinBottomU_Z;

  G4double fPhysiHoneyCombSkinBottomD_X;
  G4double fPhysiHoneyCombSkinBottomD_Y;
  G4double fPhysiHoneyCombSkinBottomD_Z;

  G4double fPhysiHoneyCombBottom_X;
  G4double fPhysiHoneyCombBottom_Y;
  G4double fPhysiHoneyCombBottom_Z;

  G4double fPhysiCFSuppTop_X;
  G4double fPhysiCFSuppTop_Y;
  G4double fPhysiCFSuppTop_Z;
  
  HEPDSWMaterial*      pMaterial;

  G4Box* fSolidCaloBox;
  G4Box* fSolidScintBox;
  G4Box* fSolidCrystalBox;
  G4Box* fSolidLayer;
  G4Box* fSolidActiveLayer; 
  G4Box* fSolidCFFront;
  G4Box* fSolidCFLat;
  G4Box* fSolidCFSuppO;
  G4Box* fSolidCFSuppV;
  G4Box* fSolidCFSuppLat;
  G4Box* fSolidCFSuppFront;
  G4UnionSolid* fSolidCFSuppStepOV;
  G4UnionSolid* fSolidCFSuppStepLat1;
  G4UnionSolid* fSolidCFSuppStepLat2;
  G4UnionSolid* fSolidCFSuppStepFront1;
  G4UnionSolid* fSolidCFSupp;
  G4Box* fSolidPoronLat;	   
  G4Box* fSolidPoronLatX;	   
  G4Box* fSolidVetoLatX;	   
  G4Box* fSolidPoronLatY;	   
  G4Box* fSolidVetoLatY;	   
  G4Box* fSolidPoronBottom;	   
  G4Box* fSolidPoronBottomLatStick;
  G4Box* fSolidVetoBottom;	   
  G4Box* fSolidCFBottomLat;	   
  G4Box* fSolidCFBottomFront;	   
  G4Box* fSolidHoneyCombLatY;	   
  G4Box* fSolidHoneyCombSkinLatY;  
  G4Box* fSolidHoneyCombLatX;	   
  G4Box* fSolidHoneyCombSkinLatX;  
  G4Box* fSolidHoneyCombBottom;	   
  G4Box* fSolidHoneyCombSkinBottom;

  G4LogicalVolume* fLogicCaloBox;
  G4LogicalVolume* fLogicScintBox;
  G4LogicalVolume* fLogicCrystalBox;
  G4LogicalVolume* fLogicLayer;
  G4LogicalVolume* fLogicScintLayer;
  G4LogicalVolume* fLogicCrystalLayer;
  G4LogicalVolume* fLogicActiveLayer;
  G4LogicalVolume* fLogicScintActiveLayer;
  G4LogicalVolume* fLogicCrystalActiveLayer;
  G4LogicalVolume* fLogicCFFront;
  G4LogicalVolume* fLogicCFLat;
  G4LogicalVolume* fLogicCFSupp;
  G4LogicalVolume* fLogicCFSuppTop;
  G4LogicalVolume* fLogicPoronLat;	   
  G4LogicalVolume* fLogicPoronLatX;	   
  G4LogicalVolume* fLogicVetoLatX;	   
  G4LogicalVolume* fLogicPoronLatY;	   
  G4LogicalVolume* fLogicVetoLatY;	   
  G4LogicalVolume* fLogicPoronBottom;	   
  G4LogicalVolume* fLogicPoronBottomLatStick;
  G4LogicalVolume* fLogicVetoBottom;	   
  G4LogicalVolume* fLogicCFBottomLat;	   
  G4LogicalVolume* fLogicCFBottomFront;	   
  G4LogicalVolume* fLogicHoneyCombLatY;	   
  G4LogicalVolume* fLogicHoneyCombSkinLatY;  
  G4LogicalVolume* fLogicHoneyCombLatX;	   
  G4LogicalVolume* fLogicHoneyCombSkinLatX;  
  G4LogicalVolume* fLogicHoneyCombBottom;	   
  G4LogicalVolume* fLogicHoneyCombSkinBottom;

  G4VPhysicalVolume* fPhysiCaloBox;
  G4VPhysicalVolume* fPhysiLayer;
  G4VPhysicalVolume* fPhysiCFSupp;
  G4VPhysicalVolume* fPhysiActiveLayer;
  G4VPhysicalVolume* fPhysiCFFrontP;
  G4VPhysicalVolume* fPhysiCFFrontM;
  G4VPhysicalVolume* fPhysiCFLatP;
  G4VPhysicalVolume* fPhysiCFLatM;
  G4VPhysicalVolume* fPhysiPoronLatP;
  G4VPhysicalVolume* fPhysiPoronLatM;
  G4VPhysicalVolume* fPhysiScintBox;
  G4VPhysicalVolume* fPhysiScintLayer;
  G4VPhysicalVolume* fPhysiScintCFSupp;
  G4VPhysicalVolume* fPhysiScintActiveLayer;
  G4VPhysicalVolume* fPhysiScintCFFrontP;
  G4VPhysicalVolume* fPhysiScintCFFrontM;
  G4VPhysicalVolume* fPhysiScintCFLatP;
  G4VPhysicalVolume* fPhysiScintCFLatM;
  G4VPhysicalVolume* fPhysiScintPoronLatP;
  G4VPhysicalVolume* fPhysiScintPoronLatM;
  G4VPhysicalVolume* fPhysiCrystalBox;
  G4VPhysicalVolume* fPhysiCrystalLayer;
  G4VPhysicalVolume* fPhysiCrystalCFSupp;
  G4VPhysicalVolume* fPhysiCrystalActiveLayer;
  G4VPhysicalVolume* fPhysiCrystalCFFrontP;
  G4VPhysicalVolume* fPhysiCrystalCFFrontM;
  G4VPhysicalVolume* fPhysiCrystalCFLatP;
  G4VPhysicalVolume* fPhysiCrystalCFLatM;
  G4VPhysicalVolume* fPhysiCrystalPoronLatP;
  G4VPhysicalVolume* fPhysiCrystalPoronLatM;
  G4VPhysicalVolume* fPhysiPoronLatXPInt;
  G4VPhysicalVolume* fPhysiPoronLatXMInt;
  G4VPhysicalVolume* fPhysiPoronLatXPExt;
  G4VPhysicalVolume* fPhysiPoronLatXMExt;
  G4VPhysicalVolume* fPhysiPoronLatYPInt;
  G4VPhysicalVolume* fPhysiPoronLatYMInt;
  G4VPhysicalVolume* fPhysiPoronLatYPExt;
  G4VPhysicalVolume* fPhysiPoronLatYMExt;
  G4VPhysicalVolume* fPhysiVetoLatXP;
  G4VPhysicalVolume* fPhysiVetoLatXM;
  G4VPhysicalVolume* fPhysiVetoLatYP;
  G4VPhysicalVolume* fPhysiVetoLatYM;
  G4VPhysicalVolume* fPhysiPoronBottomU;
  G4VPhysicalVolume* fPhysiPoronBottomD;
  G4VPhysicalVolume* fPhysiPoronBottomLatStickP;
  G4VPhysicalVolume* fPhysiPoronBottomLatStickM;
  G4VPhysicalVolume* fPhysiVetoBottom;
  G4VPhysicalVolume* fPhysiCFBottomLatP;
  G4VPhysicalVolume* fPhysiCFBottomLatM;
  G4VPhysicalVolume* fPhysiCFBottomFrontP;
  G4VPhysicalVolume* fPhysiCFBottomFrontM;
  G4VPhysicalVolume* fPhysiHoneyCombSkinLatYPInt;
  G4VPhysicalVolume* fPhysiHoneyCombLatYP;
  G4VPhysicalVolume* fPhysiHoneyCombSkinLatYPExt;
  G4VPhysicalVolume* fPhysiHoneyCombSkinLatYMInt;
  G4VPhysicalVolume* fPhysiHoneyCombLatYM;
  G4VPhysicalVolume* fPhysiHoneyCombSkinLatYMExt;
  G4VPhysicalVolume* fPhysiHoneyCombSkinLatXPInt;
  G4VPhysicalVolume* fPhysiHoneyCombLatXP;
  G4VPhysicalVolume* fPhysiHoneyCombSkinLatXPExt;
  G4VPhysicalVolume* fPhysiHoneyCombSkinLatXMInt;
  G4VPhysicalVolume* fPhysiHoneyCombLatXM;
  G4VPhysicalVolume* fPhysiHoneyCombSkinLatXMExt;
  G4VPhysicalVolume* fPhysiHoneyCombSkinBottomU;
  G4VPhysicalVolume* fPhysiHoneyCombBottom;
  G4VPhysicalVolume* fPhysiHoneyCombSkinBottomD;
  G4VPhysicalVolume* fPhysiCFSuppTop;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif














				
				
				
