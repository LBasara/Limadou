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
/// \file electromagnetic/TestEm3/include/CalorimeterConstructionConfig6.hh
/// \brief Definition of the CalorimeterConstructionConfig6 class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CalorimeterConstructionConfig6_h
#define CalorimeterConstructionConfig6_h 1

#include "globals.hh"
#include <vector>

class G4IntersectionSolid;
class G4SubtractionSolid;
class G4Box;
class G4Trap;
class G4Cons;
class G4UnionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class HEPDSWMaterial;
class G4UniformMagField;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CalorimeterConstructionConfig6 
{
public:
  
   CalorimeterConstructionConfig6();
  ~CalorimeterConstructionConfig6();

  inline void SetVetoMaterial(G4String aMat){vetoMaterial=aMat;}
  inline void SetCaloMaterial(G4String aMat1,G4String aMat2){scintMaterial=aMat1;crystalMaterial=aMat2;}
  inline void SetPoronMaterial(G4String aMat){poronMaterial=aMat;}
  inline void SetCarbonFiberMaterial(G4String aMat){cfiberMaterial=aMat;}
  inline void SetHoneyCombMaterial(G4String aMat){honeycombMaterial=aMat;}
  inline void SetTeflonMaterial(G4String aMat){teflonMaterial=aMat;}
  
  void Builder(G4VPhysicalVolume* motherVolume);
  
private:
  void ComputeObjectsPositioning();

  G4String scintMaterial;
  G4String carbonFiberMaterial;
  G4String poronMaterial;

  G4int    fS1ScintNumber;

  G4double fS1_X;
  G4double fS1_Y;
  G4double fS1_Z;

  G4double fS1SuppBottom_X;
  G4double fS1SuppBottom_Y;
  G4double fS1SuppBottom_Z;
  
  G4double fS1SuppHoleBar_X;
  G4double fS1SuppHoleBar_Y;
  G4double fS1SuppHoleBar_Z;

  G4double fS1SuppBar_X;
  G4double fS1SuppBar_Y;
  G4double fS1SuppBar_Z;

  G4double fS1SuppTopBase_X;
  G4double fS1SuppTopBase_Y;
  G4double fS1SuppTopBase_Z;
  
  G4double fS1SuppTopHole_X;
  G4double fS1SuppTopHole_Y;
  G4double fS1SuppTopHole_Z; 

  G4double fS1SuppThinBar_X;
  G4double fS1SuppThinBar_Y;
  G4double fS1SuppThinBar_Z;

  G4double fS1SuppThinBarBack_X;
  G4double fS1SuppThinBarBack_Y;
  G4double fS1SuppThinBarBack_Z;

  G4double fS1SuppPoron_X;
  G4double fS1SuppPoron_Y;
  G4double fS1SuppPoron_Z;
  
  G4double fS1ScintContainer_X;
  G4double fS1ScintContainer_Y;
  G4double fS1ScintContainer_Z;

  G4double fS1Scint_X;
  G4double fS1Scint_Y;
  G4double fS1Scint_Z;

  G4double fSolidS1SuppFront_X;
  G4double fSolidS1SuppFront_Y;
  G4double fSolidS1SuppFront_Z;

  G4int              fNbOfScintLayers;
  G4int fNbOfReplicatedScintLayers;

  G4String           vetoMaterial;
  G4String           caloMaterial;
  G4String           crystalMaterial;
  G4String           cfiberMaterial;
  G4String           honeycombMaterial;
  G4String           teflonMaterial;
  
  G4double  fCaloLayer_X;
  G4double  fCaloLayer_Y;
  G4double  fCaloLayer_Z;

  G4double  fCaloLastLayer_X;
  G4double  fCaloLastLayer_Y;
  G4double  fCaloLastLayer_Z;


  G4double  fCalo_X;
  G4double  fCalo_Y;
  G4double  fCalo_Z;

  G4double  fScint_X;
  G4double  fScint_Y;
  G4double  fScint_Z;

  G4double  fCrystal_X;
  G4double  fCrystal_Y;
  G4double  fCrystal_Z;

  G4double  fCFBlockContainerExt_X;
  G4double  fCFBlockContainerExt_Y;
  G4double  fCFBlockContainerExt_Z;

  G4double TeflonContainerExt_X;
  G4double TeflonContainerExt_Y;
  G4double TeflonContainerExt_Z;

  G4double fTeflonLYSO_X;
  G4double fTeflonLYSO_Y;
  G4double fTeflonLYSO_Z;

  G4double  fCrystalBlock_X;
  G4double  fCrystalBlock_Y;
  G4double  fCrystalBlock_Z;

  G4double  fCFCrystalPanelUp_Z;

  G4double  fCFCrystalPanelDown_X;
  G4double  fCFCrystalPanelDown_Y;
  G4double  fCFCrystalPanelDown_Z;

  G4double fCFCrystalPanelDown2_X;
  G4double fCFCrystalPanelDown2_Y;
  G4double fCFCrystalPanelDown2_Z;

  G4double fCFCrystalHole_X;
  G4double fCFCrystalHole_Y;
  G4double fCFCrystalHole_Z;
  
  G4double fCylinderRadiusMax1;
  G4double fCylinderRadiusMin1;
  G4double fCylinderRadiusMax2;
  G4double fCylinderRadiusMin2;
  G4double fCylinderZ;

  G4double  fCrystalBlockContainer_X;
  G4double  fCrystalBlockContainer_Y;
  G4double  fCrystalBlockContainer_Z;

  G4double  fCrystalBlockRawContainer_X;
  G4double  fCrystalBlockRawContainer_Y;
  G4double  fCrystalBlockRawContainer_Z;

  G4double  fCFCrystalSideX_X;
  G4double  fCFCrystalSideX_Y;
  G4double  fCFCrystalSideX_Z;

  G4double  fCFCrystalSideYSmall_X;
  G4double  fCFCrystalSideYSmall_Y;
  G4double  fCFCrystalSideYSmall_Z;

  G4double  fCFCrystalSideYBig_X;
  G4double  fCFCrystalSideYBig_Y;
  G4double  fCFCrystalSideYBig_Z;

 
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
  
  G4double fCFSuppAO_X;
  G4double fCFSuppAO_Y;
  G4double fCFSuppAO_Z;

  G4double fCFSuppAV_X;
  G4double fCFSuppAV_Y;
  G4double fCFSuppAV_Z;
  
  G4double fCFSuppPorO_X;
  G4double fCFSuppPorO_Y;
  G4double fCFSuppPorO_Z;  
  
  G4double fCFSuppPorV_X;
  G4double fCFSuppPorV_Y;
  G4double fCFSuppPorV_Z;  


  G4double fCFSuppLat_X;
  G4double fCFSuppLat_Y;
  G4double fCFSuppLat_Z;
    
  G4double fCFSuppFront_X;
  G4double fCFSuppFront_Y;
  G4double fCFSuppFront_Z;
  
  G4double fCFSuppLatA_X;
  G4double fCFSuppLatA_Y;
  G4double fCFSuppLatA_Z;
    
  G4double fCFSuppFrontA_X;
  G4double fCFSuppFrontA_Y;
  G4double fCFSuppFrontA_Z;
  

  G4double fCFLatPO_X ;
  G4double fCFLatPO_Y ;
  G4double fCFLatPO_Z ;

  G4double fCFFrontPO_X ;
  G4double fCFFrontPO_Y ;
  G4double fCFFrontPO_Z ;
  
  G4double fPoronFrontPO_X ; 
  G4double fPoronFrontPO_Y ;
  G4double fPoronFrontPO_Z ;
  
  G4double fPoronLatPO_X ;
  G4double fPoronLatPO_Y ; 
  G4double fPoronLatPO_Z ;
  
  G4double fPoronLat_X;
  G4double fPoronLat_Y;
  G4double fPoronLat_Z;

  G4double fPoronFront_X;
  G4double fPoronFront_Y;
  G4double fPoronFront_Z;

  G4double  fActiveLayer_X;
  G4double  fActiveLayer_Y;
  G4double  fActiveLayer_Z;

  G4double fActiveRectLayer_X;
  G4double fActiveRectLayer_Y;
  G4double fActiveRectLayer_Z;

  G4double fActiveTrapLayerX1;
  G4double fActiveTrapLayerX2;
  G4double fActiveTrapLayerY1;
  G4double fActiveTrapLayerY2;
  G4double fActiveTrapLayerZ;

  G4double fTrapPoronX1;
  G4double fTrapPoronX2;
  G4double fTrapPoronY1;
  G4double fTrapPoronY2;
  G4double fTrapPoronZ;

  G4double fTrapVetoX1;
  G4double fTrapVetoX2;
  G4double fTrapVetoY1;
  G4double fTrapVetoY2;
  G4double fTrapVetoZ;

  G4double fTrapCFX1;
  G4double fTrapCFX2;
  G4double fTrapCFY1;
  G4double fTrapCFY2;
  G4double fTrapCFZ;

  G4double fPoronLatX_X;
  G4double fPoronLatX_Y;
  G4double fPoronLatX_Z;

  G4double fVetoLatX_X;
  G4double fVetoLatX_Y;
  G4double fVetoLatX_Z;

  G4double fVetoLatXRight_X;
  G4double fVetoLatXRight_Y;
  G4double fVetoLatXRight_Z;

   G4double fVetoLatXLeft_X;
  G4double fVetoLatXLeft_Y;
  G4double fVetoLatXLeft_Z;

  G4double fPoronLatY_X;
  G4double fPoronLatY_Y;
  G4double fPoronLatY_Z;

  G4double fPoronLatXRight_X;
  G4double fPoronLatXRight_Y;
  G4double fPoronLatXRight_Z;

  G4double fPoronLatXHole_X;
  G4double fPoronLatXHole_Y;
  G4double fPoronLatXHole_Z;

  G4double fPoronLatXHoleLeft_X;
  G4double fPoronLatXHoleLeft_Y;
  G4double fPoronLatXHoleLeft_Z;

  G4double fPoronLatYUp_X;
  G4double fPoronLatYUp_Y;
  G4double fPoronLatYUp_Z;

  G4double fPoronLatYHole_X;
  G4double fPoronLatYHole_Y;
  G4double fPoronLatYHole_Z;

  G4double fPoronLatYHoleDown_X;
  G4double fPoronLatYHoleDown_Y;
  G4double fPoronLatYHoleDown_Z;

  G4double fVetoLatY_X;
  G4double fVetoLatY_Y;
  G4double fVetoLatY_Z;

  G4double fVetoLatYUp_X;
  G4double fVetoLatYUp_Y;
  G4double fVetoLatYUp_Z;

  G4double fVetoLatYDown_X;
  G4double fVetoLatYDown_Y;
  G4double fVetoLatYDown_Z;

  G4double fCFVetoLatX_X;
  G4double fCFVetoLatX_Y;
  G4double fCFVetoLatX_Z;
  
  G4double fCFVetoLatY_X;
  G4double fCFVetoLatY_Y;
  G4double fCFVetoLatY_Z;

  G4double fPoronPlateO_X;
  G4double fPoronPlateO_Y;
  G4double fPoronPlateO_Z;

  G4double fPoronPlateV_X;
  G4double fPoronPlateV_Y;
  G4double fPoronPlateV_Z;

  G4double  fExternalPoronSupportA_X;
  G4double  fExternalPoronSupportA_Y;
  G4double  fExternalPoronSupportA_Z;

  G4double  fExternalPoronSupportB_X;
  G4double  fExternalPoronSupportB_Y;
  G4double  fExternalPoronSupportB_Z;

  G4double  fExternalPoronSupportC_X;
  G4double  fExternalPoronSupportC_Y;
  G4double  fExternalPoronSupportC_Z;

  G4double  fExternalPoronSupportD_X;
  G4double  fExternalPoronSupportD_Y;
  G4double  fExternalPoronSupportD_Z;

  G4double  fExternalPoronSupportE_X;
  G4double  fExternalPoronSupportE_Y;
  G4double  fExternalPoronSupportE_Z;

  G4double  fExternalPoronSupport_X;
  
  G4double fSolidCFSuppStepLat_X;
  G4double fSolidCFSuppStepLat_Y;
  G4double fSolidCFSuppStepLat_Z;

  G4double fSolidCFSuppStepFront_X;
  G4double fSolidCFSuppStepFront_Y;
  G4double fSolidCFSuppStepFront_Z;

  G4double ShiftOrigin;
  
  G4double fPhysiS1_X;
  G4double fPhysiS1_Y;
  G4double fPhysiS1_Z;

  G4double fPhysiS1SuppHoleBar_X;
  G4double fPhysiS1SuppHoleBar_Y;
  G4double fPhysiS1SuppHoleBar_Z;
  
  G4double fPhysiS1SuppBar_X;
  G4double fPhysiS1SuppBar_Y;
  G4double fPhysiS1SuppBar_Z;

  G4double fPhysiS1SuppBack_X;
  G4double fPhysiS1SuppBack_Y;
  G4double fPhysiS1SuppBack_Z;

  G4double fPhysiS1SuppThinBar_X;
  G4double fPhysiS1SuppThinBar_Y;
  G4double fPhysiS1SuppThinBar_Z;
  
  G4double fPhysiS1SuppThinBarBack_X;
  G4double fPhysiS1SuppThinBarBack_Y;
  G4double fPhysiS1SuppThinBarBack_Z;

  G4double fPhysiS1SuppFrontM_X;
  G4double fPhysiS1SuppFrontM_Y;
  G4double fPhysiS1SuppFrontM_Z;

  G4double fPhysiS1SuppPoronFrontM_X;
  G4double fPhysiS1SuppPoronFrontM_Y;
  G4double fPhysiS1SuppPoronFrontM_Z;
  
  G4double fPhysiS1ScintContainer_X;
  G4double fPhysiS1ScintContainer_Y;
  G4double fPhysiS1ScintContainer_Z;
  
  G4double  fPhysiScintBox_X;
  G4double  fPhysiScintBox_Y;
  G4double  fPhysiScintBox_Z;

  G4double  fPhysiLastScintLayer_X;
  G4double  fPhysiLastScintLayer_Y;
  G4double  fPhysiLastScintLayer_Z;

  G4double  fPhysiLastSuppPor_X;
  G4double  fPhysiLastSuppPor_Y;
  G4double  fPhysiLastSuppPor_Z;

  G4double  fPhysiCrystalBox_X;
  G4double  fPhysiCrystalBox_Y;
  G4double  fPhysiCrystalBox_Z;

  G4double  fPhysiCFCrystalPanelUp_X;
  G4double  fPhysiCFCrystalPanelUp_Y;
  G4double  fPhysiCFCrystalPanelUp_Z;

  G4double  fPhysiCFCrystalPanelDown_X;
  G4double  fPhysiCFCrystalPanelDown_Y;
  G4double  fPhysiCFCrystalPanelDown_Z;

  G4double  fPhysiCFCrystalPanelDown2_X;
  G4double  fPhysiCFCrystalPanelDown2_Y;
  G4double  fPhysiCFCrystalPanelDown2_Z;

  G4double  fPhysiCFCrystalSideX_X;
  G4double  fPhysiCFCrystalSideX_Y;
  G4double  fPhysiCFCrystalSideX_Z;

  G4double  fPhysiCFCrystalSideYBig_X;
  G4double  fPhysiCFCrystalSideYBig_Y;
  G4double  fPhysiCFCrystalSideYBig_Z;

  G4double  fPhysiCrystalBlockPlaneContainer_X;
  G4double  fPhysiCrystalBlockPlaneContainer_Y;
  G4double  fPhysiCrystalBlockPlaneContainer_Z;

  G4double fPhysiCFSupp_X;
  G4double fPhysiCFSupp_Y;
  G4double fPhysiCFSupp_Z;

  G4double fPhysiCFSuppPor_X;
  G4double fPhysiCFSuppPor_Y;
  G4double fPhysiCFSuppPor_Z;

  G4double fPhysiCFSuppPor2_X;
  G4double fPhysiCFSuppPor2_Y;
  G4double fPhysiCFSuppPor2_Z;

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

  G4double fPhysiPoronFront_X;
  G4double fPhysiPoronFront_Y;
  G4double fPhysiPoronFront_Z;

  G4double fPhysiCFFrontPO_X;
  G4double fPhysiCFFrontPO_Y;
  G4double fPhysiCFFrontPO_Z;

  G4double fPhysiCFLatPO_X;
  G4double fPhysiCFLatPO_Y;
  G4double fPhysiCFLatPO_Z;

  G4double fPhysiPoronLatPO_X;
  G4double fPhysiPoronLatPO_Y;
  G4double fPhysiPoronLatPO_Z;

  G4double fPhysiPoronFrontPO_X;
  G4double fPhysiPoronFrontPO_Y;
  G4double fPhysiPoronFrontPO_Z;

  G4double fPhysiExternalPoronSupport_X;
  G4double fPhysiExternalPoronSupport_Y;
  G4double fPhysiExternalPoronSupport_Z;

  G4double fExternalPoronSupport1_X;
  G4double fExternalPoronSupport1_Y;
  G4double fExternalPoronSupport2_X;
  G4double fExternalPoronSupport3_X;
  G4double fExternalPoronSupport3_Y;
  G4double fExternalPoronSupport4_X;
  G4double fExternalPoronSupport4_Y;
  G4double fExternalPoronSupport4_Z;

  G4double fCFVetoLatXHole_X;
  G4double fCFVetoLatXHole_Y;
  G4double fCFVetoLatXHole_Z;

  G4double fCFVetoLatXHoleRight_X;
  G4double fCFVetoLatXHoleRight_Y;
  G4double fCFVetoLatXHoleRight_Z;

  G4double fCFVetoLatYHole_X;
  G4double fCFVetoLatYHole_Y;
  G4double fCFVetoLatYHole_Z;

  G4double fCFVetoLatYHoleDown_X;
  G4double fCFVetoLatYHoleDown_Y;
  G4double fCFVetoLatYHoleDown_Z;

  G4double fPhysiPoronLatX3_X;
  G4double fPhysiPoronLatX3_Y;
  G4double fPhysiPoronLatX3_Z;
  
  G4double fPhysiPoronLatY3_X;
  G4double fPhysiPoronLatY3_Y;
  G4double fPhysiPoronLatY3_Z;
  
  G4double fPhysiPoronLatX3Int_X;
  G4double fPhysiPoronLatX3Int_Y;
  G4double fPhysiPoronLatX3Int_Z;
  
  G4double fPhysiPoronLatY3Int_X;
  G4double fPhysiPoronLatY3Int_Y;
  G4double fPhysiPoronLatY3Int_Z;

  G4double fPhysiVetoLatX2_X;
  G4double fPhysiVetoLatX2_Y;
  G4double fPhysiVetoLatX2_Z;

  G4double fPhysiVetoLatY2_X;
  G4double fPhysiVetoLatY2_Y;
  G4double fPhysiVetoLatY2_Z;

  G4double fPhysiCFVetoLatX2_X;
  G4double fPhysiCFVetoLatX2_Y;
  G4double fPhysiCFVetoLatX2_Z;

  G4double fPhysiCFVetoLatY2_X;
  G4double fPhysiCFVetoLatY2_Y;
  G4double fPhysiCFVetoLatY2_Z;

  G4double fPoronLat1X_X;
  G4double fPoronLat1X_Z;

  G4double fPoronLat2X_X;
  G4double fPoronLat2X_Z;

  G4double fPoronLat3X_X;

  G4double fVetoLat1X_X;
  G4double fVetoLat1X_Z;

  G4double fVetoLat2X_X;
  G4double fVetoLat2X_Z;

  G4double fCFVetoLat1X_X;
  G4double fCFVetoLat1X_Z;

  G4double fCFVetoLat2X_X;

  G4double fPoronLat1Y_Y;
  G4double fPoronLat1Y_Z;

  G4double fPoronLat2Y_Y;
  G4double fPoronLat2Y_Z;

  G4double fPoronLat3Y_Y;

  G4double fVetoLat1Y_Y;
  G4double fVetoLat1Y_Z;

  G4double fVetoLat2Y_Y;
  G4double fVetoLat2Y_Z;

  G4double fCFVetoLat1Y_Y;
  G4double fCFVetoLat1Y_Z;

  G4double fCFVetoLat2Y_Y;

  G4double fPhysiCFSuppTop_X;
  G4double fPhysiCFSuppTop_Y;
  G4double fPhysiCFSuppTop_Z;
  
  G4double fPhysiPORSuppTop_X;
  G4double fPhysiPORSuppTop_Y;
  G4double fPhysiPORSuppTop_Z;
  
  G4double fPhysiRealTrapCF_Z;
  G4double fPhysiRealTrapCF2_Z;
  G4double fPhysiRealTrapPoron1_Z;
  G4double fPhysiRealTrapPoron2_Z;
  G4double fPhysiRealTrapVeto_Z;
  
  G4double fPhysiPoronPlateOTop_X;
  G4double fPhysiPoronPlateOTop_Y;
  G4double fPhysiPoronPlateOTop_Z;

  G4double fPhysiPoronPlateVTop_X;
  G4double fPhysiPoronPlateVTop_Y;
  G4double fPhysiPoronPlateVTop_Z;
 
  G4double fPhysiPoronPlateOBottom_X;
  G4double fPhysiPoronPlateOBottom_Y;
  G4double fPhysiPoronPlateOBottom_Z;

  G4double fPhysiPoronPlateVBottom_X;
  G4double fPhysiPoronPlateVBottom_Y;
  G4double fPhysiPoronPlateVBottom_Z;
  
  HEPDSWMaterial*      pMaterial;
  
  G4Box* fSolidS1;
  G4Box* fSolidS1SuppHoleBar;
  G4Box* fSolidS1SuppBar;
    
  G4Box* fSolidS1SuppBack;
  G4Box* fSolidS1SuppTopBase;
  G4Box* fSolidS1SuppTopHole;
  G4SubtractionSolid* fSolidS1SuppFrontTemp0;
  G4SubtractionSolid* fSolidS1SuppFrontTemp1;
  G4SubtractionSolid* fSolidS1SuppFrontTemp2;
  G4SubtractionSolid* fSolidS1SuppFront;
  G4Box* fSolidS1SuppThinBar; 
  G4Box* fSolidS1SuppThinBarBack; 
  G4Box* fSolidS1SuppPoronFront; 
  G4Box* fSolidS1ScintContainer;
  G4Box* fSolidS1Scint;

  G4Box* fSolidCFSuppPorO;
  G4Box* fSolidCFSuppPorV;
  G4Box* fSolidCaloBox;
  G4Box* fSolidScintBox;
  G4Box* fSolidCrystalBox;
  G4Trap* fSolidActiveTrapLayer;
  
  G4Trap* fSolidTrapPoronLayer;
  G4Trap* fSolidTrapVetoLayer;
  G4Trap* fSolidTrapCFLayer;

  G4Box* fSolidLayer;
  G4Box* fSolidLastLayer;
  G4Box * fSolidActiveRectLayer;
  G4UnionSolid* fSolidActiveLayer_1;
  G4UnionSolid* fSolidActiveLayer;

  G4UnionSolid* fSolidRealTrapPoronLayer;
  G4UnionSolid* fSolidRealTrapCFLayer;
  G4UnionSolid* fSolidRealTrapVetoLayer;

  G4Box* fSolidExternalPoronSupportA;
  G4Box* fSolidExternalPoronSupportB;
  G4Box* fSolidExternalPoronSupportC;
  G4Box* fSolidExternalPoronSupportD;
  G4Box* fSolidExternalPoronSupportE;
  G4UnionSolid* fSolidExternalPoronSupport1;
  G4UnionSolid* fSolidExternalPoronSupport2;
  G4UnionSolid* fSolidExternalPoronSupport3;
  G4UnionSolid* fSolidExternalPoronSupport;

  G4Box* fSolidCFBlockContainerExt;
  G4Box* fSolidCFBlockContainerInt;
  G4SubtractionSolid* fSolidCFBlockContainer;

  G4SubtractionSolid* fSolidCylinderPanelDown1;
  G4SubtractionSolid* fSolidCylinderPanelDown2;
  G4SubtractionSolid* fSolidCylinderPanelDown3;
  G4SubtractionSolid* fSolidCylinderPanelDown4;
  G4SubtractionSolid* fSolidCylinderPanelDown5;
  G4SubtractionSolid* fSolidCylinderPanelDown6;
  G4SubtractionSolid* fSolidCylinderPanelDown7;
  G4SubtractionSolid* fSolidCylinderPanelDown8;
  G4SubtractionSolid* fSolidPanelDown;
  
  G4SubtractionSolid* fSolidSquarePanelDown1;
  G4SubtractionSolid* fSolidSquarePanelDown2;
  G4SubtractionSolid* fSolidSquarePanelDown3;
  G4SubtractionSolid* fSolidSquarePanelDown4;
  G4SubtractionSolid* fSolidSquarePanelDown5;
  G4SubtractionSolid* fSolidSquarePanelDown6;
  G4SubtractionSolid* fSolidSquarePanelDown7;
  G4SubtractionSolid* fSolidSquarePanelDown8;
  G4SubtractionSolid* fSolidPanelDown2;

  G4Box* fSolidSingleCrystalBlockContainer;
  G4Box* fSolidCrystalActiveBlock;
  G4Box*  fSolidTeflonLYSO;
  G4Box*  fSolidTeflonContainerExt;
  G4SubtractionSolid* fSolidTeflonContainer;
  G4Box* fSolidCFCrystalPanelDown;
  G4Box* fSolidCFCrystalPanelDownTemp;
  G4Box* fSolidCFCrystalPanelDownHole;
  G4Cons* fSolidCylinderPanelDown;
  G4Box* fSolidCrystalBlockContainer;
  G4Box* fSolidCrystalBlockRawContainer;
  G4Box* fSolidCFCrystalSideX;
  G4Box* fSolidCFCrystalSideYBig;
  G4Box* fSolidCFFront;
  G4Box* fSolidCFLat;
  G4Box* fSolidCFFrontPO;
  G4Box* fSolidCFLatPO;
  G4Box* fSolidCFSuppO;
  G4Box* fSolidCFSuppV;
  G4Box* fSolidCFSuppLat;
  G4Box* fSolidCFSuppFront;
  G4Box* fSolidCFSuppLatA;
  G4Box* fSolidCFSuppFrontA;
  G4Box* fSolidCFSuppOA;
  G4Box* fSolidCFSuppVA;
  G4UnionSolid* fSolidCFSuppPoron;
  G4UnionSolid* fSolidCFSuppStepOV;
  G4UnionSolid* fSolidCFSuppStepLat1;
  G4UnionSolid* fSolidCFSuppStepLat2;
  G4UnionSolid* fSolidCFSuppStepFront1;
  G4UnionSolid* fSolidCFSupp;
  G4UnionSolid* fSolidCFSuppStepLat1A;
  G4UnionSolid* fSolidCFSuppStepLat2A;
  G4UnionSolid* fSolidCFSuppStepFront1A;
  G4UnionSolid* fSolidCFSuppA;  
  G4UnionSolid* fSolidCFSuppStepOVA;
  G4Box* fSolidPoronFront;
  G4Box* fSolidPoronLat; 
  G4Box* fSolidPoronFrontPO;
  G4Box* fSolidPoronLatPO;
  G4Box* fSolidPoronLatX;
  G4Box* fSolidPoronLatXRight;
  G4Box* fSolidPoronLatXHole;
  G4Box* fSolidPoronLatXHoleLeft;
  G4Box* fSolidPoronLatYUp;
  G4Box* fSolidPoronLatYHole;
  G4Box* fSolidPoronLatYHoleDown;


  G4UnionSolid* fSolidPoronLatX_1;
  G4UnionSolid* fSolidPoronLatX_2;
  G4UnionSolid* fSolidPoronLatX_3;

  G4UnionSolid* fSolidVetoLatX_1;
  G4UnionSolid* fSolidVetoLatX_2;
  
  G4UnionSolid* fSolidPoronLatY_1;
  G4UnionSolid* fSolidPoronLatY_2;
  G4UnionSolid* fSolidPoronLatY_3;

  G4UnionSolid* fSolidVetoLatY_1;
  G4UnionSolid* fSolidVetoLatY_2;
  
  G4UnionSolid* fSolidCFVetoLatX_1;
  G4UnionSolid* fSolidCFVetoLatX_2;
  
  G4UnionSolid* fSolidCFVetoLatY_1;
  G4UnionSolid* fSolidCFVetoLatY_2;

  G4Box* fSolidCFVetoLatX;
  G4Box* fSolidCFVetoLatXHole;
  G4Box* fSolidCFVetoLatXHoleRight;
  
  G4Box* fSolidCFVetoLatY;
  G4Box* fSolidCFVetoLatYHole;
  G4Box* fSolidCFVetoLatYHoleDown;
  
  G4Box* fSolidVetoLatX;
  G4Box* fSolidVetoLatXRight;
  G4Box* fSolidVetoLatXLeft;
  G4Box* fSolidPoronLatY;	   
  G4Box* fSolidVetoLatY;
  G4Box* fSolidVetoLatYUp;
  G4Box* fSolidVetoLatYDown;
  G4Box* fSolidPoronPlateO;
  G4Box* fSolidPoronPlateV;


  G4LogicalVolume* fLogicS1;
  G4LogicalVolume* fLogicS1SuppBack;
  G4LogicalVolume* fLogicS1SuppHoleBar;
  G4LogicalVolume* fLogicS1SuppBar;
  G4LogicalVolume* fLogicS1SuppFront;
  G4LogicalVolume* fLogicS1SuppThinBar;
  G4LogicalVolume* fLogicS1SuppThinBarBack;
  G4LogicalVolume* fLogicS1SuppPoronFront;
  G4LogicalVolume* fLogicS1ScintContainerP;
  G4LogicalVolume* fLogicS1ScintContainerM;
  G4LogicalVolume* fLogicS1Scint;


  G4LogicalVolume* fLogicCFSuppPoron;
  G4LogicalVolume* fLogicCaloBox;
  G4LogicalVolume* fLogicScintBox;
  G4LogicalVolume* fLogicCrystalBox;
  G4LogicalVolume* fLogicScintLayer;
  G4LogicalVolume* fLogicScintActiveLayer;
  G4LogicalVolume* fLogicCrystalBlockContainer;
  G4LogicalVolume* fLogicCrystalBlockPlaneContainer;
  G4LogicalVolume* fLogicCrystalBlockRawContainer;
  G4LogicalVolume* fLogicCrystalActiveBlock;
  G4LogicalVolume* fLogicTeflonLYSO;
  G4LogicalVolume* fLogicTeflonContainer;
  G4LogicalVolume* fLogicCFCrystalPanelDown;
  G4LogicalVolume* fLogicCFCrystalPanelDown2;
  G4LogicalVolume* fLogicCFBlockContainer;
  G4LogicalVolume* fLogicCFCrystalSideX;
  G4LogicalVolume* fLogicCFCrystalSideYBig;
  G4LogicalVolume* fLogicCFFront;
  G4LogicalVolume* fLogicCFLat; 
  G4LogicalVolume* fLogicCFFrontPO;
  G4LogicalVolume* fLogicCFLatPO;
  G4LogicalVolume* fLogicCFSupp;
  G4LogicalVolume* fLogicCFSuppTop;
  G4LogicalVolume* fLogicPoronLat;
  G4LogicalVolume* fLogicPoronFront; 
  G4LogicalVolume* fLogicPoronLatPO;
  G4LogicalVolume* fLogicPoronFrontPO;

  G4LogicalVolume*  fLogicPoronLatX;
  G4LogicalVolume*  fLogicPoronLatY;

  G4LogicalVolume*  fLogicVetoLatX;
  G4LogicalVolume*  fLogicVetoLatX2;
  G4LogicalVolume*  fLogicVetoLatY;
  G4LogicalVolume*  fLogicVetoLatY2;
  
  G4LogicalVolume*  fLogicCFVetoLatX;
  G4LogicalVolume*  fLogicCFVetoLatY;

  G4LogicalVolume*  fLogicExternalPoronSupport;
  G4LogicalVolume*  fLogicRealTrapPoronLayer;
  G4LogicalVolume*  fLogicRealTrapVetoLayer;
  G4LogicalVolume*  fLogicRealTrapCFLayer;
  G4LogicalVolume*  fLogicPoronPlateO;
  G4LogicalVolume*  fLogicPoronPlateV;
  G4LogicalVolume* fLogicLastScintLayer;

  G4VPhysicalVolume* fPhysiS1;
  G4VPhysicalVolume* fPhysiS1SuppBack;
  G4VPhysicalVolume* fPhysiS1SuppHoleBar;
  G4VPhysicalVolume* fPhysiS1SuppBar1;
  G4VPhysicalVolume* fPhysiS1SuppBar2;
  G4VPhysicalVolume* fPhysiS1SuppThinBar1;
  G4VPhysicalVolume* fPhysiS1SuppThinBarBack1;
  G4VPhysicalVolume* fPhysiS1SuppFrontM;
  G4VPhysicalVolume* fPhysiS1SuppFrontPoronM;
  G4VPhysicalVolume* fPhysiS1SuppThinBar2;
  G4VPhysicalVolume* fPhysiS1SuppThinBarBack2;
  G4VPhysicalVolume* fPhysiS1ScintContainerP;
  G4VPhysicalVolume* fPhysiS1ScintContainerM;
  G4VPhysicalVolume* fPhysiS1ScintP;
  G4VPhysicalVolume* fPhysiS1ScintM;

  G4VPhysicalVolume* fPhysiCaloBox;
  G4VPhysicalVolume* fPhysiScintBox;
  G4VPhysicalVolume* fPhysiScintLayer;
  G4VPhysicalVolume* fPhysiScintCFSupp;
  G4VPhysicalVolume* fPhysiScintActiveLayer;
  G4VPhysicalVolume* fPhysiLastScintLayer;
  G4VPhysicalVolume* fPhysiScintCFFrontP;
  G4VPhysicalVolume* fPhysiScintCFFrontM;
  G4VPhysicalVolume* fPhysiScintCFLatP;
  G4VPhysicalVolume* fPhysiScintCFLatM;
  G4VPhysicalVolume* fPhysiScintPoronLatP;
  G4VPhysicalVolume* fPhysiScintPoronLatM;
  G4VPhysicalVolume* fPhysiScintPoronLatUp;
  G4VPhysicalVolume* fPhysiScintPoronLatDown;
  G4VPhysicalVolume* fPhysiScintCFFrontPPO;
  G4VPhysicalVolume* fPhysiScintCFFrontMPO;
  G4VPhysicalVolume* fPhysiScintCFLatPPO;
  G4VPhysicalVolume* fPhysiScintCFLatMPO;
  G4VPhysicalVolume* fPhysiScintPoronLatPPO;
  G4VPhysicalVolume* fPhysiScintPoronLatMPO;
  G4VPhysicalVolume* fPhysiScintPoronLatUpPO;
  G4VPhysicalVolume* fPhysiScintPoronLatDownPO;
  G4VPhysicalVolume* fPhysiCrystalBox;
  G4VPhysicalVolume* fPhysiCFCrystalPanelDown;
  G4VPhysicalVolume* fPhysiCFCrystalPanelDown2;
  G4VPhysicalVolume* fPhysiCFCrystalSideXP;
  G4VPhysicalVolume* fPhysiCFCrystalSideXM;
  G4VPhysicalVolume* fPhysiCFCrystalSideYBigP;
  G4VPhysicalVolume* fPhysiCFCrystalSideYBigM;
  G4VPhysicalVolume* fPhysiCrystalBlockPlaneContainer;
  G4VPhysicalVolume* fPhysiCrystalBlockRaw;
  G4VPhysicalVolume* fPhysiCrystalActiveBlock;
  G4VPhysicalVolume* fPhysiTeflonLYSO;
  G4VPhysicalVolume* fPhysiTeflonContainer;
  
  G4VPhysicalVolume* fPhysiPoronLatXRight_3;
  G4VPhysicalVolume* fPhysiPoronLatXLeft_3;
  G4VPhysicalVolume* fPhysiPoronLatYUp_3;
  G4VPhysicalVolume* fPhysiPoronLatYDown_3;
  
  G4VPhysicalVolume* fPhysiPoronLatXRightInt_3;
  G4VPhysicalVolume* fPhysiPoronLatXLeftInt_3;
  G4VPhysicalVolume* fPhysiPoronLatYUpInt_3;
  G4VPhysicalVolume* fPhysiPoronLatYDownInt_3;
  
  G4VPhysicalVolume* fPhysiVetoLatXRight_2;
  G4VPhysicalVolume* fPhysiVetoLatXLeft_2;
  G4VPhysicalVolume* fPhysiVetoLatYUp_2;
  G4VPhysicalVolume* fPhysiVetoLatYDown_2;

  G4VPhysicalVolume* fPhysiCFVetoLatXRight_2;
  G4VPhysicalVolume* fPhysiCFVetoLatXLeft_2;
  G4VPhysicalVolume* fPhysiCFVetoLatYUp_2;
  G4VPhysicalVolume* fPhysiCFVetoLatYDown_2;

  G4VPhysicalVolume* fPhysiRealTrapPoronLayer1;
  G4VPhysicalVolume* fPhysiRealTrapPoronLayer2;
  G4VPhysicalVolume* fPhysiRealTrapVetoLayer;
  G4VPhysicalVolume* fPhysiRealTrapCFLayer;
  G4VPhysicalVolume* fPhysiRealTrapCFLayer2;
  
  G4VPhysicalVolume* fPhysiExternalPoronSupport1;
  G4VPhysicalVolume* fPhysiExternalPoronSupport2;
  G4VPhysicalVolume* fPhysiExternalPoronSupport3;
  G4VPhysicalVolume* fPhysiExternalPoronSupport4;
    
  G4VPhysicalVolume* fPhysiPoronPlateO1_Top;
  G4VPhysicalVolume* fPhysiPoronPlateV1_Top;
  G4VPhysicalVolume* fPhysiPoronPlateO2_Top;
  G4VPhysicalVolume* fPhysiPoronPlateV2_Top;
  G4VPhysicalVolume* fPhysiPoronPlateO1_Bottom;
  G4VPhysicalVolume* fPhysiPoronPlateV1_Bottom;
  G4VPhysicalVolume* fPhysiPoronPlateO2_Bottom;
  G4VPhysicalVolume* fPhysiPoronPlateV2_Bottom;

  G4VPhysicalVolume* fPhysiCFSuppTop;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif














				
				
				
