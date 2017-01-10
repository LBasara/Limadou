
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

#include "CalorimeterConstructionConfig6.hh"

#include "G4NistManager.hh"
#include "HEPDSWMaterial.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4IntersectionSolid.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"

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
CalorimeterConstructionConfig6::CalorimeterConstructionConfig6()
  :fSolidS1(0),
   fSolidS1SuppHoleBar(0),  
   fSolidS1SuppBack(0),
   fSolidS1SuppTopBase(0),
   fSolidS1SuppTopHole(0),
   fSolidS1SuppFrontTemp0(0),
   fSolidS1SuppFrontTemp1(0),
   fSolidS1SuppFrontTemp2(0),
   fSolidS1SuppFront(0),
   fSolidS1ScintContainer(0),
   fSolidS1Scint(0),
  
   fSolidCaloBox(0),fSolidScintBox(0),fSolidCrystalBox(0),fSolidLayer(0),fSolidActiveLayer(0),
   fSolidCFBlockContainerExt(0),fSolidCFBlockContainerInt(0),fSolidCFBlockContainer(0),
   fSolidSingleCrystalBlockContainer(0),
   fSolidCrystalActiveBlock(0),fSolidCFCrystalPanelDown(0),
   fSolidCrystalBlockContainer(0),fSolidCrystalBlockRawContainer(0),fSolidCFCrystalSideX(0),
   fSolidCFCrystalSideYBig(0),
   fSolidCFFront(0),fSolidCFLat(0),fSolidCFSuppO(0),fSolidCFSuppV(0),fSolidCFSuppLat(0),fSolidCFSuppFront(0),fSolidCFSuppLatA(0),fSolidCFSuppFrontA(0),
   fSolidCFSuppPoron(0), fSolidCFSuppStepOV(0),fSolidCFSuppStepLat1(0),fSolidCFSuppStepLat2(0),fSolidCFSuppStepFront1(0),fSolidCFSupp(0),fSolidCFSuppStepLat1A(0),
   fSolidCFSuppStepLat2A(0),fSolidCFSuppStepFront1A(0),fSolidCFSuppA(0), fSolidPoronFront(0), fSolidPoronLat(0), fSolidPoronFrontPO(0),fSolidPoronLatPO(0),
   fSolidPoronLatX(0),fSolidVetoLatX(0),fSolidPoronLatY(0),fSolidVetoLatY(0),
   // fSolidHoneyCombLatY(0),
   //fSolidHoneyCombSkinLatY(0),fSolidHoneyCombLatX(0),fSolidHoneyCombSkinLatX(0),fSolidHoneyCombBottom(0),fSolidHoneyCombSkinBottom(0),
   
   fLogicS1(0),
   fLogicS1SuppBack(0),
   fLogicS1SuppHoleBar(0),
   fLogicS1SuppFront(0),
   fLogicS1SuppThinBar(0),
   fLogicS1ScintContainerP(0),
   fLogicS1ScintContainerM(0),
   fLogicS1Scint(0),

   fLogicCaloBox(0),fLogicScintBox(0),fLogicCrystalBox(0),fLogicScintLayer(0),fLogicScintActiveLayer(0),
   fLogicCrystalBlockContainer(0),fLogicCrystalBlockPlaneContainer(0),fLogicCrystalBlockRawContainer(0),fLogicCrystalActiveBlock(0),
   fLogicCFCrystalPanelDown(0),fLogicCFBlockContainer(0),
   fLogicCFCrystalSideX(0),fLogicCFCrystalSideYBig(0),
   fLogicCFFront(0),fLogicCFLat(0),fLogicCFFrontPO(0),fLogicCFLatPO(0),fLogicCFSupp(0),fLogicCFSuppTop(0),fLogicPoronLat(0), fLogicPoronFront(0),fLogicPoronLatPO(0),fLogicPoronFrontPO(0), fLogicPoronLatX(0),	   
   fLogicPoronLatY(0),fLogicVetoLatX(0),fLogicVetoLatY(0),  

   fPhysiS1(0),
   fPhysiS1SuppBack(0),
   fPhysiS1SuppThinBar1(0),
   fPhysiS1SuppFrontM(0),
   fPhysiS1SuppThinBar2(0),
   fPhysiS1ScintContainerP(0),
   fPhysiS1ScintContainerM(0),
   fPhysiS1ScintP(0),
   fPhysiS1ScintM(0),
   
   fPhysiCaloBox(0),fPhysiScintBox(0),fPhysiScintLayer(0),fPhysiScintCFSupp(0),fPhysiScintActiveLayer(0),fPhysiScintCFFrontP(0),fPhysiScintCFFrontM(0),fPhysiScintCFLatP(0),fPhysiScintCFLatM(0),fPhysiScintPoronLatP(0),fPhysiScintPoronLatM(0),fPhysiScintPoronLatUp(0),fPhysiScintPoronLatDown(0),
   fPhysiCrystalBox(0),fPhysiCFCrystalPanelDown(0),fPhysiCFCrystalSideXP(0),fPhysiCFCrystalSideXM(0),fPhysiCFCrystalSideYBigP(0),fPhysiCFCrystalSideYBigM(0),fPhysiCrystalBlockPlaneContainer(0),fPhysiCrystalBlockRaw(0),
   fPhysiCrystalActiveBlock(0),
   fPhysiCFSuppTop(0)
{

  pMaterial           = new HEPDSWMaterial();


  // S1 Scintillator
  fS1ScintNumber     = 6;
  
  fS1_X = 242*mm;//224*mm;
  fS1_Y = 241.12;//223.72*mm;
  fS1_Z = 10.5*mm;
 
  fS1SuppBottom_X  = 200*mm;
  fS1SuppBottom_Y  = 241.12*mm;  
  fS1SuppBottom_Z  = 0.5*mm;
  
  fS1SuppHoleBar_X  = 242*mm;
  fS1SuppHoleBar_Y  = 6*mm;
  fS1SuppHoleBar_Z  = 5*mm; 

  fS1SuppBar_X  = 242*mm;
  fS1SuppBar_Y  = 27.56*mm;
  fS1SuppBar_Z  = 6.14*mm;

  fS1SuppTopBase_X  = 6*mm;
  fS1SuppTopBase_Y  = 241.12;//202*mm;
  fS1SuppTopBase_Z  = 3.86*mm;
  
  fS1SuppTopHole_Y  = 90*mm;
  fS1SuppTopHole_X  = 2*mm;
  fS1SuppTopHole_Z  = 3*mm;
  
  fS1SuppThinBar_Y  = 186*mm;
  fS1SuppThinBar_X  = 20*mm;         // piccola approsimazione. In realtà la sovrapposizione tra fS1SuppThinBar e fS1SuppBar è più complessa 
  fS1SuppThinBar_Z  = .5*mm;         // e agli angoli lo spessore è un po diverso
  
  fS1SuppThinBarBack_Y  = 241.12*mm;
  fS1SuppThinBarBack_X  = 20*mm;
  fS1SuppThinBarBack_Z  = .5*mm;

  fS1SuppPoron_Y = 185*mm;
  fS1SuppPoron_X = 6*mm;
  fS1SuppPoron_Z = 1.14*mm; 
 
  fS1ScintContainer_Y = 90*mm;
  fS1ScintContainer_X = 242*mm;
  fS1ScintContainer_Z = 5*mm;

  fS1Scint_Y = 30*mm;
  fS1Scint_X = 242*mm; //210*mm
  fS1Scint_Z = 5*mm;
  /////////////////////////////////////////////


  // Calorimeter
  fNbOfScintLayers   = 16;
  fNbOfReplicatedScintLayers   = 15;
    
  fCaloLayer_X    = 177*mm;
  fCaloLayer_Y    = 177*mm;
  fCaloLayer_Z    = 14.78*mm;

  fCaloLastLayer_X    = 177*mm;
  fCaloLastLayer_Y    = 177*mm;
  fCaloLastLayer_Z    = 13.64*mm;

  fCalo_X = 177*mm;
  fCalo_Y = 177*mm; 
  fCalo_Z = 280.34*mm;   // dimensione dell'intero calorimetro
  
  fScint_X = 177*mm;
  fScint_Y = 177*mm;
  fScint_Z = fCaloLayer_Z*fNbOfReplicatedScintLayers;

  fCrystal_X = 177*mm;
  fCrystal_Y = 177*mm;
  fCrystal_Z = 45*mm;

  fCFBlockContainerExt_X = 50.8*mm;
  fCFBlockContainerExt_Y = 50.8*mm;
  fCFBlockContainerExt_Z = 40.0*mm;
 
  TeflonContainerExt_X = 48.2*mm;
  TeflonContainerExt_Y = 48.2*mm;
  TeflonContainerExt_Z = 40.2*mm;
  
  fTeflonLYSO_X = 48.2*mm;
  fTeflonLYSO_Y = 48.2*mm;
  fTeflonLYSO_Z = 40.2*mm;
  
  // cubi di LYSO
  fCrystalBlock_X = 48*mm;
  fCrystalBlock_Y = 48*mm;
  fCrystalBlock_Z = 40*mm;

  fCFCrystalPanelDown_X = 177*mm;
  fCFCrystalPanelDown_Y = 177*mm;
  fCFCrystalPanelDown_Z = 3*mm;

  fCylinderRadiusMax1 = 9*mm;
  fCylinderRadiusMin1 = 0*mm;
  fCylinderRadiusMax2 = 9*mm;
  fCylinderRadiusMin2 = 0*mm;
  fCylinderZ = 3*mm;
  
  fCFCrystalPanelDown2_X = 158.4*mm;
  fCFCrystalPanelDown2_Y = 158.4*mm;
  fCFCrystalPanelDown2_Z = 5. *mm;

  fCFCrystalHole_X = 44.8*mm;
  fCFCrystalHole_Y = 44.8*mm;
  fCFCrystalHole_Z = 5*mm;
  
  fCrystalBlockContainer_X = 152.4*mm;
  fCrystalBlockContainer_Y = 152.4*mm;
  fCrystalBlockContainer_Z = 40*mm;

  fCrystalBlockRawContainer_X = 152.4*mm;
  fCrystalBlockRawContainer_Y = 50.8*mm;
  fCrystalBlockRawContainer_Z = 40*mm;

  // costola laterale LYSO
  fCFCrystalSideX_X = 10.9*mm;
  fCFCrystalSideX_Y = 122*mm;
  fCFCrystalSideX_Z = 42*mm;

  // costola superiore LYSO
  fCFCrystalSideYSmall_X = 0.1*mm;
  fCFCrystalSideYSmall_Y = 0.1*mm;
  fCFCrystalSideYSmall_Z = 0.1*mm;

  fCFCrystalSideYBig_X = 122*mm;
  fCFCrystalSideYBig_Y = 10.9*mm;
  fCFCrystalSideYBig_Z = 42*mm;


  fCFLat_X = 10*mm;            //
  fCFLat_Y = 60*mm;            //     
  fCFLat_Z = 11.14*mm;         //
                               // costole laterali e frontali delle croci
  fCFFront_X = 60*mm;          //
  fCFFront_Y = 10*mm;          // 
  fCFFront_Z = 11.14*mm;       //

  
  fPoronFront_X = 60*mm;
  fPoronFront_Y = 3.18*mm;
  fPoronFront_Z = 11.14*mm;
  
  fPoronLat_X = 3.18*mm;
  fPoronLat_Y = 60*mm;       // poron sotto blocchi croce destra-sinistra
  fPoronLat_Z = 11.14*mm;
  
  fCFLatPO_X =  10*mm;
  fCFLatPO_Y =  60*mm;
  fCFLatPO_Z =  1.14*mm;

  fCFFrontPO_X =  60*mm;
  fCFFrontPO_Y =  10*mm;
  fCFFrontPO_Z =  1.14*mm;

  
  fPoronFrontPO_X =  60*mm; 
  fPoronFrontPO_Y =  3.18*mm;
  fPoronFrontPO_Z = 1.14*mm;
  
  fPoronLatPO_X =  3.18*mm;
  fPoronLatPO_Y = 60*mm; 
  fPoronLatPO_Z =  1.14*mm;
  
  fCFSuppO_X = 150.64*mm;    //
  fCFSuppO_Y = 3*mm;         //  braccio N-S della croce carbonio
  fCFSuppO_Z = 2.5*mm;       //
  
  fCFSuppV_X = 3*mm;          //
  fCFSuppV_Y = 150.64*mm;     //  braccio E-O della croce carbonio
  fCFSuppV_Z = 2.5*mm;        //
  
  fCFSuppPorO_X = 149*mm;    //
  fCFSuppPorO_Y = 7*mm;      // braccio N-S croce poron
  fCFSuppPorO_Z = 1.14*mm;   //

  
  fCFSuppPorV_X = 7*mm;      //
  fCFSuppPorV_Y = 149*mm;    // braccio E-O croce poron
  fCFSuppPorV_Z = 1.14*mm;   //
  

  fCFSuppAO_X = 150.64*mm;       //
  fCFSuppAO_Y = 3*mm;            //
  fCFSuppAO_Z = 2.9*mm;          //
                                 //
  fCFSuppAV_X = 3*mm;            //  
  fCFSuppAV_Y = 150.64*mm;       //  croce carbonio trigger S2
  fCFSuppAV_Z = 2.9*mm;          //
                                 //
  fCFSuppLatA_X = 13.18*mm;      //
  fCFSuppLatA_Y = 14*mm;         //
  fCFSuppLatA_Z = 2.9*mm;        //
                                 //
  fCFSuppFrontA_X = 14*mm;       //
  fCFSuppFrontA_Y = 13.18*mm;    //
  fCFSuppFrontA_Z = 2.9*mm;      //
                                 //

  fCFSuppLat_X = 13.18*mm;   //
  fCFSuppLat_Y = 60*mm;      // blocchi rettangolari attaccati alla croce a destra e sinistra
  fCFSuppLat_Z = 2.5*mm;     //
    
  fCFSuppFront_X = 60*mm;      //
  fCFSuppFront_Y = 13.18*mm;   // blocchi rettangolari attaccati alla croce sopra e sotto
  fCFSuppFront_Z = 2.5*mm;     //

  fActiveLayer_X = 15*cm;
  fActiveLayer_Y = 15*cm;
  fActiveLayer_Z = 1*cm;

  fActiveRectLayer_X = 15*cm;
  fActiveRectLayer_Y = 12*cm;  // parte rettangolare centrale del piano di scintillatore
  fActiveRectLayer_Z = 1*cm;
  
  fActiveTrapLayerX1 = 150*mm;
  fActiveTrapLayerX2 = 120*mm;
  fActiveTrapLayerY1 = 10*mm;   // parte trapezoidale del piano di scintillatore
  fActiveTrapLayerY2 = 10*mm;
  fActiveTrapLayerZ = 15*mm;

  fPoronPlateO_X = 147.17*mm; 
  fPoronPlateO_Y = 24.06*mm;
  fPoronPlateO_Z = 3.18*mm;
  
  fPoronPlateV_X = 24.06*mm; 
  fPoronPlateV_Y = 173.76*mm;
  fPoronPlateV_Z = 3.18*mm;


  ////////////// TRAPEZI FINALI BOTTOM VETO  ///////////////

  fTrapPoronX1 = 312.48*mm ;
  fTrapPoronX2 = 71.05*mm;
  fTrapPoronY1 = 3.18*mm;
  fTrapPoronY2 = 3.18*mm; 
  fTrapPoronZ =  241.43/2.*mm;

  fTrapVetoX1 = 309.20*mm ;
  fTrapVetoX2 = 67.77*mm;
  fTrapVetoY1 = 5*mm;
  fTrapVetoY2 = 5*mm; 
  fTrapVetoZ =  241.43/2.*mm;

  fTrapCFX1 = 312.48*mm ;
  fTrapCFX2 = 71.05*mm;
  fTrapCFY1 = 0.7*mm;
  fTrapCFY2 = 0.7*mm; 
  fTrapCFZ =  241.43/2.*mm;
  ///////////////////////////////////////////////////////////
 
 
  fPoronLatX_X = 118.34*mm;
  fPoronLatX_Y = 9.53*mm;
  fPoronLatX_Z = 303.64*mm;

  fPoronLatXRight_X = 34.83*mm;  //parte laterale destra del poron esterno attorno al calorimetro
  fPoronLatXRight_Y = 9.53*mm;
  fPoronLatXRight_Z = 280.14*mm;

  fPoronLatXHole_X = 20*mm;  //parte laterale sinistra col BUCO  del poron esterno attorno al calorimetro
  fPoronLatXHole_Y = 9.53*mm;
  fPoronLatXHole_Z = 290.64*mm;

  fPoronLatXHoleLeft_X = 5.83*mm;  //parte laterale a sinistra del BUCO  del poron esterno attorno al calorimetro
  fPoronLatXHoleLeft_Y = 9.53*mm;
  fPoronLatXHoleLeft_Z = 303.64*mm;



  //////////////////////////////////   VETO ////////////////////////////
  fVetoLatX_X = 118.34*mm;
  fVetoLatX_Y = 5*mm;
  fVetoLatX_Z = 303.64*mm;

  fVetoLatXRight_X = 42*mm;
  fVetoLatXRight_Y = 5*mm;
  fVetoLatXRight_Z = 280.14*mm;

  fVetoLatXLeft_X = 33*mm;
  fVetoLatXLeft_Y = 5*mm;
  fVetoLatXLeft_Z = 290.64*mm;

  fVetoLatY_X = 5*mm;
  fVetoLatY_Y = 128.4*mm;
  fVetoLatY_Z = 303.64*mm;

  fVetoLatYUp_X = 5*mm;
  fVetoLatYUp_Y = 51*mm;
  fVetoLatYUp_Z = 280.14*mm;

  fVetoLatYDown_X = 5*mm;
  fVetoLatYDown_Y = 43*mm;
  fVetoLatYDown_Z = 290.64*mm;
  /////////////////////////////////////////////////////////////////////////////////


  fPoronLatY_X = 9.53*mm;
  fPoronLatY_Y = 128.4*mm;
  fPoronLatY_Z = 303.64*mm;

  fPoronLatYUp_X = 9.53*mm;  //parte superiore/inferiore del poron esterno attorno al calorimetro
  fPoronLatYUp_Y = 38.83*mm;
  fPoronLatYUp_Z = 280.14*mm;

  fPoronLatYHole_X = 9.53*mm;  //parte laterale sinistra col BUCO  del poron esterno attorno al calorimetro
  fPoronLatYHole_Y = 20*mm;
  fPoronLatYHole_Z = 290.64*mm;

  fPoronLatYHoleDown_X = 9.53*mm;  //parte laterale a sinistra del BUCO  del poron esterno attorno al calorimetro
  fPoronLatYHoleDown_Y = 10.83*mm;
  fPoronLatYHoleDown_Z = 303.64*mm;



  fCFVetoLatX_X = 152.17*mm;   //parte frontale carbonio attorno al calorimetro
  fCFVetoLatX_Y = 1*mm;
  fCFVetoLatX_Z = 310*mm;
  
  fCFVetoLatXHole_X = 20*mm;  
  fCFVetoLatXHole_Y = 1*mm;
  fCFVetoLatXHole_Z = 293.82*mm;

  fCFVetoLatXHoleRight_X = 4.83*mm;  
  fCFVetoLatXHoleRight_Y = 1*mm;
  fCFVetoLatXHoleRight_Z = 310*mm;

  fCFVetoLatY_X = 1*mm;  //parte laterale carbonio attorno al calorimetro
  fCFVetoLatY_Y = 157.70*mm;
  fCFVetoLatY_Z = 310*mm;

  fCFVetoLatYHole_X = 1*mm;  
  fCFVetoLatYHole_Y = 20*mm;
  fCFVetoLatYHole_Z = 293.82*mm;

  fCFVetoLatYHoleDown_X = 1*mm;  
  fCFVetoLatYHoleDown_Y = 1.3*mm;
  fCFVetoLatYHoleDown_Z = 310*mm;


  ////////////////////////////////
  //   fHoneyCombLatY_X = 10*mm;
  //   fHoneyCombLatY_Y = 192.72*mm;
  //   fHoneyCombLatY_Z = 322.24*mm;

  //   fHoneyCombSkinLatY_X = 1*mm;
  //   fHoneyCombSkinLatY_Y = 192.72*mm;
  //   fHoneyCombSkinLatY_Z = 322.24*mm;

  //   fHoneyCombLatX_X = 199.72*mm;
  //   fHoneyCombLatX_Y = 10*mm;
  //   fHoneyCombLatX_Z = 322.24*mm;

  //   fHoneyCombSkinLatX_X = 199.72*mm;
  //   fHoneyCombSkinLatX_Y = 1*mm;
  //   fHoneyCombSkinLatX_Z = 322.24*mm;

  //   fHoneyCombBottom_X = 300*mm;
  //   fHoneyCombBottom_Y = 300*mm;
  //   fHoneyCombBottom_Z = 10*mm;

  //   fHoneyCombSkinBottom_X = 300*mm;
  //   fHoneyCombSkinBottom_Y = 300*mm;
  //   fHoneyCombSkinBottom_Z = 1*mm;



  /*############### Costole poron angoli calorimetro ################### */
  fExternalPoronSupportA_X = 2.5*mm;
  fExternalPoronSupportA_Y = 9.53*mm;
  fExternalPoronSupportA_Z = 291.64*mm;

  fExternalPoronSupportB_X = 4.67*mm;
  fExternalPoronSupportB_Y = 9.53*mm;
  fExternalPoronSupportB_Z = 303.64*mm;

  fExternalPoronSupportC_X = 2.36*mm;
  fExternalPoronSupportC_Y = 14.53*mm;
  fExternalPoronSupportC_Z = 303.64*mm;

  fExternalPoronSupportD_X = 5*mm;
  fExternalPoronSupportD_Y = 2.36*mm;
  fExternalPoronSupportD_Z = 303.64*mm;

  fExternalPoronSupportE_X = 9.53*mm;
  fExternalPoronSupportE_Y = 14.53*mm;
  fExternalPoronSupportE_Z = 303.64*mm;

  
  fExternalPoronSupport_X = fExternalPoronSupportA_X + fExternalPoronSupportB_X + fExternalPoronSupportC_X + fExternalPoronSupportD_X + fExternalPoronSupportE_X; 
  ///////////////////////////////////////////////////////////////




  ComputeObjectsPositioning();
  
  // materials
  vetoMaterial      = "Scintillator";
  caloMaterial      = "Scintillator";
  scintMaterial     = "Scintillator";
  crystalMaterial   = "LYSO";

  poronMaterial     = "Poron";
  cfiberMaterial    = "CarbonFiber";
  honeycombMaterial = "nomex";
  teflonMaterial    = "Teflon";
  /*
    poronMaterial     = "Galactic";
    cfiberMaterial    = "Galactic";
    honeycombMaterial = "Galactic";
  */  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterConstructionConfig6::~CalorimeterConstructionConfig6()
{
  if (pMaterial) delete pMaterial;
}

void CalorimeterConstructionConfig6::ComputeObjectsPositioning(){

 
  ShiftOrigin = fCalo_Z/2. + 36.16*mm + 0.7*mm;
  
  fPhysiS1SuppBack_Y = 0;//fS1SuppBottom_X/2.-18.860*mm-fS1SuppBottomHole_X/2.;
  fPhysiS1SuppBack_X = 0;
  fPhysiS1SuppBack_Z = -fS1_Z/2. +fS1SuppBottom_Z/2. ; //fS1SuppBottom_Z/2.-0.5*mm-fS1SuppBottomHole_Z/2.;

  fSolidS1SuppFront_Y = 3*mm + fS1SuppTopHole_Y/2;//fS1SuppTopBase_Y/2.-8*mm-fS1SuppTopHole_Y/2.;
  fSolidS1SuppFront_X = fS1SuppTopBase_X/2.-fS1SuppTopHole_X/2.;
  fSolidS1SuppFront_Z = fS1SuppTopBase_Z/2.-fS1SuppTopHole_Z/2.;

  fPhysiS1_Y = 0;
  fPhysiS1_X = 0;
  fPhysiS1_Z = fCalo_Z/2. + fS1_Z/2. + fCFSuppAV_Z +fCFSuppPorO_Z +7.22*mm;//   147.8*mm+18.28*mm+fS1_Z/2.;

  fPhysiS1SuppHoleBar_Y  = 0;
  fPhysiS1SuppHoleBar_X  = 0;
  fPhysiS1SuppHoleBar_Z  = -fS1_Z/2. + fS1SuppHoleBar_Z/2. +fS1SuppBottom_Z ;

  fPhysiS1SuppBar_Y  = fS1_Y/2. -fS1SuppBar_Y/2. ;
  fPhysiS1SuppBar_X  = 0;
  fPhysiS1SuppBar_Z  = -fS1_Z/2. + fS1SuppBar_Z/2. +fS1SuppBottom_Z ;

  fPhysiS1SuppThinBar_Y  = 0;
  fPhysiS1SuppThinBar_X  = fS1_X/2. -fS1SuppThinBar_X/2. ;
  fPhysiS1SuppThinBar_Z  = -fS1_Z/2. + fS1SuppThinBar_Z/2. +fS1SuppBottom_Z + fS1ScintContainer_Z  ;
  
  fPhysiS1SuppThinBarBack_Y  = 0;
  fPhysiS1SuppThinBarBack_X  = fS1_X/2. -fS1SuppThinBarBack_X/2. ;
  fPhysiS1SuppThinBarBack_Z  = -fS1_Z/2. + fS1SuppThinBarBack_Z/2.;
  
  fPhysiS1SuppFrontM_Y = 0;
  fPhysiS1SuppFrontM_X = 0;
  fPhysiS1SuppFrontM_Z = fS1_Z/2.-fS1SuppTopBase_Z/2.;

  fPhysiS1SuppPoronFrontM_Y = 0;
  fPhysiS1SuppPoronFrontM_X = 0;
  fPhysiS1SuppPoronFrontM_Z = fS1_Z/2. -fS1SuppTopBase_Z - fS1SuppPoron_Z/2.;
 
  fPhysiS1ScintContainer_Y = 0 + 3*mm + fS1ScintContainer_Y/2.;            //fS1_Y/2.-18.86*mm-fS1ScintContainer_Y/2.;
  fPhysiS1ScintContainer_X = 0;
  fPhysiS1ScintContainer_Z = -fS1_Z/2.+0.5*mm + fS1ScintContainer_Z/2.;

  fPhysiScintBox_X = 0;
  fPhysiScintBox_Y = 0;
  fPhysiScintBox_Z = fCalo_Z/2.-fScint_Z/2.;

  fPhysiLastScintLayer_X = 0;
  fPhysiLastScintLayer_Y = 0;
  fPhysiLastScintLayer_Z = +fCalo_Z/2. - fCaloLastLayer_Z/2.-fScint_Z;
 
  fSolidCFSuppStepLat_X = fCaloLayer_X/2.-fCFLat_X/2. -fPoronLat_X/2.;
  fSolidCFSuppStepLat_Y = 0;
  fSolidCFSuppStepLat_Z = 0;

  fSolidCFSuppStepFront_X = 0;
  fSolidCFSuppStepFront_Y = fCaloLayer_Y/2.-fCFFront_Y/2. -fPoronFront_Y/2.;
  fSolidCFSuppStepFront_Z = 0;
  
  fPhysiCFSuppPor_X = 0;
  fPhysiCFSuppPor_Y = 0;
  fPhysiCFSuppPor_Z = fCaloLayer_Z/2. - fCFSuppPorV_Z/2. -fActiveLayer_Z ;
  
  fPhysiCFSupp_X = 0;
  fPhysiCFSupp_Y = 0;
  fPhysiCFSupp_Z =  fCaloLayer_Z/2. - fCFSuppV_Z/2. -fActiveLayer_Z  - fCFSuppPorV_Z ; 

  fPhysiCFSuppPor2_X = 0;
  fPhysiCFSuppPor2_Y = 0;
  fPhysiCFSuppPor2_Z = fCaloLayer_Z/2. - fCFSuppPorV_Z/2. -fActiveLayer_Z -fCFSuppPorV_Z - fCFSuppV_Z  ;

  fPhysiActiveLayer_X = 0;
  fPhysiActiveLayer_Y = 0;
  fPhysiActiveLayer_Z = fCaloLayer_Z/2. -fActiveLayer_Z/2.; 

  fPhysiCFFront_X = 0;
  fPhysiCFFront_Y = fCaloLayer_Y/2.-fCFFront_Y/2.;
  fPhysiCFFront_Z = fCaloLayer_Z/2. -fActiveLayer_Z/2.-fCFSuppPorV_Z/2.; 

  fPhysiCFLat_X = fCaloLayer_X/2.-fCFLat_X/2.;
  fPhysiCFLat_Y = 0;
  fPhysiCFLat_Z =  fCaloLayer_Z/2. -fActiveLayer_Z/2. -fCFSuppPorV_Z/2.;

  fPhysiPoronLat_X = fCaloLayer_X/2. -fCFLat_X -fPoronLat_X/2.;
  fPhysiPoronLat_Y = 0;
  fPhysiPoronLat_Z = fCaloLayer_Z/2. -fActiveLayer_Z/2. -fCFSuppPorV_Z/2.;
 
  fPhysiPoronFront_X = 0;
  fPhysiPoronFront_Y = fCaloLayer_Y/2. -fCFFront_Y- fPoronFront_Y/2.;
  fPhysiPoronFront_Z =  fCaloLayer_Z/2. -fActiveLayer_Z/2. -fCFSuppPorV_Z/2.;

  //////////////////////////////
  
  fPhysiCFFrontPO_X = 0;
  fPhysiCFFrontPO_Y = fCaloLayer_Y/2.-fCFFront_Y/2.;
  fPhysiCFFrontPO_Z = fCaloLayer_Z/2. - fCFSuppPorV_Z/2. -fActiveLayer_Z -fCFSuppPorV_Z - fCFSuppV_Z  ;

  fPhysiCFLatPO_X = fCaloLayer_X/2.-fCFLat_X/2.;
  fPhysiCFLatPO_Y = 0;
  fPhysiCFLatPO_Z = fCaloLayer_Z/2. - fCFSuppPorV_Z/2. -fActiveLayer_Z -fCFSuppPorV_Z - fCFSuppV_Z  ;


  fPhysiPoronLatPO_X = fCaloLayer_X/2. -fCFLat_X -fPoronLat_X/2.;
  fPhysiPoronLatPO_Y = 0;
  fPhysiPoronLatPO_Z = fCaloLayer_Z/2. - fCFSuppPorV_Z/2. -fActiveLayer_Z -fCFSuppPorV_Z - fCFSuppV_Z  ;
 
  fPhysiPoronFrontPO_X = 0;
  fPhysiPoronFrontPO_Y = fCaloLayer_Y/2. -fCFFront_Y- fPoronFront_Y/2.;
  fPhysiPoronFrontPO_Z = fCaloLayer_Z/2. - fCFSuppPorV_Z/2. -fActiveLayer_Z -fCFSuppPorV_Z - fCFSuppV_Z  ;  

  fPhysiCrystalBox_X = 0;
  fPhysiCrystalBox_Y = 0;
  fPhysiCrystalBox_Z = -fCalo_Z/2. + fCrystal_Z/2.;  
  
  fPhysiCFCrystalPanelUp_X = 0;
  fPhysiCFCrystalPanelUp_Y = 0;
  fPhysiCFCrystalPanelUp_Z = fCrystal_Z/2.-fCFCrystalPanelUp_Z/2.;

  
  fPhysiCFCrystalPanelDown_X = 0;
  fPhysiCFCrystalPanelDown_Y = 0;
  fPhysiCFCrystalPanelDown_Z = -fCrystal_Z/2. + fCFCrystalPanelDown_Z/2.;

  fPhysiCFCrystalPanelDown2_X = 0;
  fPhysiCFCrystalPanelDown2_Y = 0;
  fPhysiCFCrystalPanelDown2_Z = -fCalo_Z/2. - fCFCrystalPanelDown2_Z/2.;

  fPhysiCFCrystalSideX_X = fCrystal_X/2. -fCFCrystalSideX_X/2.;
  fPhysiCFCrystalSideX_Y = 0;
  fPhysiCFCrystalSideX_Z = -fCrystal_Z/2. + fCFCrystalPanelDown_Z + fCFCrystalSideYBig_Z/2.;
  
  fPhysiCFCrystalSideYBig_X = 0;
  fPhysiCFCrystalSideYBig_Y = fCrystal_Y/2.-fCFCrystalSideYBig_Y/2.;
  fPhysiCFCrystalSideYBig_Z = -fCrystal_Z/2. + fCFCrystalPanelDown_Z + fCFCrystalSideYBig_Z/2.;

  fPhysiCrystalBlockPlaneContainer_X = 0;
  fPhysiCrystalBlockPlaneContainer_Y = 0;
  fPhysiCrystalBlockPlaneContainer_Z = -fCrystal_Z/2. + fCFCrystalPanelDown_Z + fCrystalBlockContainer_Z/2.;

  fPhysiLastSuppPor_X = 0;
  fPhysiLastSuppPor_Y = 0;
  fPhysiLastSuppPor_Z = +fCrystal_Z/2. -fCFSuppPorO_Z/2.;


  

  //////////  Bottom VETO
  fPhysiRealTrapCF_Z = -fCalo_Z/2.-23.4*mm-fTrapCFY1/2.;
  fPhysiRealTrapCF2_Z = -fCalo_Z/2.-23.4*mm -fTrapCFY1 -fTrapPoronY1-fTrapVetoY1-fTrapPoronY1 -fTrapCFY1/2.;
  fPhysiRealTrapPoron1_Z = -fCalo_Z/2.-23.4*mm-fTrapCFY1-fTrapPoronY1-fTrapVetoY1-fTrapPoronY1/2. ;
  fPhysiRealTrapPoron2_Z = -fCalo_Z/2.-23.4*mm-fTrapCFY1-fTrapPoronY1/2.;
  fPhysiRealTrapVeto_Z = -fCalo_Z/2.-23.4*mm-fTrapCFY1-fTrapPoronY1-fTrapVetoY1/2.;
  //////////////////////////////////////
 
  fPhysiPoronLatX3_X = -4.5*mm;
  fPhysiPoronLatX3_Y = fCalo_Y/2.+fCFVetoLatX_Y+fPoronLatX_Y+fVetoLatX_Y+fPoronLatX_Y/2.;
  fPhysiPoronLatX3_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;

  fPhysiPoronLatY3_X =  fCalo_X/2.+fCFVetoLatY_X+fPoronLatY_X+fVetoLatY_X+fPoronLatY_X/2.;
  fPhysiPoronLatY3_Y =  -4*mm;
  fPhysiPoronLatY3_Z =  -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;

  fPhysiPoronLatX3Int_X =  -4.5*mm;
  fPhysiPoronLatX3Int_Y =  fCalo_Y/2.+fCFVetoLatX_Y+fPoronLatX_Y/2.;
  fPhysiPoronLatX3Int_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;

  fPhysiPoronLatY3Int_X = fCalo_X/2.+fCFVetoLatY_X+fPoronLatY_X/2.;
  fPhysiPoronLatY3Int_Y =  -4*mm;
  fPhysiPoronLatY3Int_Z =  -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;

  fPhysiVetoLatX2_X = (fVetoLatXRight_X - fVetoLatXLeft_X)/2.;    /// 4.5*mm;
  fPhysiVetoLatX2_Y = fCalo_Y/2. + fCFVetoLatX_Y +  fPoronLatX_Y + fVetoLatX_Y/2. ;
  fPhysiVetoLatX2_Z = -fCalo_Z/2. -15.22*mm +fVetoLatY_Z/2.;

  fPhysiVetoLatY2_X = fCalo_X/2.+ fPoronLatY_X +fVetoLatY_X/2.+fCFVetoLatY_X;
  fPhysiVetoLatY2_Y = (fVetoLatYUp_Y - fVetoLatYDown_Y)/2.; //4*mm;
  fPhysiVetoLatY2_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;

  fPhysiCFVetoLatX2_X = (fCFVetoLatXHole_X + fCFVetoLatXHoleRight_X)/2.;
  fPhysiCFVetoLatX2_Y =  fCalo_Y/2. + fCFVetoLatX_Y/2.;
  fPhysiCFVetoLatX2_Z = -fCalo_Z/2. -18.4*mm +fVetoLatY_Z/2.;

  fPhysiCFVetoLatY2_X =  fCalo_X/2.+fCFVetoLatY_X/2.;
  fPhysiCFVetoLatY2_Y = (fCFVetoLatYHole_Y + fCFVetoLatYHoleDown_Y)/2.; // 10.65*mm;
  fPhysiCFVetoLatY2_Z =   -fCalo_Z/2. -18.4*mm+fVetoLatY_Z/2.;

  fPoronLat1X_X = fPoronLatX_X/2. +  fPoronLatXRight_X/2. ;
  fPoronLat1X_Z = -fPoronLatX_Z/2. + fPoronLatXRight_Z/2. ;
   
  fPoronLat2X_X = -fPoronLatX_X/2. -fPoronLatXHole_X/2. ; //-fPoronLatXRight_X/2.
  fPoronLat2X_Z = fPoronLatX_Z/2. - fPoronLatXHole_Z/2. ;


  fPoronLat3X_X = -fPoronLatX_X/2. -fPoronLatXHole_X -fPoronLatXHoleLeft_X/2.; //-fPoronLatXRight_X/2.


  fVetoLat1X_X = fVetoLatX_X/2. + fVetoLatXRight_X/2.;
  fVetoLat1X_Z = -fVetoLatX_Z/2. + fVetoLatXRight_Z/2.;
  
  fVetoLat2X_X = -fVetoLatX_X/2. - fVetoLatXLeft_X/2.;
  fVetoLat2X_Z = fVetoLatX_Z/2. - fVetoLatXLeft_Z/2.;

  fCFVetoLat1X_X = - fCFVetoLatX_X/2. - fCFVetoLatXHole_X/2.;
  fCFVetoLat1X_Z =  fCFVetoLatX_Z/2. - fCFVetoLatXHole_Z/2.;

  
  fCFVetoLat2X_X = - fCFVetoLatX_X/2. - fCFVetoLatXHole_X -  fCFVetoLatXHoleRight_X/2.;


  fPoronLat1Y_Y = fPoronLatY_Y/2. +  fPoronLatYUp_Y/2. ;
  fPoronLat1Y_Z = -fPoronLatY_Z/2. + fPoronLatYUp_Z/2. ;
   
  fPoronLat2Y_Y = -fPoronLatY_Y/2. -fPoronLatYHole_Y/2. ;
  fPoronLat2Y_Z = fPoronLatY_Z/2. - fPoronLatYHole_Z/2. ;


  fPoronLat3Y_Y = -fPoronLatY_Y/2. -  fPoronLatYHole_Y -  fPoronLatYHoleDown_Y/2.;


  fVetoLat1Y_Y = fVetoLatY_Y/2. + fVetoLatYUp_Y/2.;
  fVetoLat1Y_Z = -fVetoLatY_Z/2. + fVetoLatYUp_Z/2.;
  
  fVetoLat2Y_Y = -fVetoLatY_Y/2. - fVetoLatYDown_Y/2.;
  fVetoLat2Y_Z = fVetoLatY_Z/2. - fVetoLatYDown_Z/2.;

  fCFVetoLat1Y_Y = - fCFVetoLatY_Y/2. - fCFVetoLatYHole_Y/2.;
  fCFVetoLat1Y_Z =  fCFVetoLatY_Z/2. - fCFVetoLatYHole_Z/2.;

  
  fCFVetoLat2Y_Y = - fCFVetoLatY_Y/2. - fCFVetoLatYHole_Y -  fCFVetoLatYHoleDown_Y/2.;

  // Sbarre poron angoli calorimetro
  fPhysiExternalPoronSupport_X = fCalo_X/2.+ 25.06*mm - fExternalPoronSupport_X/2.+7.265*mm;
  fPhysiExternalPoronSupport_Y =  fCalo_Y/2.+fCFVetoLatX_Y+fPoronLatX_Y+fVetoLatX_Y+fPoronLatX_Y - fExternalPoronSupportE_Y/2.;
  fPhysiExternalPoronSupport_Z =  -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;

 
  /*
  fPhysiPoronLatXInt_X = 0*mm;
  fPhysiPoronLatXInt_Y = fCalo_Y/2.+fPoronLatX_Y/2.;
  fPhysiPoronLatXInt_Z = 1.96*mm;

  
  fPhysiPoronLatXExt_X = -4.5*mm;
  fPhysiPoronLatXExt_Y = fCalo_Y/2.+fCFVetoLatX_Y+fPoronLatX_Y+fVetoLatX_Y+fPoronLatX_Y/2.;
  fPhysiPoronLatXExt_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;

  fPhysiPoronLatXExtRight_X = -4.5*mm+fPoronLatX_X/2.+ fPoronLatXRight_X/2.;
  fPhysiPoronLatXExtRight_Y = fCalo_Y/2.+fCFVetoLatX_Y+fPoronLatX_Y+fVetoLatX_Y+fPoronLatX_Y/2.;
  fPhysiPoronLatXExtRight_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.-fPoronLatX_Z/2+fPoronLatXRight_Z/2.;

  fPhysiPoronLatXExtHole_X = -4.5*mm- fPoronLatX_X/2. - fPoronLatXHole_X/2.;
  fPhysiPoronLatXExtHole_Y = fCalo_Y/2.+fCFVetoLatX_Y+fPoronLatX_Y+fVetoLatX_Y+fPoronLatX_Y/2.;
  fPhysiPoronLatXExtHole_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.+fPoronLatX_Z/2-fPoronLatXHole_Z/2.;

  fPhysiPoronLatXExtHoleLeft_X =-4.5*mm - fPoronLatX_X/2. - fPoronLatXHole_X-fPoronLatXHoleLeft_X/2.;
  fPhysiPoronLatXExtHoleLeft_Y = fCalo_Y/2.+fCFVetoLatX_Y+fPoronLatX_Y+fVetoLatX_Y+fPoronLatX_Y/2.;
  fPhysiPoronLatXExtHoleLeft_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;
  */

  


  ////////  poron attorno Calorimetro interni - frontale
  /*
  fPhysiPoronLatXInt_X = -4.5*mm;
  fPhysiPoronLatXInt_Y = fCalo_Y/2.+fCFVetoLatX_Y+fPoronLatX_Y/2.;
  fPhysiPoronLatXInt_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;

  fPhysiPoronLatXIntRight_X = -4.5*mm+fPoronLatX_X/2.+ fPoronLatXRight_X/2.;
  fPhysiPoronLatXIntRight_Y = fCalo_Y/2.+fPoronLatX_Y/2.+fCFVetoLatX_Y;
  fPhysiPoronLatXIntRight_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.-fPoronLatX_Z/2+fPoronLatXRight_Z/2.;

  fPhysiPoronLatXIntHole_X = -4.5*mm- fPoronLatX_X/2. - fPoronLatXHole_X/2.;
  fPhysiPoronLatXIntHole_Y = fCalo_Y/2.+fCFVetoLatX_Y+fPoronLatX_Y/2.;
  fPhysiPoronLatXIntHole_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.+fPoronLatX_Z/2-fPoronLatXHole_Z/2.;

  fPhysiPoronLatXIntHoleLeft_X =-4.5*mm - fPoronLatX_X/2. - fPoronLatXHole_X-fPoronLatXHoleLeft_X/2.;
  fPhysiPoronLatXIntHoleLeft_Y = fCalo_Y/2.+fCFVetoLatX_Y+fPoronLatX_Y/2.;
  fPhysiPoronLatXIntHoleLeft_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;
     
  fPhysiPoronLatYInt_X = fCalo_X/2.+fPoronLatY_X/2.;
  fPhysiPoronLatYInt_Y = 0;
  fPhysiPoronLatYInt_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;
    
  fPhysiPoronLatYExt_X = fCalo_X/2.+fCFVetoLatY_X+fPoronLatY_X+fVetoLatY_X+fPoronLatY_X/2.;
  fPhysiPoronLatYExt_Y = -4*mm;
  fPhysiPoronLatYExt_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;

  fPhysiPoronLatYExtUp_X = fCalo_X/2.+fCFVetoLatY_X+fPoronLatY_X+fVetoLatY_X+fPoronLatY_X/2.;
  fPhysiPoronLatYExtUp_Y = fPoronLatY_Y/2.+fPoronLatYUp_Y/2.-4*mm;
  fPhysiPoronLatYExtUp_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.-fPoronLatY_Z/2+fPoronLatYUp_Z/2.;

  fPhysiPoronLatYExtHole_X = fCalo_X/2.+fCFVetoLatY_X+fPoronLatY_X+fVetoLatY_X+fPoronLatY_X/2.;
  fPhysiPoronLatYExtHole_Y = -fPoronLatY_Y/2-fPoronLatYHole_Y/2.-4*mm;
  fPhysiPoronLatYExtHole_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.+fPoronLatY_Z/2-fPoronLatYHole_Z/2.;

  fPhysiPoronLatYExtHoleDown_X = fCalo_X/2.+fCFVetoLatY_X+fPoronLatY_X+fVetoLatY_X+fPoronLatY_X/2.;
  fPhysiPoronLatYExtHoleDown_Y = -fPoronLatY_Y/2-fPoronLatYHole_Y - fPoronLatYHoleDown_Y/2.-4*mm;
  fPhysiPoronLatYExtHoleDown_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;
  */
  ///////////////////// poron attorno Calorimetro interni - laterale

  /*
  fPhysiPoronLatYInt_X = fCalo_X/2.+fCFVetoLatY_X+fPoronLatY_X/2.;
  fPhysiPoronLatYInt_Y = -4*mm;
  fPhysiPoronLatYInt_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;

  fPhysiPoronLatYIntUp_X = fCalo_X/2.+fCFVetoLatY_X+fPoronLatY_X/2.;
  fPhysiPoronLatYIntUp_Y = fPoronLatY_Y/2.+fPoronLatYUp_Y/2.-4*mm;
  fPhysiPoronLatYIntUp_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.-fPoronLatY_Z/2+fPoronLatYUp_Z/2.;

  fPhysiPoronLatYIntHole_X = fCalo_X/2.+fCFVetoLatY_X+fPoronLatY_X/2.;
  fPhysiPoronLatYIntHole_Y = -fPoronLatY_Y/2-fPoronLatYHole_Y/2.-4*mm;
  fPhysiPoronLatYIntHole_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.+fPoronLatY_Z/2-fPoronLatYHole_Z/2.;

  fPhysiPoronLatYIntHoleDown_X = fCalo_X/2.+fCFVetoLatY_X+fPoronLatY_X/2.;
  fPhysiPoronLatYIntHoleDown_Y = -fPoronLatY_Y/2-fPoronLatYHole_Y - fPoronLatYHoleDown_Y/2.-4*mm;
  fPhysiPoronLatYIntHoleDown_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;

  ///////////////////////////// carbonio attorno al calorimetro

  fPhysiCFVetoLatX_X = +12.415*mm;
  fPhysiCFVetoLatX_Y = fCalo_Y/2. + fCFVetoLatX_Y/2.;
  fPhysiCFVetoLatX_Z = -fCalo_Z/2. -18.4*mm +fVetoLatY_Z/2.;

  fPhysiCFVetoLatXHole_X =  +12.415*mm - fCFVetoLatX_X/2. - fCFVetoLatXHole_X/2.;
  fPhysiCFVetoLatXHole_Y =  fCalo_Y/2. + fCFVetoLatX_Y/2.;
  fPhysiCFVetoLatXHole_Z = -fCalo_Z/2. -18.4*mm +fVetoLatY_Z/2.+ fCFVetoLatX_Z/2. - fCFVetoLatXHole_Z/2. ;


  fPhysiCFVetoLatXHoleRight_X = +12.415*mm - fCFVetoLatX_X/2. - fCFVetoLatXHole_X - fCFVetoLatXHoleRight_X/2.; 
  fPhysiCFVetoLatXHoleRight_Y =  fCalo_Y/2. + fCFVetoLatX_Y/2.;
  fPhysiCFVetoLatXHoleRight_Z = -fCalo_Z/2. -18.4*mm+fVetoLatY_Z/2.;

  fPhysiCFVetoLatY_X = fCalo_X/2.+fCFVetoLatY_X/2.;
  fPhysiCFVetoLatY_Y = 10.65*mm;
  fPhysiCFVetoLatY_Z =  -fCalo_Z/2. -18.4*mm+fVetoLatY_Z/2.;

  fPhysiCFVetoLatYHole_X = fCalo_X/2.+fCFVetoLatY_X/2.; 
  fPhysiCFVetoLatYHole_Y = 10.65*mm - fCFVetoLatY_Y/2. - fCFVetoLatYHole_Y/2.;
  fPhysiCFVetoLatYHole_Z = -fCalo_Z/2. -18.4*mm +fVetoLatY_Z/2.+ fCFVetoLatY_Z/2. -  fCFVetoLatYHole_Z/2.;


  fPhysiCFVetoLatYHoleDown_X = fCalo_X/2.+fCFVetoLatY_X/2.;
  fPhysiCFVetoLatYHoleDown_Y = 10.65*mm - fCFVetoLatY_Y/2. - fCFVetoLatYHole_Y - fCFVetoLatYHoleDown_Y/2.;
  fPhysiCFVetoLatYHoleDown_Z = -fCalo_Z/2. -18.4*mm+fVetoLatY_Z/2.;

  ////////////////////////////

  fPhysiVetoLatX_X = -4.5*mm;
  fPhysiVetoLatX_Y = fCalo_Y/2.+fPoronLatX_Y+fVetoLatX_Y/2.+fCFVetoLatX_Y;
  fPhysiVetoLatX_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;

  fPhysiVetoLatXRight_X = +fVetoLatX_X/2.+fVetoLatXRight_X/2.-4.5*mm;
  fPhysiVetoLatXRight_Y = fCalo_Y/2.+fPoronLatX_Y+fVetoLatX_Y/2.+fCFVetoLatX_Y;
  fPhysiVetoLatXRight_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.-fVetoLatX_Z/2.+fVetoLatXRight_Z/2.;

  fPhysiVetoLatXLeft_X = -fVetoLatX_X/2.-fVetoLatXLeft_X/2.-4.5*mm ;
  fPhysiVetoLatXLeft_Y = fCalo_Y/2.+fPoronLatX_Y+fVetoLatXLeft_Y/2.+fCFVetoLatX_Y;
  fPhysiVetoLatXLeft_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.-fVetoLatXLeft_Z/2.+fVetoLatX_Z/2.;

  fPhysiVetoLatY_X = fCalo_X/2.+fPoronLatY_X+fVetoLatY_X/2.+fCFVetoLatY_X;
  fPhysiVetoLatY_Y = -4*mm;
  fPhysiVetoLatY_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.;

  fPhysiVetoLatYUp_X = fCalo_X/2.+fPoronLatY_X+fVetoLatY_X/2.+fCFVetoLatY_X;
  fPhysiVetoLatYUp_Y = -4*mm+fVetoLatY_Y/2.+fVetoLatYUp_Y/2.;
  fPhysiVetoLatYUp_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.-fVetoLatY_Z/2+fVetoLatYUp_Z/2;

  fPhysiVetoLatYDown_X = fCalo_X/2.+fPoronLatY_X+fVetoLatY_X/2.+fCFVetoLatY_X;
  fPhysiVetoLatYDown_Y = -4*mm-fVetoLatY_Y/2.-fVetoLatYDown_Y/2.;
  fPhysiVetoLatYDown_Z = -fCalo_Z/2. -15.22*mm + fVetoLatY_Z/2.+fVetoLatY_Z/2-fVetoLatYDown_Z/2;;
  */
  ///////////////////////////////////////////////

  
  

  
  fPhysiPoronPlateOTop_X = -fCFCrystalPanelDown_X/2. -1*mm +fPoronPlateO_X/2.;
  fPhysiPoronPlateOTop_Y = fCFCrystalPanelDown_Y/2. + 1*mm + fPoronPlateO_Y/2.;
  fPhysiPoronPlateOTop_Z = +fCalo_Z/2. +8.08*mm +fPoronPlateO_Z/2.;
  
  fPhysiPoronPlateVTop_X = -fCFCrystalPanelDown_X/2. -1*mm -fPoronPlateV_X/2.;
  fPhysiPoronPlateVTop_Y = fCFCrystalPanelDown_Y/2. +25.06*mm - fPoronPlateV_Y/2.;
  fPhysiPoronPlateVTop_Z = +fCalo_Z/2. +8.08*mm +fPoronPlateO_Z/2.;



  fPhysiPoronPlateOBottom_X = +fCFCrystalPanelDown_X/2. +1*mm -fPoronPlateO_X/2.;
  fPhysiPoronPlateOBottom_Y = fCFCrystalPanelDown_Y/2. + 1*mm + fPoronPlateO_Y/2.;
  fPhysiPoronPlateOBottom_Z = -fCalo_Z/2. -15.22*mm -fPoronPlateO_Z/2. ;
  
  fPhysiPoronPlateVBottom_X = fCFCrystalPanelDown_X/2. +1*mm +fPoronPlateV_X/2.;
  fPhysiPoronPlateVBottom_Y = fCFCrystalPanelDown_Y/2. +25.06*mm - fPoronPlateV_Y/2.;
  fPhysiPoronPlateVBottom_Z = -fCalo_Z/2. -15.22*mm -fPoronPlateV_Z/2.;
  

  fPhysiPORSuppTop_X = 0;
  fPhysiPORSuppTop_Y = 0;
  fPhysiPORSuppTop_Z = fCalo_Z/2.+ fCFSuppPorO_Z/2.;


  fPhysiCFSuppTop_X = 0;
  fPhysiCFSuppTop_Y = 0;
  fPhysiCFSuppTop_Z = fCalo_Z/2.+fCFSuppAO_Z/2. + fCFSuppPorO_Z;  // aggiungere spessore croce poron

  fExternalPoronSupport1_X = -fExternalPoronSupportE_X/2.- fExternalPoronSupportD_X/2.;
  fExternalPoronSupport1_Y = +fExternalPoronSupportE_Y/2.- fExternalPoronSupportD_Y/2.;

  fExternalPoronSupport2_X = -fExternalPoronSupportE_X/2.- fExternalPoronSupportD_X - fExternalPoronSupportC_X/2.;

  fExternalPoronSupport3_X = -fExternalPoronSupportE_X/2.- fExternalPoronSupportD_X - fExternalPoronSupportC_X - fExternalPoronSupportB_X/2;
  fExternalPoronSupport3_Y = +fExternalPoronSupportC_Y/2.- fExternalPoronSupportB_Y/2;

  fExternalPoronSupport4_X = -fExternalPoronSupportE_X/2.- fExternalPoronSupportD_X - fExternalPoronSupportC_X - fExternalPoronSupportB_X - fExternalPoronSupportA_X/2.;
  fExternalPoronSupport4_Y = +fExternalPoronSupportC_Y/2.- fExternalPoronSupportA_Y/2;
  fExternalPoronSupport4_Z = +fExternalPoronSupportB_Z/2.- +fExternalPoronSupportA_Z/2.;



}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterConstructionConfig6::Builder(G4VPhysicalVolume* motherVolume)
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
  G4Material* scintLayerMat   = pMaterial->GetMaterial(scintMaterial);
  G4Material* crystalLayerMat = pMaterial->GetMaterial(crystalMaterial);
  G4Material* vetoLayerMat    = pMaterial->GetMaterial(vetoMaterial); 
  G4Material* cfMat           = pMaterial->GetMaterial(cfiberMaterial);
  G4Material* porMat          = pMaterial->GetMaterial(poronMaterial);
  //  G4Material* hcMat           = pMaterial->GetMaterial(honeycombMaterial);
  G4Material* teflon          = pMaterial->GetMaterial(teflonMaterial);



  // S1 scintillator
  fSolidS1 = new G4Box("S1",fS1_X/2.,fS1_Y/2.,fS1_Z/2.);

  fSolidS1SuppHoleBar = new G4Box("CFSupportHoleBar",fS1SuppHoleBar_X/2.,fS1SuppHoleBar_Y/2.,fS1SuppHoleBar_Z/2.);

  fSolidS1SuppBar = new G4Box("CFSupportBar",fS1SuppBar_X/2.,fS1SuppBar_Y/2.,fS1SuppBar_Z/2.);

  fSolidS1SuppThinBar = new G4Box("CFSupportThinBar", fS1SuppThinBar_X/2., fS1SuppThinBar_Y/2.,fS1SuppThinBar_Z/2.);

  fSolidS1SuppThinBarBack = new G4Box("CFSupportThinBar", fS1SuppThinBarBack_X/2., fS1SuppThinBarBack_Y/2.,fS1SuppThinBarBack_Z/2.);

  fSolidS1SuppBack = new G4Box("CFSupportBack", fS1SuppBottom_X/2., fS1SuppBottom_Y/2., fS1SuppBottom_Z/2.);

  fSolidS1SuppTopBase = new G4Box("CFSupport",fS1SuppTopBase_X/2.,fS1SuppTopBase_Y/2.,fS1SuppTopBase_Z/2.);

  fSolidS1SuppTopHole = new G4Box("CFSupport",fS1SuppTopHole_X/2.,fS1SuppTopHole_Y/2.,fS1SuppTopHole_Z/2.);

  fSolidS1SuppFrontTemp0 = new G4SubtractionSolid("CFSupport",fSolidS1SuppTopBase,fSolidS1SuppTopHole,0,
						  G4ThreeVector(fSolidS1SuppFront_X,fSolidS1SuppFront_Y,fSolidS1SuppFront_Z));

  fSolidS1SuppFrontTemp1 = new G4SubtractionSolid("CFSupport",fSolidS1SuppFrontTemp0,fSolidS1SuppTopHole,0,
						  G4ThreeVector(-fSolidS1SuppFront_X,fSolidS1SuppFront_Y,fSolidS1SuppFront_Z));

  fSolidS1SuppFrontTemp2 = new G4SubtractionSolid("CFSupport",fSolidS1SuppFrontTemp1,fSolidS1SuppTopHole,0,
						  G4ThreeVector(fSolidS1SuppFront_X,-fSolidS1SuppFront_Y,fSolidS1SuppFront_Z));

  fSolidS1SuppFront = new G4SubtractionSolid("CFSupport",fSolidS1SuppFrontTemp2,fSolidS1SuppTopHole,0,
					     G4ThreeVector(-fSolidS1SuppFront_X,-fSolidS1SuppFront_Y,fSolidS1SuppFront_Z));

  fSolidS1SuppPoronFront = new G4Box("CFSupportPoron", fS1SuppPoron_X/2., fS1SuppPoron_Y/2., fS1SuppPoron_Z/2.);
					
 
  fSolidS1ScintContainer = new G4Box("ScintContainer",fS1ScintContainer_X/2.,fS1ScintContainer_Y/2.,fS1ScintContainer_Z/2.);

  fSolidS1Scint = new G4Box("S1Scint",fS1Scint_X/2.,fS1Scint_Y/2.,fS1Scint_Z/2.);
  ///////////////////////////////////////////////////////



  // Calorimeter
  fSolidCaloBox = new G4Box("Calorimeter",fCalo_X/2.,fCalo_Y/2.,fCalo_Z/2.);
   
  fSolidScintBox = new G4Box("CalorimeterScint",fScint_X/2.,fScint_Y/2.,fScint_Z/2.);

  fSolidCrystalBox = new G4Box("CalorimeterCrystal",fCrystal_X/2.,fCrystal_Y/2.,fCrystal_Z/2.);

  fSolidLayer = new G4Box("Layer",fCaloLayer_X/2.,fCaloLayer_Y/2.,fCaloLayer_Z/2.);
  
  fSolidLastLayer = new G4Box("LastActiveLayer", fCaloLastLayer_X/2., fCaloLastLayer_Y/2., fCaloLastLayer_Z/2.);

  G4RotationMatrix * _rot90 = new G4RotationMatrix;
  G4RotationMatrix * _rot45 = new G4RotationMatrix;
  G4RotationMatrix * _rot180 = new G4RotationMatrix;
  G4RotationMatrix * _rotminus90 = new G4RotationMatrix;
  G4RotationMatrix * _rot180Z = new G4RotationMatrix;
  G4RotationMatrix * _rot180Y = new G4RotationMatrix;
  G4RotationMatrix * _rot180ZY = new G4RotationMatrix;

  _rot90->rotateX(90 * deg);
  _rot45->rotateX(90 * deg);
  _rot45->rotateY(-45 * deg);
  _rotminus90->rotateX(-90 * deg);
  _rot180->rotateX(180 * deg);
  _rot180Z->rotateZ(180 * deg);
  _rot180Y->rotateY(180 * deg);
  _rot180ZY->rotateY(180 * deg);
  _rot180ZY->rotateZ(180 * deg);
  
  
  fSolidActiveTrapLayer = new G4Trap ("TrapezoidalLayer", fActiveTrapLayerX1/2., fActiveTrapLayerX2/2., fActiveTrapLayerY1/2., fActiveTrapLayerY2/2., fActiveTrapLayerZ/2.); 

  fSolidActiveRectLayer = new G4Box("RectangularLayer", fActiveRectLayer_X/2.,fActiveRectLayer_Y/2.,fActiveRectLayer_Z/2.);  
  
  fSolidActiveLayer_1 = new G4UnionSolid("ActiveLayer1",  fSolidActiveRectLayer,  fSolidActiveTrapLayer, _rot90, G4ThreeVector(0, fActiveLayer_Y/2. -fActiveTrapLayerZ/2. ,0));     

  fSolidActiveLayer = new G4UnionSolid("ActiveLayer", fSolidActiveLayer_1,  fSolidActiveTrapLayer, _rotminus90, G4ThreeVector(0, -fActiveLayer_Y/2. + fActiveTrapLayerZ/2. ,0)); 


  fSolidTrapPoronLayer = new G4Trap ("TrapPoronLayer", fTrapPoronX1/2., fTrapPoronX2/2.,fTrapPoronY1/2.,fTrapPoronY2/2.,fTrapPoronZ/2.);

  fSolidRealTrapPoronLayer = new G4UnionSolid("RealTrapPoronLayer",fSolidTrapPoronLayer,fSolidTrapPoronLayer, _rot180,G4ThreeVector(0,0,-fTrapPoronZ));

  fSolidTrapVetoLayer = new G4Trap ("TrapPoronLayer", fTrapVetoX1/2., fTrapVetoX2/2.,fTrapVetoY1/2.,fTrapVetoY2/2.,fTrapVetoZ/2.);

  fSolidRealTrapVetoLayer = new G4UnionSolid("RealTrapVetoLayer",fSolidTrapVetoLayer,fSolidTrapVetoLayer, _rot180,G4ThreeVector(0,0,-fTrapVetoZ));

  fSolidTrapCFLayer = new G4Trap ("TrapPoronLayer", fTrapCFX1/2., fTrapCFX2/2.,fTrapCFY1/2.,fTrapCFY2/2.,fTrapCFZ/2.);

  fSolidRealTrapCFLayer = new G4UnionSolid("RealTrapCFLayer",fSolidTrapCFLayer,fSolidTrapCFLayer, _rot180,G4ThreeVector(0,0,-fTrapCFZ));
 
  fSolidExternalPoronSupportA = new G4Box("ExternalPoronSupportA", fExternalPoronSupportA_X/2.,fExternalPoronSupportA_Y/2.,fExternalPoronSupportA_Z/2.);

  fSolidExternalPoronSupportB = new G4Box("ExternalPoronSupportB", fExternalPoronSupportB_X/2.,fExternalPoronSupportB_Y/2.,fExternalPoronSupportB_Z/2.);

  fSolidExternalPoronSupportC = new G4Box("ExternalPoronSupportC", fExternalPoronSupportC_X/2.,fExternalPoronSupportC_Y/2.,fExternalPoronSupportC_Z/2.);

  fSolidExternalPoronSupportD = new G4Box("ExternalPoronSupportD", fExternalPoronSupportD_X/2.,fExternalPoronSupportD_Y/2.,fExternalPoronSupportD_Z/2.);
  
  fSolidExternalPoronSupportE = new G4Box("ExternalPoronSupportE", fExternalPoronSupportE_X/2.,fExternalPoronSupportE_Y/2.,fExternalPoronSupportE_Z/2.);

  fSolidExternalPoronSupport1 = new G4UnionSolid("ExternalPoronSupport1",fSolidExternalPoronSupportE,fSolidExternalPoronSupportD,0,G4ThreeVector(fExternalPoronSupport1_X,fExternalPoronSupport1_Y,0));

  fSolidExternalPoronSupport2 = new G4UnionSolid("ExternalPoronSupport2",fSolidExternalPoronSupport1,fSolidExternalPoronSupportC,0,G4ThreeVector(fExternalPoronSupport2_X,0,0));

  fSolidExternalPoronSupport3 = new G4UnionSolid("ExternalPoronSupport3",fSolidExternalPoronSupport2,fSolidExternalPoronSupportB,0,G4ThreeVector(fExternalPoronSupport3_X,fExternalPoronSupport3_Y,0));

  fSolidExternalPoronSupport =  new G4UnionSolid("ExternalPoronSupport",fSolidExternalPoronSupport3,fSolidExternalPoronSupportA,0,G4ThreeVector(fExternalPoronSupport4_X,fExternalPoronSupport4_Y,fExternalPoronSupport4_Z));

  fSolidTeflonLYSO = new G4Box("TeflonLYSO",fTeflonLYSO_X/2.,fTeflonLYSO_Y/2.,fTeflonLYSO_Z/2.); 
  
  fSolidCFBlockContainer = new G4SubtractionSolid("CFBlockContainer",fSolidCFBlockContainerExt,fSolidCFBlockContainerInt);

  fSolidSingleCrystalBlockContainer = new G4Box("SingleCrystalBlockContainer",fCFBlockContainerExt_X/2.,fCFBlockContainerExt_Y/2.,fCFBlockContainerExt_Z/2.);

  fSolidCrystalActiveBlock = new G4Box("CrystalActiveBlock",fCrystalBlock_X/2.,fCrystalBlock_Y/2.,fCrystalBlock_Z/2.);
  
  fSolidTeflonContainerExt = new G4Box("TeflonContainerExt",TeflonContainerExt_X/2.,TeflonContainerExt_Y/2.,TeflonContainerExt_Z/2.);
  
  fSolidTeflonContainer = new G4SubtractionSolid("TeflonContainer",fSolidTeflonContainerExt,fSolidCrystalActiveBlock);
  
  fSolidCFCrystalPanelDown = new G4Box("CFCrystalPanelDown",fCFCrystalPanelDown_X/2.,fCFCrystalPanelDown_Y/2.,fCFCrystalPanelDown_Z/2);
  
  fSolidCylinderPanelDown = new G4Cons("CylinderPanelDown", fCylinderRadiusMin1,  fCylinderRadiusMax1, fCylinderRadiusMin2, fCylinderRadiusMax2, fCylinderZ/2., 8 ,8);
  
  fSolidCFCrystalPanelDownTemp = new G4Box("CFCrystalPanelDown2",fCFCrystalPanelDown2_X/2., fCFCrystalPanelDown2_Y/2.,fCFCrystalPanelDown2_Z/2.);

  fSolidCFCrystalPanelDownHole = new G4Box("CFCrystalHole",fCFCrystalHole_X/2.,fCFCrystalHole_Y/2.,fCFCrystalHole_Z/2.);
  
  fSolidPoronPlateO = new G4Box("PoronPlateO", fPoronPlateO_X/2., fPoronPlateO_Y/2.,  fPoronPlateO_Z/2. );

  fSolidPoronPlateV = new G4Box("PoronPlateV", fPoronPlateV_X/2. ,fPoronPlateV_Y/2. ,fPoronPlateV_Z/2. );
  
  
  G4double nholes = 3;
  G4double posX[3];
  G4double posY[3];
  G4double holeDist = 37.7*mm;
  G4double holeDist2 = 50.8*mm;

  posX[0] = -fCFCrystalPanelDown_X/2.+ holeDist ;
  posY[0] = -fCFCrystalPanelDown_Y/2.+ holeDist ;
  for (G4int i = 1; i < nholes; i++) {
    posX[i] = posX[0]+ i*holeDist2;
    posY[i] = posY[0]+ i*holeDist2;
  }
  
  fSolidCylinderPanelDown1 = new G4SubtractionSolid("Hole0", fSolidCFCrystalPanelDown,  fSolidCylinderPanelDown, 0, G4ThreeVector(posX[0], posY[0], 0));
  fSolidCylinderPanelDown2 = new G4SubtractionSolid("Hole1", fSolidCylinderPanelDown1,  fSolidCylinderPanelDown, 0, G4ThreeVector(posX[0], posY[1], 0));
  fSolidCylinderPanelDown3 = new G4SubtractionSolid("Hole2", fSolidCylinderPanelDown2,  fSolidCylinderPanelDown, 0, G4ThreeVector(posX[0], posY[2], 0));
  fSolidCylinderPanelDown4 = new G4SubtractionSolid("Hole3", fSolidCylinderPanelDown3,  fSolidCylinderPanelDown, 0, G4ThreeVector(posX[1], posY[0], 0));
  fSolidCylinderPanelDown5 = new G4SubtractionSolid("Hole4", fSolidCylinderPanelDown4,  fSolidCylinderPanelDown, 0, G4ThreeVector(posX[1], posY[1], 0));
  fSolidCylinderPanelDown6 = new G4SubtractionSolid("Hole5", fSolidCylinderPanelDown5,  fSolidCylinderPanelDown, 0, G4ThreeVector(posX[1], posY[2], 0));
  fSolidCylinderPanelDown7 = new G4SubtractionSolid("Hole6", fSolidCylinderPanelDown6,  fSolidCylinderPanelDown, 0, G4ThreeVector(posX[2], posY[0], 0));
  fSolidCylinderPanelDown8 = new G4SubtractionSolid("Hole7", fSolidCylinderPanelDown7,  fSolidCylinderPanelDown, 0, G4ThreeVector(posX[2], posY[1], 0));
  fSolidPanelDown          = new G4SubtractionSolid("Hole8", fSolidCylinderPanelDown8,  fSolidCylinderPanelDown, 0, G4ThreeVector(posX[2], posY[2], 0));
  
  
  G4double posSquareX[3];
  G4double posSquareY[3];
  G4double squareDist = 6*mm;

  posSquareX[0] = -fCFCrystalPanelDown2_X/2. + squareDist + fCFCrystalHole_X/2;
  posSquareY[0] = -fCFCrystalPanelDown2_Y/2. + squareDist + fCFCrystalHole_Y/2;

  for (G4int j = 1; j < nholes; j++) {
    posSquareX[j] = posSquareX[0]+ j*(squareDist+ fCFCrystalHole_X/2);
    posSquareY[j] = posSquareY[0]+ j*(squareDist+ fCFCrystalHole_Y/2 );
  }
  
  fSolidSquarePanelDown1 = new G4SubtractionSolid("SquadreHole0", fSolidCFCrystalPanelDownTemp, fSolidCFCrystalPanelDownHole, 0, G4ThreeVector(posX[0], posY[0], 0));
  fSolidSquarePanelDown2 = new G4SubtractionSolid("SquadreHole1", fSolidSquarePanelDown1, fSolidCFCrystalPanelDownHole, 0, G4ThreeVector(posX[0], posY[1], 0));
  fSolidSquarePanelDown3 = new G4SubtractionSolid("SquadreHole2", fSolidSquarePanelDown2, fSolidCFCrystalPanelDownHole, 0, G4ThreeVector(posX[0], posY[2], 0));
  fSolidSquarePanelDown4 = new G4SubtractionSolid("SquadreHole3", fSolidSquarePanelDown3, fSolidCFCrystalPanelDownHole, 0, G4ThreeVector(posX[1], posY[0], 0));
  fSolidSquarePanelDown5 = new G4SubtractionSolid("SquadreHole4", fSolidSquarePanelDown4, fSolidCFCrystalPanelDownHole, 0, G4ThreeVector(posX[1], posY[1], 0));
  fSolidSquarePanelDown6 = new G4SubtractionSolid("SquadreHole5", fSolidSquarePanelDown5, fSolidCFCrystalPanelDownHole, 0, G4ThreeVector(posX[1], posY[2], 0));
  fSolidSquarePanelDown7 = new G4SubtractionSolid("SquadreHole6", fSolidSquarePanelDown6, fSolidCFCrystalPanelDownHole, 0, G4ThreeVector(posX[2], posY[0], 0));
  fSolidSquarePanelDown8 = new G4SubtractionSolid("SquadreHole7", fSolidSquarePanelDown7, fSolidCFCrystalPanelDownHole, 0, G4ThreeVector(posX[2], posY[1], 0));
  fSolidPanelDown2          = new G4SubtractionSolid("SquadreHole8", fSolidSquarePanelDown8, fSolidCFCrystalPanelDownHole, 0, G4ThreeVector(posX[2], posY[2], 0));
  

  
  fSolidCrystalBlockContainer = new G4Box("CrystalBlockContainer",fCrystalBlockContainer_X/2.,fCrystalBlockContainer_Y/2.,fCrystalBlockContainer_Z/2.);
  
  fSolidCrystalBlockRawContainer = new G4Box("CrystalBlockRawContainer",fCrystalBlockRawContainer_X/2.,fCrystalBlockRawContainer_Y/2.,fCrystalBlockRawContainer_Z/2.);
  fSolidCFCrystalSideX = new G4Box("CFCrystalSideX",fCFCrystalSideX_X/2.,fCFCrystalSideX_Y/2.,fCFCrystalSideX_Z/2.);
  
  fSolidCFCrystalSideYBig = new G4Box("CFCrystalSideYBig",fCFCrystalSideYBig_X/2.,fCFCrystalSideYBig_Y/2.,fCFCrystalSideYBig_Z/2.);
  
  fSolidCFFront = new  G4Box("CFSupport",fCFFront_X/2.,fCFFront_Y/2.,fCFFront_Z/2.); 
  
  fSolidCFLat = new  G4Box("CFSupport",fCFLat_X/2.,fCFLat_Y/2.,fCFLat_Z/2.); 

  fSolidCFFrontPO = new  G4Box("CFSupportPO",fCFFrontPO_X/2.,fCFFrontPO_Y/2.,fCFFrontPO_Z/2.); 

  fSolidCFLatPO = new  G4Box("CFSupportPO",fCFLatPO_X/2.,fCFLatPO_Y/2.,fCFLatPO_Z/2.); 

  fSolidCFSuppO = new G4Box("CFSupport",fCFSuppO_X/2.,fCFSuppO_Y/2.,fCFSuppO_Z/2.);   

  fSolidCFSuppV = new G4Box("CFSupport",fCFSuppV_X/2.,fCFSuppV_Y/2.,fCFSuppV_Z/2.);
 
  fSolidCFSuppOA = new G4Box("CFSupportTop",fCFSuppAO_X/2.,fCFSuppAO_Y/2.,fCFSuppAO_Z/2.);   

  fSolidCFSuppVA = new G4Box("CFSupportTop",fCFSuppAV_X/2.,fCFSuppAV_Y/2.,fCFSuppAV_Z/2.);
 
  fSolidCFSuppLat = new G4Box("CFSupport",fCFSuppLat_X/2.,fCFSuppLat_Y/2.,fCFSuppLat_Z/2.);
  
  fSolidCFSuppLatA = new G4Box("CFSupportTop",fCFSuppLatA_X/2.,fCFSuppLatA_Y/2.,fCFSuppLatA_Z/2.);
  
  fSolidCFSuppFront = new G4Box("CFSupport",fCFSuppFront_X/2.,fCFSuppFront_Y/2.,fCFSuppFront_Z/2.);

  fSolidCFSuppFrontA = new G4Box("CFSupportTop",fCFSuppFrontA_X/2.,fCFSuppFrontA_Y/2.,fCFSuppFrontA_Z/2.);
  
  fSolidCFSuppPorO = new G4Box("PORSupport",fCFSuppPorO_X/2.,fCFSuppPorO_Y/2.,fCFSuppPorO_Z/2.);   
  
  fSolidCFSuppPorV = new G4Box("PORSupport",fCFSuppPorV_X/2.,fCFSuppPorV_Y/2.,fCFSuppPorV_Z/2.);
  
  fSolidCFSuppPoron = new G4UnionSolid("PORSupport",fSolidCFSuppPorO, fSolidCFSuppPorV);    // creazione croce poron

  fSolidCFSuppStepOV = new G4UnionSolid("CFSupport",fSolidCFSuppO,fSolidCFSuppV);    // creazione croce carbonio
  
  fSolidCFSuppStepOVA = new G4UnionSolid("CFSupportTop",fSolidCFSuppOA,fSolidCFSuppVA);    // creazione croce carbonio
  
  fSolidCFSuppStepLat1 = new G4UnionSolid("CFSupport",fSolidCFSuppStepOV,fSolidCFSuppLat,0,    // aggiunta blocchi ai lati della croce
					  G4ThreeVector(fSolidCFSuppStepLat_X,fSolidCFSuppStepLat_Y,fSolidCFSuppStepLat_Z));
  
  fSolidCFSuppStepLat2 = new G4UnionSolid("CFSupport",fSolidCFSuppStepLat1,fSolidCFSuppLat,0,
					  G4ThreeVector(-fSolidCFSuppStepLat_X,fSolidCFSuppStepLat_Y,fSolidCFSuppStepLat_Z));
  
  fSolidCFSuppStepFront1 = new G4UnionSolid("CFSupport",fSolidCFSuppStepLat2,fSolidCFSuppFront,0,
					    G4ThreeVector(fSolidCFSuppStepFront_X,fSolidCFSuppStepFront_Y,fSolidCFSuppStepFront_Z));
  
  fSolidCFSupp = new G4UnionSolid("CFSupport",fSolidCFSuppStepFront1,fSolidCFSuppFront,0,
				  G4ThreeVector(fSolidCFSuppStepFront_X,-fSolidCFSuppStepFront_Y,fSolidCFSuppStepFront_Z));

  fSolidCFSuppStepLat1A = new G4UnionSolid("CFSupportTop",fSolidCFSuppStepOVA,fSolidCFSuppLatA,0,
					   G4ThreeVector(fSolidCFSuppStepLat_X,fSolidCFSuppStepLat_Y,fSolidCFSuppStepLat_Z));
  
  fSolidCFSuppStepLat2A = new G4UnionSolid("CFSupportTop",fSolidCFSuppStepLat1A,fSolidCFSuppLatA,0,
					   G4ThreeVector(-fSolidCFSuppStepLat_X,fSolidCFSuppStepLat_Y,fSolidCFSuppStepLat_Z));
  
  fSolidCFSuppStepFront1A = new G4UnionSolid("CFSupportTop",fSolidCFSuppStepLat2A,fSolidCFSuppFrontA,0,
					     G4ThreeVector(fSolidCFSuppStepFront_X,fSolidCFSuppStepFront_Y,fSolidCFSuppStepFront_Z));
  
  fSolidCFSuppA = new G4UnionSolid("CFSupporTop",fSolidCFSuppStepFront1A,fSolidCFSuppFrontA,0,
				   G4ThreeVector(fSolidCFSuppStepFront_X,-fSolidCFSuppStepFront_Y,fSolidCFSuppStepFront_Z));

  fSolidPoronLat = new G4Box("Poron",fPoronLat_X/2.,fPoronLat_Y/2.,fPoronLat_Z/2.);

  fSolidPoronFront = new G4Box("Poron",fPoronFront_X/2.,fPoronFront_Y/2.,fPoronFront_Z/2.);

  fSolidPoronLatPO = new G4Box("PoronPO",fPoronLatPO_X/2.,fPoronLatPO_Y/2.,fPoronLatPO_Z/2.);

  fSolidPoronFrontPO = new G4Box("PoronPO",fPoronFrontPO_X/2.,fPoronFrontPO_Y/2.,fPoronFrontPO_Z/2.);
 





  // VETO STRUCTURES
  fSolidPoronLatX = new G4Box("Poron",fPoronLatX_X/2.,fPoronLatX_Y/2.,fPoronLatX_Z/2.);

  fSolidPoronLatXRight = new G4Box("Poron",fPoronLatXRight_X/2.,fPoronLatXRight_Y/2.,fPoronLatXRight_Z/2.);

  fSolidPoronLatXHole = new G4Box("Poron",fPoronLatXHole_X/2.,fPoronLatXHole_Y/2.,fPoronLatXHole_Z/2.);

  fSolidPoronLatXHoleLeft = new G4Box("Poron",fPoronLatXHoleLeft_X/2.,fPoronLatXHoleLeft_Y/2.,fPoronLatXHoleLeft_Z/2.);

  fSolidPoronLatX_1 = new G4UnionSolid("SolidPoronLatX_1",fSolidPoronLatX ,fSolidPoronLatXRight,0,G4ThreeVector(fPoronLat1X_X,0,fPoronLat1X_Z));
  
  fSolidPoronLatX_2 = new G4UnionSolid("SolidPoronLatX_2",fSolidPoronLatX_1 ,fSolidPoronLatXHole,0,G4ThreeVector(fPoronLat2X_X,0,fPoronLat2X_Z));

  fSolidPoronLatX_3 = new G4UnionSolid("SolidPoronLatX_3",fSolidPoronLatX_2 ,fSolidPoronLatXHoleLeft,0,G4ThreeVector(fPoronLat3X_X,0,0));

  fSolidVetoLatX = new G4Box("Veto",fVetoLatX_X/2.,fVetoLatX_Y/2.,fVetoLatX_Z/2.);

  fSolidVetoLatXRight = new G4Box("Veto",fVetoLatXRight_X/2.,fVetoLatXRight_Y/2.,fVetoLatXRight_Z/2.);

  fSolidVetoLatXLeft = new G4Box("Veto",fVetoLatXLeft_X/2.,fVetoLatXLeft_Y/2.,fVetoLatXLeft_Z/2.);

  fSolidVetoLatX_1 = new G4UnionSolid("SolidVetoLatX_1",fSolidVetoLatX ,fSolidVetoLatXRight,0,G4ThreeVector(fVetoLat1X_X,0,fVetoLat1X_Z));

  fSolidVetoLatX_2 = new G4UnionSolid("SolidVetoLatX_2",fSolidVetoLatX_1 ,fSolidVetoLatXLeft,0,G4ThreeVector(fVetoLat2X_X,0,fVetoLat2X_Z));

  fSolidPoronLatY = new G4Box("Poron",fPoronLatY_X/2.,fPoronLatY_Y/2.,fPoronLatY_Z/2.);

  fSolidPoronLatYUp = new G4Box("Poron",fPoronLatYUp_X/2.,fPoronLatYUp_Y/2.,fPoronLatYUp_Z/2.);

  fSolidPoronLatYHole = new G4Box("Poron",fPoronLatYHole_X/2.,fPoronLatYHole_Y/2.,fPoronLatYHole_Z/2.);

  fSolidPoronLatYHoleDown = new G4Box("Poron",fPoronLatYHoleDown_X/2.,fPoronLatYHoleDown_Y/2.,fPoronLatYHoleDown_Z/2.);
 
  fSolidPoronLatY_1 = new G4UnionSolid("SolidPoronLatY_1",fSolidPoronLatY ,fSolidPoronLatYUp,0,G4ThreeVector(0, fPoronLat1Y_Y, fPoronLat1Y_Z));
  
  fSolidPoronLatY_2 = new G4UnionSolid("SolidPoronLatY_2",fSolidPoronLatY_1 ,fSolidPoronLatYHole,0,G4ThreeVector(0, fPoronLat2Y_Y, fPoronLat2Y_Z));

  fSolidPoronLatY_3 = new G4UnionSolid("SolidPoronLatY_3",fSolidPoronLatY_2 ,fSolidPoronLatYHoleDown,0,G4ThreeVector(0,fPoronLat3Y_Y,0));

  fSolidVetoLatY = new G4Box("Veto",fVetoLatY_X/2.,fVetoLatY_Y/2.,fVetoLatY_Z/2.);

  fSolidVetoLatYUp = new G4Box("Veto",fVetoLatYUp_X/2.,fVetoLatYUp_Y/2.,fVetoLatYUp_Z/2.);

  fSolidVetoLatYDown = new G4Box("Veto",fVetoLatYDown_X/2.,fVetoLatYDown_Y/2.,fVetoLatYDown_Z/2.);

  fSolidVetoLatY_1 = new G4UnionSolid("SolidVetoLatX_1",fSolidVetoLatY ,fSolidVetoLatYUp,0,G4ThreeVector(0,fVetoLat1Y_Y, fVetoLat1Y_Z));

  fSolidVetoLatY_2 = new G4UnionSolid("SolidVetoLatX_2",fSolidVetoLatY_1 ,fSolidVetoLatYDown,0,G4ThreeVector(0,fVetoLat2Y_Y,fVetoLat2Y_Z));

  fSolidCFVetoLatX = new G4Box("CFSupport",fCFVetoLatX_X/2.,fCFVetoLatX_Y/2.,fCFVetoLatX_Z/2.);

  fSolidCFVetoLatXHole = new G4Box("CFSupport",fCFVetoLatXHole_X/2.,fCFVetoLatXHole_Y/2.,fCFVetoLatXHole_Z/2.);

  fSolidCFVetoLatXHoleRight = new G4Box("CFSupport",fCFVetoLatXHoleRight_X/2.,fCFVetoLatXHoleRight_Y/2.,fCFVetoLatXHoleRight_Z/2.);

  fSolidCFVetoLatX_1 = new G4UnionSolid("SolidCFLatX_1",fSolidCFVetoLatX , fSolidCFVetoLatXHole,0,G4ThreeVector(fCFVetoLat1X_X,0,fCFVetoLat1X_Z));

  fSolidCFVetoLatX_2 = new G4UnionSolid("SolidCFLatX_2",fSolidCFVetoLatX_1 , fSolidCFVetoLatXHoleRight,0,G4ThreeVector(fCFVetoLat2X_X,0,0));

  fSolidCFVetoLatY = new G4Box("CFSupport",fCFVetoLatY_X/2.,fCFVetoLatY_Y/2.,fCFVetoLatY_Z/2.);

  fSolidCFVetoLatYHole = new G4Box("CFSupport",fCFVetoLatYHole_X/2.,fCFVetoLatYHole_Y/2.,fCFVetoLatYHole_Z/2.);

  fSolidCFVetoLatYHoleDown = new G4Box("CFSupport",fCFVetoLatYHoleDown_X/2.,fCFVetoLatYHoleDown_Y/2.,fCFVetoLatYHoleDown_Z/2.);
   
  fSolidCFVetoLatY_1 = new G4UnionSolid("SolidCFLatY_1",fSolidCFVetoLatY , fSolidCFVetoLatYHole,0,G4ThreeVector(0, fCFVetoLat1Y_Y ,fCFVetoLat1Y_Z));

  fSolidCFVetoLatY_2 = new G4UnionSolid("SolidCFLatX_2",fSolidCFVetoLatY_1 , fSolidCFVetoLatYHoleDown,0,G4ThreeVector(0, fCFVetoLat2Y_Y, 0));
  //////////////////////////////////////






  //fSolidHoneyCombLatY = new G4Box("HoneyCombCore",fHoneyCombLatY_X/2.,fHoneyCombLatY_Y/2.,fHoneyCombLatY_Z/2.);

  //fSolidHoneyCombSkinLatY = new G4Box("HoneyCombSkin",fHoneyCombSkinLatY_X/2.,fHoneyCombSkinLatY_Y/2.,fHoneyCombSkinLatY_Z/2.);

  //fSolidHoneyCombLatX = new G4Box("HoneyCombCore",fHoneyCombLatX_X/2.,fHoneyCombLatX_Y/2.,fHoneyCombLatX_Z/2.);

  //fSolidHoneyCombSkinLatX = new G4Box("HoneyCombSkin",fHoneyCombSkinLatX_X/2.,fHoneyCombSkinLatX_Y/2.,fHoneyCombSkinLatX_Z/2.);

  //fSolidHoneyCombBottom = new G4Box("HoneyCombCore",fHoneyCombBottom_X/2.,fHoneyCombBottom_Y/2.,fHoneyCombBottom_Z/2.);

  //  fSolidHoneyCombSkinBottom = new G4Box("HoneyCombSkin",fHoneyCombSkinBottom_X/2.,fHoneyCombSkinBottom_Y/2.,fHoneyCombSkinBottom_Z/2.);
  




  // S1 Scintillator
  fLogicS1 = new G4LogicalVolume(fSolidS1,vacuum,"S1");

  fLogicS1SuppBack = new G4LogicalVolume(fSolidS1SuppBack,cfMat,"S1SupportBack");

  fLogicS1SuppHoleBar = new G4LogicalVolume(fSolidS1SuppHoleBar,cfMat,"S1SupportHoleBar");

  fLogicS1SuppBar = new G4LogicalVolume(fSolidS1SuppBar,cfMat,"S1SupportBar");

  fLogicS1SuppThinBar = new G4LogicalVolume(fSolidS1SuppThinBar,cfMat,"S1SupportThinBar");

  fLogicS1SuppThinBarBack = new G4LogicalVolume(fSolidS1SuppThinBarBack,cfMat,"S1SupportThinBarBack");

  fLogicS1SuppFront = new G4LogicalVolume(fSolidS1SuppFront,cfMat,"S1SupportFront");

  fLogicS1SuppPoronFront = new G4LogicalVolume(fSolidS1SuppPoronFront,porMat,"S1SupportPoronFront");

  fLogicS1ScintContainerP = new G4LogicalVolume(fSolidS1ScintContainer,vacuum,"S1ScintContainer1");

  fLogicS1ScintContainerM = new G4LogicalVolume(fSolidS1ScintContainer,vacuum,"S1ScintContainer2");

  fLogicS1Scint = new G4LogicalVolume(fSolidS1Scint,scintLayerMat,"S1Scint");
  /////////////////////////////


  // Calorimeter
  fLogicCaloBox = new G4LogicalVolume(fSolidCaloBox,vacuum,"Calorimeter");

  fLogicScintBox = new G4LogicalVolume(fSolidScintBox,vacuum,"CalorimeterScint");

  fLogicCrystalBox = new G4LogicalVolume(fSolidCrystalBox,vacuum,"CalorimeterCrystal");

  fLogicScintLayer = new G4LogicalVolume(fSolidLayer,vacuum,"LayerScint");

  fLogicLastScintLayer = new G4LogicalVolume(fSolidLastLayer,vacuum,"LastLayerScint");
  
  fLogicScintActiveLayer = new G4LogicalVolume(fSolidActiveLayer,scintLayerMat,"ActiveLayerScint");

  fLogicRealTrapPoronLayer = new G4LogicalVolume(fSolidRealTrapPoronLayer,porMat,"RealTrapPoronLayer");

  fLogicRealTrapVetoLayer = new G4LogicalVolume(fSolidRealTrapVetoLayer,vetoLayerMat,"RealTrapVetoLayer");

  fLogicRealTrapCFLayer = new G4LogicalVolume(fSolidRealTrapCFLayer,cfMat,"RealTrapCFLayer");
 
  fLogicExternalPoronSupport = new G4LogicalVolume(fSolidExternalPoronSupport,porMat,"ExternalPoronSupport");

  fLogicTeflonLYSO =  new G4LogicalVolume(fSolidTeflonLYSO,vacuum,"TeflonLYSO");
 
  fLogicCrystalBlockContainer = new G4LogicalVolume(fSolidSingleCrystalBlockContainer,vacuum,"SingleCrystalBlockContainer");
  
  fLogicCrystalBlockPlaneContainer = new G4LogicalVolume(fSolidCrystalBlockContainer,vacuum,"CrystalBlockPlane");
   
  fLogicCrystalBlockRawContainer = new G4LogicalVolume(fSolidCrystalBlockRawContainer,vacuum,"CrystalBlockRaw");

  fLogicCrystalActiveBlock = new G4LogicalVolume(fSolidCrystalActiveBlock,crystalLayerMat,"CrystalActiveBlock");

  fLogicCFCrystalPanelDown = new G4LogicalVolume(fSolidPanelDown, cfMat, "CFSupport");

  fLogicCFCrystalPanelDown2 = new G4LogicalVolume(fSolidPanelDown2, cfMat, "CFSupport");
  
  fLogicCFBlockContainer = new G4LogicalVolume(fSolidCFBlockContainer,cfMat,"CFSupport");

  fLogicTeflonContainer = new G4LogicalVolume(fSolidTeflonContainer,teflon,"TeflonLayer");///////////////////////////////////////////

  fLogicCFCrystalSideX = new G4LogicalVolume(fSolidCFCrystalSideX,cfMat,"CFSupport");

  fLogicCFCrystalSideYBig = new G4LogicalVolume(fSolidCFCrystalSideYBig,cfMat,"CFSupport");

  fLogicCFFront = new G4LogicalVolume(fSolidCFFront,cfMat,"CFSupport");

  fLogicCFLat = new G4LogicalVolume(fSolidCFLat,cfMat,"CFSupport");

  fLogicCFFrontPO = new G4LogicalVolume(fSolidCFFrontPO,cfMat,"CFSupportPO");

  fLogicCFLatPO = new G4LogicalVolume(fSolidCFLatPO,cfMat,"CFSupportPO");
  
  fLogicCFSuppPoron = new G4LogicalVolume(fSolidCFSuppPoron, porMat,"PORSupport");
  
  fLogicCFSupp = new G4LogicalVolume(fSolidCFSupp,cfMat,"CFSupport");

  fLogicCFSuppTop = new G4LogicalVolume(fSolidCFSuppA,cfMat,"CFSupportTop"); 

  fLogicPoronLat = new G4LogicalVolume(fSolidPoronLat,porMat,"Poron");

  fLogicPoronFront = new G4LogicalVolume(fSolidPoronFront,porMat,"Poron"); 

  fLogicPoronLatPO = new G4LogicalVolume(fSolidPoronLatPO,porMat,"PoronPO");

  fLogicPoronFrontPO = new G4LogicalVolume(fSolidPoronFrontPO,porMat,"PoronPO");

 
  // vecchi VETO
  //   fLogicPoronLatX = new G4LogicalVolume(fSolidPoronLatX,porMat,"Poron");

  //   fLogicPoronLatXRight = new G4LogicalVolume(fSolidPoronLatXRight,porMat,"Poron");

  //   fLogicPoronLatXHole = new G4LogicalVolume(fSolidPoronLatXHole,porMat,"Poron");

  //   fLogicPoronLatXHoleLeft = new G4LogicalVolume(fSolidPoronLatXHoleLeft,porMat,"Poron");

 



  // GOOD VETO 
  fLogicPoronLatX = new G4LogicalVolume(fSolidPoronLatX_3,porMat,"Poron");
 
  fLogicPoronLatY = new G4LogicalVolume(fSolidPoronLatY_3,porMat,"Poron");

  fLogicVetoLatX = new G4LogicalVolume(fSolidVetoLatX_2,vetoLayerMat,"Veto");
  
  fLogicVetoLatX2 = new G4LogicalVolume(fSolidVetoLatX_2,vetoLayerMat,"Veto");

  fLogicVetoLatY = new G4LogicalVolume(fSolidVetoLatY_2,vetoLayerMat,"Veto");

  fLogicVetoLatY2 = new G4LogicalVolume(fSolidVetoLatY_2,vetoLayerMat,"Veto");

  fLogicCFVetoLatX = new G4LogicalVolume(fSolidCFVetoLatX_2,cfMat,"CFSupport");

  fLogicCFVetoLatY = new G4LogicalVolume(fSolidCFVetoLatY_2, cfMat,"CFSupport");
  ///////////////////



  // vecchi veto
  //fLogicVetoLatX = new G4LogicalVolume(fSolidVetoLatX,vetoLayerMat,"Veto");

  //fLogicVetoLatXRight = new G4LogicalVolume(fSolidVetoLatXRight,vetoLayerMat,"Veto");

  //fLogicVetoLatXLeft = new G4LogicalVolume(fSolidVetoLatXLeft,vetoLayerMat,"Veto");

  //fLogicPoronLatY = new G4LogicalVolume(fSolidPoronLatY,porMat,"Poron");

  //fLogicPoronLatYUp = new G4LogicalVolume(fSolidPoronLatYUp,porMat,"Poron");

  //fLogicPoronLatYHole = new G4LogicalVolume(fSolidPoronLatYHole,porMat,"Poron");

  //fLogicPoronLatYHoleDown = new G4LogicalVolume(fSolidPoronLatYHoleDown,porMat,"Poron");

  //fLogicVetoLatY = new G4LogicalVolume(fSolidVetoLatY,vetoLayerMat,"Veto");

  //fLogicVetoLatYUp = new G4LogicalVolume(fSolidVetoLatYUp,vetoLayerMat,"Veto");

  //fLogicVetoLatYDown = new G4LogicalVolume(fSolidVetoLatYDown,vetoLayerMat,"Veto");

  //fLogicCFVetoLatX = new G4LogicalVolume(fSolidCFVetoLatX,cfMat,"CFSupport");

  //fLogicCFVetoLatXHole = new G4LogicalVolume(fSolidCFVetoLatXHole,cfMat,"CFSupport");

  //fLogicCFVetoLatXHoleRight = new G4LogicalVolume(fSolidCFVetoLatXHoleRight,cfMat,"CFSupport");

  //  fLogicCFVetoLatY = new G4LogicalVolume(fSolidCFVetoLatY,cfMat,"CFSupport");

  //  fLogicCFVetoLatYHole = new G4LogicalVolume(fSolidCFVetoLatYHole,cfMat,"CFSupport");

  //fLogicCFVetoLatYHoleDown = new G4LogicalVolume(fSolidCFVetoLatYHoleDown,cfMat,"CFSupport");

  fLogicPoronPlateO = new G4LogicalVolume(fSolidPoronPlateO, porMat, "PoronPlateO");

  fLogicPoronPlateV = new G4LogicalVolume(fSolidPoronPlateV, porMat, "PorPlateV");
 
  //fLogicHoneyCombLatY = new G4LogicalVolume(fSolidHoneyCombLatY,hcMat,"HoneyCombCore"); 

  //fLogicHoneyCombSkinLatY = new G4LogicalVolume(fSolidHoneyCombSkinLatY,cfMat,"HoneyCombSkin"); 

  //fLogicHoneyCombLatX = new G4LogicalVolume(fSolidHoneyCombLatX,hcMat,"HoneyCombCore"); 

  //fLogicHoneyCombSkinLatX = new G4LogicalVolume(fSolidHoneyCombSkinLatX,cfMat,"HoneyCombSkin"); 

  //  fLogicHoneyCombBottom = new G4LogicalVolume(fSolidHoneyCombBottom,hcMat,"HoneyCombCore"); 

  //fLogicHoneyCombSkinBottom = new G4LogicalVolume(fSolidHoneyCombSkinBottom,cfMat,"HoneyCombSkin"); 


  fLogicScintActiveLayer->SetSensitiveDetector(caloSD);
  fLogicCrystalActiveBlock->SetSensitiveDetector(caloSD);
  fLogicS1Scint->SetSensitiveDetector(caloSD);

  fLogicRealTrapVetoLayer->SetSensitiveDetector(vetoSD);
  fLogicVetoLatX->SetSensitiveDetector(vetoSD);
  fLogicVetoLatX2->SetSensitiveDetector(vetoSD);
  fLogicVetoLatY->SetSensitiveDetector(vetoSD);
  fLogicVetoLatY2->SetSensitiveDetector(vetoSD);
  // bisogna separare i due piani



  

  fPhysiS1 = new G4PVPlacement(0,
			       G4ThreeVector(fPhysiS1_X,fPhysiS1_Y,fPhysiS1_Z + ShiftOrigin ),
			       "S1",
			       fLogicS1,                
			       motherVolume,                
			       false,                       
			       0,true);       

  
  fPhysiS1SuppBack = new G4PVPlacement(0,
				       G4ThreeVector(fPhysiS1SuppBack_X, fPhysiS1SuppBack_Y, fPhysiS1SuppBack_Z),
				       "S1SupportBack",
				       fLogicS1SuppBack,                
				       fPhysiS1,                
				       false,                       
				       0,true); 
  
  fPhysiS1SuppHoleBar = new G4PVPlacement(0,
					  G4ThreeVector( fPhysiS1SuppHoleBar_X , fPhysiS1SuppHoleBar_Y, fPhysiS1SuppHoleBar_Z),
					  "S1SupportHoleBar",
					  fLogicS1SuppHoleBar,                
					  fPhysiS1,                
					  false,                       
					  0,true);   
  
  fPhysiS1SuppBar1 = new G4PVPlacement(0,
				       G4ThreeVector( fPhysiS1SuppBar_X , fPhysiS1SuppBar_Y, fPhysiS1SuppBar_Z),
				       "S1SupportBar1",
				       fLogicS1SuppBar,                
				       fPhysiS1,                
				       false,                       
				       0,true);  
 
  fPhysiS1SuppBar2 = new G4PVPlacement(0,
				       G4ThreeVector( fPhysiS1SuppBar_X , -fPhysiS1SuppBar_Y, fPhysiS1SuppBar_Z),
				       "S1SupportBar2",
				       fLogicS1SuppBar,                
				       fPhysiS1,                
				       false,                       
				       0,true);     

  fPhysiS1SuppThinBar1 = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiS1SuppThinBar_X, fPhysiS1SuppThinBar_Y, fPhysiS1SuppThinBar_Z),
					   "S1SupportThinBar1",
					   fLogicS1SuppThinBar,                
					   fPhysiS1,                
					   false,                       
					   0,true);     

  fPhysiS1SuppThinBarBack1 = new G4PVPlacement(0,
					       G4ThreeVector(fPhysiS1SuppThinBarBack_X, fPhysiS1SuppThinBarBack_Y, fPhysiS1SuppThinBarBack_Z),
					       "S1SupportThinBarBack1",
					       fLogicS1SuppThinBarBack,                
					       fPhysiS1,                
					       false,                       
					       0,true);     
  
  
  fPhysiS1SuppFrontM = new G4PVPlacement(0,
					 G4ThreeVector(fPhysiS1SuppFrontM_X,fPhysiS1SuppFrontM_Y,fPhysiS1SuppFrontM_Z),
					 "S1SupportM",
					 fLogicS1SuppFront,                
					 fPhysiS1,                
					 false,                       
					 0,true);       
  
  fPhysiS1SuppFrontPoronM = new G4PVPlacement(0,
					      G4ThreeVector(fPhysiS1SuppPoronFrontM_X,fPhysiS1SuppPoronFrontM_Y,fPhysiS1SuppPoronFrontM_Z),
					      "S1SupportPoronM",
					      fLogicS1SuppPoronFront,                
					      fPhysiS1,                
					      false,                       
					      0,true);       
 

  
  fPhysiS1SuppThinBar2 = new G4PVPlacement(0,
					   G4ThreeVector(-fPhysiS1SuppThinBar_X, fPhysiS1SuppThinBar_Y, fPhysiS1SuppThinBar_Z),
					   "S1SupportThinBar2",
					   fLogicS1SuppThinBar,                
					   fPhysiS1,                
					   false,                       
					   0,true);     
  
  fPhysiS1SuppThinBarBack2 = new G4PVPlacement(0,
					       G4ThreeVector(-fPhysiS1SuppThinBarBack_X, fPhysiS1SuppThinBarBack_Y, fPhysiS1SuppThinBarBack_Z),
					       "S1SupportThinBarBack2",
					       fLogicS1SuppThinBarBack,                
					       fPhysiS1,                
					       false,                       
					       0,true);     
 
  fPhysiS1ScintContainerP = new G4PVPlacement(0,
					      G4ThreeVector(fPhysiS1ScintContainer_X,fPhysiS1ScintContainer_Y,fPhysiS1ScintContainer_Z),
					      "S1ContainerP",
					      fLogicS1ScintContainerP,
					      fPhysiS1,
					      false,0,true);
					      
  fPhysiS1ScintContainerM = new G4PVPlacement(0,
					      G4ThreeVector(fPhysiS1ScintContainer_X,-fPhysiS1ScintContainer_Y,fPhysiS1ScintContainer_Z),
					      "S1ContainerM",
					      fLogicS1ScintContainerM,
					      fPhysiS1,
					      false,0,true);
					      
  fPhysiS1ScintP = new G4PVReplica("S1ScintillatorP",
				   fLogicS1Scint,
				   fPhysiS1ScintContainerP,
				   kYAxis,
				   fS1ScintNumber/2.,
				   fS1Scint_Y);
  
  fPhysiS1ScintM = new G4PVReplica("S1ScintillatorM",
				   fLogicS1Scint,
				   fPhysiS1ScintContainerM,
				   kYAxis,
				   fS1ScintNumber/2.,
				   fS1Scint_Y);
  ///////////////////////////////////////////////////////


  // CALO 
  
  fPhysiCaloBox = new G4PVPlacement(0,                        
				    G4ThreeVector(0,0,ShiftOrigin),   
				    "Calorimeter",               
				    fLogicCaloBox,                
				    motherVolume,                
				    false,                       
				    0,true);       
  
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
				     fNbOfReplicatedScintLayers,        
				     fCaloLayer_Z);

  fPhysiLastScintLayer = new G4PVPlacement(0,                        
					   G4ThreeVector(fPhysiLastScintLayer_X,fPhysiLastScintLayer_Y,fPhysiLastScintLayer_Z),
					   "LastScintLayer",               
					   fLogicLastScintLayer,                
					   fPhysiCaloBox,                
					   false,                       
					   0,true); 
  

  
  // piani scintillatore
  fPhysiScintActiveLayer = new G4PVPlacement(0,                   
					     G4ThreeVector(fPhysiActiveLayer_X,fPhysiActiveLayer_Y,fPhysiActiveLayer_Z),     
					     "ActiveLayerScint",              
					     fLogicScintActiveLayer,       
					     fPhysiScintLayer,        
					     false,               
					     0,true); 

 

  /*######################  Costole Poron angoli calorimetro ################################*/
  fPhysiExternalPoronSupport1 = new G4PVPlacement(0,                   
						  G4ThreeVector(fPhysiExternalPoronSupport_X,fPhysiExternalPoronSupport_Y,fPhysiExternalPoronSupport_Z + ShiftOrigin),     
						  "ExternalPoronBar1",
						  fLogicExternalPoronSupport,
						  motherVolume,        
						  false,               
						  0,true); 


  fPhysiExternalPoronSupport2 = new G4PVPlacement(_rot180Y,                   
						  G4ThreeVector(-fPhysiExternalPoronSupport_X,fPhysiExternalPoronSupport_Y,fPhysiExternalPoronSupport_Z + ShiftOrigin),     
						  "ExternalPoronBar2",
						  fLogicExternalPoronSupport,
						  motherVolume,        
						  false,               
						  0,true); 


  fPhysiExternalPoronSupport3 = new G4PVPlacement(_rot180ZY,                   
						  G4ThreeVector(fPhysiExternalPoronSupport_X,-fPhysiExternalPoronSupport_Y,fPhysiExternalPoronSupport_Z + ShiftOrigin),     
						  "ExternalPoronBar3",
						  fLogicExternalPoronSupport,
						  motherVolume,        
						  false,               
						  0,true); 


  fPhysiExternalPoronSupport4 = new G4PVPlacement(_rot180Z,                   
						  G4ThreeVector(-fPhysiExternalPoronSupport_X,-fPhysiExternalPoronSupport_Y,fPhysiExternalPoronSupport_Z + ShiftOrigin),     
						  "ExternalPoronBar4",
						  fLogicExternalPoronSupport,
						  motherVolume,        
						  false,               
						  0,true); 


  //////////////////////////////////////


  // croce di poron
  fPhysiScintCFSupp = new G4PVPlacement(0,    
					G4ThreeVector(fPhysiCFSuppPor_X, fPhysiCFSuppPor_Y, fPhysiCFSuppPor_Z ),     
					"PORSupport_1",              
					fLogicCFSuppPoron,       
					fPhysiScintLayer,        
					false,               
					0,true); 
  
  // croce di carbonio
  fPhysiScintCFSupp = new G4PVPlacement(0,    
					G4ThreeVector(fPhysiCFSupp_X,fPhysiCFSupp_Y,fPhysiCFSupp_Z),     
					"CFSupport",              
					fLogicCFSupp,       
					fPhysiScintLayer,        
					false,               
					0,true); 

  // croce di poron 2
  fPhysiScintCFSupp = new G4PVPlacement(0,
					G4ThreeVector(fPhysiCFSuppPor2_X, fPhysiCFSuppPor2_Y, fPhysiCFSuppPor2_Z),     
					"PORSupport_2",              
					fLogicCFSuppPoron,       
					fPhysiScintLayer,        
					false,               
					0,true); 
  
 
  
  // blocco carbonio sopra piano scint
  fPhysiScintCFFrontP = new G4PVPlacement(0,                   
					  G4ThreeVector(fPhysiCFFront_X,fPhysiCFFront_Y,fPhysiCFFront_Z),     
					  "CFSupportFP",              
					  fLogicCFFront,       
					  fPhysiScintLayer,        
					  false,               
					  0,true); 
  // blocco carbonio sotto piano scint
  fPhysiScintCFFrontM = new G4PVPlacement(0,                   
					  G4ThreeVector(fPhysiCFFront_X,-fPhysiCFFront_Y,fPhysiCFFront_Z),     
					  "CFSupportFM",              
					  fLogicCFFront,       
					  fPhysiScintLayer,        
					  false,               
					  0,true); 
  // blocco carbonio destra piano scint
  fPhysiScintCFLatP = new G4PVPlacement(0,                   
					G4ThreeVector(fPhysiCFLat_X,fPhysiCFLat_Y,fPhysiCFLat_Z),     
					"CFSupportLP",              
					fLogicCFLat,       
					fPhysiScintLayer,        
					false,               
					0,true); 
  // blocco carbonio sinistra piano scint
  fPhysiScintCFLatM = new G4PVPlacement(0,                   
					G4ThreeVector(-fPhysiCFLat_X,fPhysiCFLat_Y,fPhysiCFLat_Z),     
					"CFSupportLM",              
					fLogicCFLat,       
					fPhysiScintLayer,        
					false,               
					0,true); 
  
  // poron parete laterale, sinistra piano scint
  fPhysiScintPoronLatP = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiPoronLat_X,fPhysiPoronLat_Y,fPhysiPoronLat_Z),     
					   "PoronLP",              
					   fLogicPoronLat,       
					   fPhysiScintLayer,        
					   false,               
					   0,true); 

  // poron parete laterla, destra piano scint					  
  fPhysiScintPoronLatM = new G4PVPlacement(0,                   
					   G4ThreeVector(-fPhysiPoronLat_X,fPhysiPoronLat_Y,fPhysiPoronLat_Z),     
					   "PoronLM",              
					   fLogicPoronLat,       
					   fPhysiScintLayer,        
					   false,               
					   0,true); 
  
  // poron parete sopra piano scint
  fPhysiScintPoronLatUp = new G4PVPlacement(0,
					    G4ThreeVector(fPhysiPoronFront_X,-fPhysiPoronFront_Y,fPhysiPoronFront_Z),     
					    "PoronFP",              
					    fLogicPoronFront,     
					    fPhysiScintLayer,        
					    false,               
					    0,true);
 
  //poron parete sotto piano scint
  fPhysiScintPoronLatDown = new G4PVPlacement(0,
					      G4ThreeVector(fPhysiPoronFront_X,fPhysiPoronFront_Y,fPhysiPoronFront_Z),     
					      "PoronFM",              
					      fLogicPoronFront,     
					      fPhysiScintLayer,        
					      false,               
					      0,true);

  // blocco carbonio sopra croce poron
  fPhysiScintCFFrontPPO = new G4PVPlacement(0,                   
					    G4ThreeVector(fPhysiCFFrontPO_X,fPhysiCFFrontPO_Y,fPhysiCFFrontPO_Z),     
					    "CFSupportFPsmall",              
					    fLogicCFFrontPO,       
					    fPhysiScintLayer,        
					    false,               
					    0,true); 
  // blocco carbonio sotto croce  poron
  fPhysiScintCFFrontMPO = new G4PVPlacement(0,                   
					    G4ThreeVector(fPhysiCFFrontPO_X,-fPhysiCFFrontPO_Y,fPhysiCFFrontPO_Z),     
					    "CFSupportFMsmall",              
					    fLogicCFFrontPO,       
					    fPhysiScintLayer,        
					    false,               
					    0,true); 
  // blocco carbonio destra croce poron
  fPhysiScintCFLatPPO = new G4PVPlacement(0,                   
					  G4ThreeVector(fPhysiCFLatPO_X,fPhysiCFLatPO_Y,fPhysiCFLatPO_Z),     
					  "CFSupportLPsmall",              
					  fLogicCFLatPO,       
					  fPhysiScintLayer,        
					  false,               
					  0,true); 
  // blocco carbonio sinistra croce poron
  fPhysiScintCFLatMPO = new G4PVPlacement(0,                   
					  G4ThreeVector(-fPhysiCFLatPO_X,fPhysiCFLatPO_Y,fPhysiCFLatPO_Z),     
					  "CFSupportLMsmall",              
					  fLogicCFLatPO,       
					  fPhysiScintLayer,        
					  false,               
					  0,true); 
  
  // poron parete laterale sinistra croce poron
  fPhysiScintPoronLatPPO = new G4PVPlacement(0,
					     G4ThreeVector(fPhysiPoronLatPO_X,fPhysiPoronLatPO_Y,fPhysiPoronLatPO_Z),     
					     "PoronLPsmall",              
					     fLogicPoronLatPO,       
					     fPhysiScintLayer,        
					     false,               
					     0,true); 

  // poron parete laterla destra croce poron				  
  fPhysiScintPoronLatMPO = new G4PVPlacement(0,                   
					     G4ThreeVector(-fPhysiPoronLatPO_X,fPhysiPoronLatPO_Y,fPhysiPoronLatPO_Z),     
					     "PoronLMsmall",              
					     fLogicPoronLatPO,       
					     fPhysiScintLayer,        
					     false,               
					     0,true); 
  // poron parete sopra croce poron
  fPhysiScintPoronLatUpPO = new G4PVPlacement(0,
					      G4ThreeVector(fPhysiPoronFrontPO_X,-fPhysiPoronFrontPO_Y,fPhysiPoronFrontPO_Z),     
					      "PoronFPsmall",              
					      fLogicPoronFrontPO,     
					      fPhysiScintLayer,        
					      false,               
					      0,true);
  // poron parete sotto croce poron
  fPhysiScintPoronLatDownPO = new G4PVPlacement(0,
						G4ThreeVector(fPhysiPoronFrontPO_X,fPhysiPoronFrontPO_Y,fPhysiPoronFrontPO_Z),     
						"PoronFMsmall",              
						fLogicPoronFrontPO,     
						fPhysiScintLayer,        
						false,               
						0,true);
  


  /*###################################    Ultimo piano di scintillatore plastico    ###################################################*/
  
  fPhysiScintActiveLayer = new G4PVPlacement(0,                   
					     G4ThreeVector(fPhysiActiveLayer_X,fPhysiActiveLayer_Y,fPhysiActiveLayer_Z-fPoronFrontPO_Z/2.),     
					     // G4ThreeVector(0 , 0 , 0),     
					     "ActiveLastLayerScint",              
					     fLogicScintActiveLayer,       
					     fPhysiLastScintLayer,        
					     false,               
					     0,true); 
  // croce di poron
  fPhysiScintCFSupp = new G4PVPlacement(0,    
					G4ThreeVector(fPhysiCFSuppPor_X, fPhysiCFSuppPor_Y, fPhysiCFSuppPor_Z-fPoronFrontPO_Z/2. ),     
					"PORSupport_1",              
					fLogicCFSuppPoron,       
					fPhysiLastScintLayer,        
					false,               
					0,true); 
  // croce di carbonio
  fPhysiScintCFSupp = new G4PVPlacement(0,    
					G4ThreeVector(fPhysiCFSupp_X,fPhysiCFSupp_Y,fPhysiCFSupp_Z-fPoronFrontPO_Z/2.),     
					"CFSupport",              
					fLogicCFSupp,       
					fPhysiLastScintLayer,        
					false,               
					0,true); 

  // blocco carbonio sopra piano scint
  fPhysiScintCFFrontP = new G4PVPlacement(0,                   
					  G4ThreeVector(fPhysiCFFront_X,fPhysiCFFront_Y,fPhysiCFFront_Z-fPoronFrontPO_Z/2.),     
					  "CFSupportFP",              
					  fLogicCFFront,       
					  fPhysiLastScintLayer,        
					  false,               
					  0,true); 
  // blocco carbonio sotto piano scint
  fPhysiScintCFFrontM = new G4PVPlacement(0,                   
					  G4ThreeVector(fPhysiCFFront_X,-fPhysiCFFront_Y,fPhysiCFFront_Z-fPoronFrontPO_Z/2.),     
					  "CFSupportFM",              
					  fLogicCFFront,       
					  fPhysiLastScintLayer,        
					  false,               
					  0,true); 
  // blocco carbonio destra piano scint
  fPhysiScintCFLatP = new G4PVPlacement(0,                   
					G4ThreeVector(fPhysiCFLat_X,fPhysiCFLat_Y,fPhysiCFLat_Z-fPoronFrontPO_Z/2.),     
					"CFSupportLP",              
					fLogicCFLat,       
					fPhysiLastScintLayer,        
					false,               
					0,true); 

  // blocco carbonio sinistra piano scint
  fPhysiScintCFLatM = new G4PVPlacement(0,                   
					G4ThreeVector(-fPhysiCFLat_X,fPhysiCFLat_Y,fPhysiCFLat_Z-fPoronFrontPO_Z/2.),     
					"CFSupportLM",              
					fLogicCFLat,       
					fPhysiLastScintLayer,        
					0,true); 
  
  // poron parete laterale, destra piano scint 					  
  fPhysiScintPoronLatM = new G4PVPlacement(0,                   
					   G4ThreeVector(-fPhysiPoronLat_X, fPhysiPoronLat_Y, fPhysiPoronLat_Z-fPoronFrontPO_Z/2.),     
					   "PoronLM",              
					   fLogicPoronLat,       
					   fPhysiLastScintLayer,        
					   false,               
					   0,true); 
  // poron parete sotto piano scint
  fPhysiScintPoronLatDown = new G4PVPlacement(0,
					      G4ThreeVector(fPhysiPoronFront_X, fPhysiPoronFront_Y, fPhysiPoronFront_Z-fPoronFrontPO_Z/2.),     
					      "PoronFP",              
					      fLogicPoronFront,     
					      fPhysiLastScintLayer,        
					      false,               
					      0,true);
  // poron parete laterale, sinistra piano scint 					  
  fPhysiScintPoronLatP = new G4PVPlacement(0,                   
					   G4ThreeVector(fPhysiPoronLat_X, fPhysiPoronLat_Y, fPhysiPoronLat_Z-fPoronFrontPO_Z/2.),     
					   "PoronLP",              
					   fLogicPoronLat,       
					   fPhysiLastScintLayer,        
					   false,               
					   0,true); 
  // poron parete sopra piano scint
  fPhysiScintPoronLatUp = new G4PVPlacement(0,
					    G4ThreeVector(-fPhysiPoronFront_X, -fPhysiPoronFront_Y, fPhysiPoronFront_Z-fPoronFrontPO_Z/2.),     
					    "PoronFM",              
					    fLogicPoronFront,     
					    fPhysiLastScintLayer,        
					    false,               
					    0,true);
 

  /*###############################################################################################################################*/

  /* ################################     Cristalli di LYSO     ##################################*/
  
  fPhysiCrystalBox = new G4PVPlacement(0,                        
				       G4ThreeVector(fPhysiCrystalBox_X,fPhysiCrystalBox_Y,fPhysiCrystalBox_Z),   
				       "CalorimeterCrystal",               
				       fLogicCrystalBox,                
				       fPhysiCaloBox,                
				       false,                       
				       0,true);      
 
  fPhysiCFCrystalPanelDown = new G4PVPlacement(0, // piano carbonio dopo LYSO
					       G4ThreeVector(fPhysiCFCrystalPanelDown_X,fPhysiCFCrystalPanelDown_Y,fPhysiCFCrystalPanelDown_Z),   
					       "CFSupportPlaneDown",
					       fLogicCFCrystalPanelDown,
					       fPhysiCrystalBox,
					       false,
					       0,true);

  fPhysiCFCrystalPanelDown2 = new G4PVPlacement(0, // secondo piano carbonio dopo LYSO
						G4ThreeVector(fPhysiCFCrystalPanelDown2_X,fPhysiCFCrystalPanelDown2_Y, fPhysiCFCrystalPanelDown2_Z + ShiftOrigin),   
						"CFSupportPlaneDown2",
						fLogicCFCrystalPanelDown2,
						motherVolume,
						false,
						0,true);
 

  fPhysiCFCrystalSideXP = new G4PVPlacement(0,  // costole intorno al LYSO
					    G4ThreeVector(fPhysiCFCrystalSideX_X,fPhysiCFCrystalSideX_Y,fPhysiCFCrystalSideX_Z),
					    "CFSuportSideXP",
					    fLogicCFCrystalSideX,
					    fPhysiCrystalBox,
					    false,
					    0,true);
					    
  fPhysiCFCrystalSideXM = new G4PVPlacement(0,  // costole intorno al LYSO
					    G4ThreeVector(-fPhysiCFCrystalSideX_X,fPhysiCFCrystalSideX_Y,fPhysiCFCrystalSideX_Z),
					    "CFSuportSideXM",
					    fLogicCFCrystalSideX,
					    fPhysiCrystalBox,
					    false,
					    0,true);

  fPhysiCFCrystalSideYBigP = new G4PVPlacement(0, // costole intorno al LYSO
					       G4ThreeVector(fPhysiCFCrystalSideYBig_X,fPhysiCFCrystalSideYBig_Y,fPhysiCFCrystalSideYBig_Z),
					       "CFSuportSideYBigP",
					       fLogicCFCrystalSideYBig,
					       fPhysiCrystalBox,
					       false,
					       0,true);
					    
  fPhysiCFCrystalSideYBigM = new G4PVPlacement(0,  // costole intorno al LYSO
					       G4ThreeVector(-fPhysiCFCrystalSideYBig_X,-fPhysiCFCrystalSideYBig_Y,fPhysiCFCrystalSideYBig_Z),
					       "CFSuportSideYBigM",
					       fLogicCFCrystalSideYBig,
					       fPhysiCrystalBox,
					       false,
					       0,true);
					      
  fPhysiCrystalBlockPlaneContainer = new G4PVPlacement(0,
						       G4ThreeVector(fPhysiCrystalBlockPlaneContainer_X,fPhysiCrystalBlockPlaneContainer_Y,fPhysiCrystalBlockPlaneContainer_Z),
						       "CristalBlockPlaneContainer",
						       fLogicCrystalBlockPlaneContainer,
						       fPhysiCrystalBox,
						       false,
						       0,true);
  fPhysiScintCFSupp = new G4PVPlacement(0,  // croce poron attaccata al LYSO
					G4ThreeVector(fPhysiLastSuppPor_X, fPhysiLastSuppPor_Y, fPhysiLastSuppPor_Z ),     
					"PORLastSupport",              
					fLogicCFSuppPoron,       
					fPhysiCrystalBox,        
					false,               
					0,true); 
  				        

  fPhysiCrystalBlockRaw = new G4PVReplica("CrystalBlockRaw", 
					  fLogicCrystalBlockRawContainer,            
					  fPhysiCrystalBlockPlaneContainer,
					  kYAxis,                 
					  3,
					  fCrystalBlockRawContainer_Y);   

  fPhysiTeflonLYSO = new G4PVReplica("CrystalBlock", 
				     fLogicTeflonLYSO,            
				     fPhysiCrystalBlockRaw,
				     kXAxis,                 
				     3,
				     fCFBlockContainerExt_X);   

  
  fPhysiCrystalActiveBlock = new G4PVPlacement(0,                   
					       G4ThreeVector(0 , 0 ,0),
					       "ActiveBlockCrystal",
					       fLogicCrystalActiveBlock,       
					       fPhysiTeflonLYSO,        
					       false,               
					       0,true); 
  fPhysiTeflonContainer = new G4PVPlacement(0,                   
					    G4ThreeVector(0,0,0),
					    "Teflon",              
					    fLogicTeflonContainer,       
					    fPhysiTeflonLYSO,        
					    false,               
					    0,true); 

  

  /*#####################################################################################################################*/
  

    
  /*#####################################################################################################################*/
  
  /*#########################       VETO Detector      #############################*/
 
  
  fPhysiPoronLatXRight_3 = new G4PVPlacement(0,
					     G4ThreeVector(fPhysiPoronLatX3_X,fPhysiPoronLatX3_Y,fPhysiPoronLatX3_Z+ ShiftOrigin),
					     "VetoPoronXRight_Ext",              
					     fLogicPoronLatX,       
					     motherVolume,        
					     false,               
					     0,true);
  
  
  fPhysiPoronLatXLeft_3 = new G4PVPlacement(_rot180Z,
					    G4ThreeVector(-fPhysiPoronLatX3_X, -fPhysiPoronLatX3_Y,fPhysiPoronLatX3_Z+ ShiftOrigin),
					    "VetoPoronXLeft_Ext",              
					    fLogicPoronLatX,       
					    motherVolume,        
					    false,               
					    0,true);


  fPhysiPoronLatYUp_3 = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiPoronLatY3_X,fPhysiPoronLatY3_Y,fPhysiPoronLatY3_Z+ ShiftOrigin),
					  "VetoPoronYUp_Ext",              
					  fLogicPoronLatY,       
					  motherVolume,        
					  false,               
					  0,true);


  fPhysiPoronLatYDown_3 = new G4PVPlacement(_rot180Z,
					    G4ThreeVector(-fPhysiPoronLatY3_X,-fPhysiPoronLatY3_Y,fPhysiPoronLatY3_Z+ ShiftOrigin),
					    "VetoPoronYDown_Ext",              
					    fLogicPoronLatY,       
					    motherVolume,        
					    false,               
					    0,true);




  ///


  fPhysiPoronLatXRightInt_3 = new G4PVPlacement(0,
						G4ThreeVector(fPhysiPoronLatX3Int_X,fPhysiPoronLatX3Int_Y,fPhysiPoronLatX3_Z+ ShiftOrigin),
						"VetoPoronXRight_Int",              
						fLogicPoronLatX,       
						motherVolume,        
						false,               
						0,true);
 
 
  fPhysiPoronLatXLeftInt_3 = new G4PVPlacement(_rot180Z,
					       G4ThreeVector(-fPhysiPoronLatX3Int_X,-fPhysiPoronLatX3Int_Y,fPhysiPoronLatX3_Z+ ShiftOrigin),
					       "VetoPoronXLeft_Int",              
					       fLogicPoronLatX,       
					       motherVolume,        
					       false,               
					       0,true);


  fPhysiPoronLatYUpInt_3 = new G4PVPlacement(0,
					     G4ThreeVector(fPhysiPoronLatY3Int_X,fPhysiPoronLatY3Int_Y,fPhysiPoronLatY3_Z+ ShiftOrigin),
					     "VetoPoronYUp_Int",              
					     fLogicPoronLatY,       
					     motherVolume,        
					     false,               
					     0,true);
 
  fPhysiPoronLatYDownInt_3 = new G4PVPlacement(_rot180Z,
					       G4ThreeVector(-fPhysiPoronLatY3Int_X,-fPhysiPoronLatY3Int_Y,fPhysiPoronLatY3_Z+ ShiftOrigin),
					       "VetoPoronYDown_Int",              
					       fLogicPoronLatY,       
					       motherVolume,        
					       false,               
					       0,true);
 
 
  ///
 
 
  fPhysiVetoLatXRight_2 = new G4PVPlacement(0,
					    G4ThreeVector(-fPhysiVetoLatX2_X, fPhysiVetoLatX2_Y,fPhysiVetoLatX2_Z+ ShiftOrigin),
					    "VetoXRight",              
					    fLogicVetoLatX,       
					    motherVolume,        
					    false,               
					    0,true);
 
  fPhysiVetoLatXLeft_2 = new G4PVPlacement(_rot180Z,
					   G4ThreeVector(fPhysiVetoLatX2_X, -fPhysiVetoLatX2_Y, fPhysiVetoLatX2_Z+ ShiftOrigin),
					   "VetoXLeft",              
					   fLogicVetoLatX2,       
					   motherVolume,        
					   false,               
					   0,true);



  fPhysiVetoLatYUp_2 = new G4PVPlacement(0,
					 G4ThreeVector(fPhysiVetoLatY2_X, -fPhysiVetoLatY2_Y,fPhysiVetoLatY2_Z+ ShiftOrigin),
					 "VetoYUp",              
					 fLogicVetoLatY,       
					 motherVolume,        
					 false,               
					 0,true);

  fPhysiVetoLatYDown_2 = new G4PVPlacement(_rot180Z,
					   G4ThreeVector(-fPhysiVetoLatY2_X, fPhysiVetoLatY2_Y,fPhysiVetoLatY2_Z+ ShiftOrigin),
					   "VetoYDown",              
					   fLogicVetoLatY2,       
					   motherVolume,        
					   false,               
					   0,true);


  fPhysiCFVetoLatXRight_2 = new G4PVPlacement(0,
					      G4ThreeVector(fPhysiCFVetoLatX2_X, fPhysiCFVetoLatX2_Y, fPhysiCFVetoLatX2_Z+ ShiftOrigin),
					      "VetoCFXRight",              
					      fLogicCFVetoLatX,       
					      motherVolume,        
					      false,               
					      0,true);
 
  fPhysiCFVetoLatXLeft_2 = new G4PVPlacement(_rot180Z,
					     G4ThreeVector(-fPhysiCFVetoLatX2_X, -fPhysiCFVetoLatX2_Y, fPhysiCFVetoLatX2_Z+ ShiftOrigin),
					     "VetoCFXLeft",              
					     fLogicCFVetoLatX,       
					     motherVolume,        
					     false,               
					     0,true);
 

  fPhysiCFVetoLatYUp_2 = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiCFVetoLatY2_X,fPhysiCFVetoLatY2_Y,fPhysiCFVetoLatY2_Z+ ShiftOrigin),
					   "VetoCFYUp",              
					   fLogicCFVetoLatY,       
					   motherVolume,        
					   false,               
					   0,true);

  fPhysiCFVetoLatYDown_2 = new G4PVPlacement(_rot180Z,
					     G4ThreeVector(-fPhysiCFVetoLatY2_X,-fPhysiCFVetoLatY2_Y,fPhysiCFVetoLatY2_Z+ ShiftOrigin),
					     "VetoCFYDown",              
					     fLogicCFVetoLatY,       
					     motherVolume,        
					     false,               
					     0,true);
  

  ////////////////////////////// Bottom VETO
  fPhysiRealTrapPoronLayer1 = new G4PVPlacement(_rot45,                   
						G4ThreeVector(85.3583951/2.*mm,85.3583951/2.*mm,fPhysiRealTrapPoron1_Z+ ShiftOrigin),     
						"VETOBottPoronLayer1",              
						fLogicRealTrapPoronLayer,       
						motherVolume,        
						false,               
						0,true); 

  fPhysiRealTrapVetoLayer = new G4PVPlacement(_rot45,                
					      G4ThreeVector(85.3583951/2.*mm,85.3583951/2.*mm,fPhysiRealTrapVeto_Z+ ShiftOrigin),     
					      "VETOBottScintLayer",              
					      fLogicRealTrapVetoLayer,       
					      motherVolume,        
					      false,               
					      0,true); 

  fPhysiRealTrapPoronLayer2 = new G4PVPlacement(_rot45,                   
						G4ThreeVector(85.3583951/2.*mm,85.3583951/2.*mm,fPhysiRealTrapPoron2_Z+ ShiftOrigin),     
						"VETOBottPoronLayer2",              
						fLogicRealTrapPoronLayer,       
						motherVolume,        
						false,               
						0,true); 

  fPhysiRealTrapCFLayer = new G4PVPlacement(_rot45,                   
					    G4ThreeVector(85.3583951/2.*mm,85.3583951/2.*mm,fPhysiRealTrapCF_Z+ ShiftOrigin),     
					    "VETOBottCFLayer",              
					    fLogicRealTrapCFLayer,       
					    motherVolume,        
					    false,               
					    0,true); 

  fPhysiRealTrapCFLayer2 = new G4PVPlacement(_rot45,                   
					    G4ThreeVector(85.3583951/2.*mm,85.3583951/2.*mm,fPhysiRealTrapCF2_Z+ ShiftOrigin),     
					    "VETOBottCFLayer2",              
					    fLogicRealTrapCFLayer,       
					    motherVolume,        
					    false,               
					    0,true); 


  ////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////


  
  fPhysiPoronPlateV1_Top = new G4PVPlacement(0,
					     G4ThreeVector(fPhysiPoronPlateVTop_X, fPhysiPoronPlateVTop_Y, fPhysiPoronPlateVTop_Z+ ShiftOrigin),
					     "PoronPlateV_Top",              
					     fLogicPoronPlateV,
					     motherVolume,        
					     false,               
					     0,true);  
  
  fPhysiPoronPlateO1_Top = new G4PVPlacement(0,
					     G4ThreeVector(fPhysiPoronPlateOTop_X, fPhysiPoronPlateOTop_Y, fPhysiPoronPlateOTop_Z+ ShiftOrigin),
					     "PoronPlateO_Top",              
					     fLogicPoronPlateO,
					     motherVolume,        
					     false,               
					     0,true);  
  

  fPhysiPoronPlateV2_Top = new G4PVPlacement(0,
					     G4ThreeVector(-fPhysiPoronPlateVTop_X, -fPhysiPoronPlateVTop_Y, fPhysiPoronPlateVTop_Z+ ShiftOrigin),
					     "PoronPlateV_Top",              
					     fLogicPoronPlateV,
					     motherVolume,        
					     false,               
					     0,true);  
  
  fPhysiPoronPlateO2_Top = new G4PVPlacement(0,
					     G4ThreeVector(-fPhysiPoronPlateOTop_X, -fPhysiPoronPlateOTop_Y, fPhysiPoronPlateOTop_Z+ ShiftOrigin),
					     "PoronPlateO_Top",              
					     fLogicPoronPlateO,
					     motherVolume,        
					     false,               
					     0,true);  


  fPhysiPoronPlateV1_Bottom = new G4PVPlacement(0,
						G4ThreeVector(fPhysiPoronPlateVBottom_X, fPhysiPoronPlateVBottom_Y, fPhysiPoronPlateVBottom_Z+ ShiftOrigin),
						"PoronPlateV_Bottom",              
						fLogicPoronPlateV,
						motherVolume,        
						false,               
						0,true);  
  
  fPhysiPoronPlateO1_Bottom = new G4PVPlacement(0,
						G4ThreeVector(fPhysiPoronPlateOBottom_X, fPhysiPoronPlateOBottom_Y, fPhysiPoronPlateOBottom_Z+ ShiftOrigin),
						"PoronPlateO_Bottom",              
						fLogicPoronPlateO,
						motherVolume,        
						false,               
						0,true);  
  
  fPhysiPoronPlateV2_Bottom = new G4PVPlacement(0,
						G4ThreeVector(-fPhysiPoronPlateVBottom_X, -fPhysiPoronPlateVBottom_Y, fPhysiPoronPlateVBottom_Z+ ShiftOrigin),
						"PoronPlateV_Bottom",              
						fLogicPoronPlateV,
						motherVolume,        
						false,               
						0,true);  
  
  fPhysiPoronPlateO2_Bottom = new G4PVPlacement(0,
						G4ThreeVector(-fPhysiPoronPlateOBottom_X, -fPhysiPoronPlateOBottom_Y, fPhysiPoronPlateOBottom_Z+ ShiftOrigin),
						"PoronPlateO_Bottom",              
						fLogicPoronPlateO,
						motherVolume,        
						false,               
						0,true);  


  
  /*
  
 
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
  G4ThreeVector(fPhysiHoneyCombLatX_X,-fPhysiHoneyCombLatX_Y,fPhysifLogicCrystalActiveBlock->SetVisAttributes(attGreen);HoneyCombLatX_Z),
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
  fLogicCrystalActiveBlock->SetVisAttributes(attGreen);
  fPhysiHoneyCombSkinBottomD = new G4PVPlacement(0,
  G4ThreeVector(fPhysiHoneyCombSkinBottomD_X,fPhysiHoneyCombSkinBottomD_Y,fPhysiHoneyCombSkinBottomD_Z),
  "HoneyCombSkinBD",
  fLogicHoneyCombSkinBottom,
  motherVolume,
  false,
  0,true);
  */

 
  fPhysiCFSuppTop = new G4PVPlacement(0,      // prima croce all'inizio del calorimetro
				      G4ThreeVector(fPhysiCFSuppTop_X,fPhysiCFSuppTop_Y,fPhysiCFSuppTop_Z+ ShiftOrigin),
				      "CFSupportTop",
				      fLogicCFSuppTop, 
 				      motherVolume,
				      false,
				      0,true);


  fPhysiCFSuppTop = new G4PVPlacement(0,      // prima croce poron all'inizio del calorimetro
				      G4ThreeVector(fPhysiPORSuppTop_X, fPhysiPORSuppTop_Y, fPhysiPORSuppTop_Z+ ShiftOrigin),
				      "PORSupportTop",
				      fLogicCFSuppPoron, 
 				      motherVolume,
				      false,
				      0,true);

  // blocco carbonio sopra croce poron
  fPhysiScintCFFrontPPO = new G4PVPlacement(0,                   
					    G4ThreeVector(fPhysiCFFrontPO_X, fPhysiCFFrontPO_Y, fPhysiPORSuppTop_Z+ ShiftOrigin), 
					    "CFSupportFPsmall",              
					    fLogicCFFrontPO,       
					    motherVolume,        
					    false,               
					    0,true); 
  // blocco carbonio sotto croce  poron
  fPhysiScintCFFrontMPO = new G4PVPlacement(0,                   
					    G4ThreeVector(fPhysiCFFrontPO_X,-fPhysiCFFrontPO_Y, fPhysiPORSuppTop_Z+ ShiftOrigin), 
					    "CFSupportFMsmall",              
					    fLogicCFFrontPO,       
					    motherVolume,        
					    false,               
					    0,true); 
  // blocco carbonio destra croce poron
  fPhysiScintCFLatPPO = new G4PVPlacement(0,                   
					  G4ThreeVector(fPhysiCFLatPO_X,fPhysiCFLatPO_Y, fPhysiPORSuppTop_Z+ ShiftOrigin), 
					  "CFSupportLPsmall",              
					  fLogicCFLatPO,       
					  motherVolume,        
					  false,               
					  0,true); 
  // blocco carbonio sinistra croce poron
  fPhysiScintCFLatMPO = new G4PVPlacement(0,                   
					  G4ThreeVector(-fPhysiCFLatPO_X,fPhysiCFLatPO_Y, fPhysiPORSuppTop_Z+ ShiftOrigin), 
					  "CFSupportLMsmall",              
					  fLogicCFLatPO,       
					  motherVolume,        
					  false,               
					  0,true); 
  
  // poron parete laterale sinistra croce poron
  fPhysiScintPoronLatPPO = new G4PVPlacement(0,
					     G4ThreeVector(fPhysiPoronLatPO_X,fPhysiPoronLatPO_Y, fPhysiPORSuppTop_Z+ ShiftOrigin), 
					     "PoronLPsmall",              
					     fLogicPoronLatPO,       
					     motherVolume,        
					     false,               
					     0,true); 
  
  // poron parete laterla destra croce poron				  
  fPhysiScintPoronLatMPO = new G4PVPlacement(0,                   
					     G4ThreeVector(-fPhysiPoronLatPO_X,fPhysiPoronLatPO_Y, fPhysiPORSuppTop_Z+ ShiftOrigin), 
					     "PoronLMsmall",              
					     fLogicPoronLatPO,       
					     motherVolume,        
					     false,               
					     0,true); 
  // poron parete sopra croce poron
  fPhysiScintPoronLatUpPO = new G4PVPlacement(0,
					      G4ThreeVector(fPhysiPoronFrontPO_X,-fPhysiPoronFrontPO_Y, fPhysiPORSuppTop_Z+ ShiftOrigin), 
					      "PoronFPsmall",              
					      fLogicPoronFrontPO,     
					      motherVolume,        
					      false,               
					      0,true);
  // poron parete sotto croce poron
  fPhysiScintPoronLatDownPO = new G4PVPlacement(0,
						G4ThreeVector(fPhysiPoronFrontPO_X,fPhysiPoronFrontPO_Y, fPhysiPORSuppTop_Z+ ShiftOrigin), 
						"PoronFMsmall",              
						fLogicPoronFrontPO,     
						motherVolume,        
						false,               
						0,true);

  //  ///////////////////////////////////////////////////////////////////////////////////



  //Visualization Attribute
  G4VisAttributes* attCyan = new G4VisAttributes(G4Colour::Cyan());
  attCyan->SetVisibility(true);
  attCyan->SetForceAuxEdgeVisible(true);
  
  
  G4VisAttributes* attMagenta = new G4VisAttributes(G4Colour::Magenta());
  attMagenta->SetVisibility(true);
  attMagenta->SetForceAuxEdgeVisible(true);
  fLogicS1Scint->SetVisAttributes(attMagenta);
  fLogicScintActiveLayer->SetVisAttributes(attMagenta);

    
  G4VisAttributes* attYellow = new G4VisAttributes(G4Colour::Yellow());
  attYellow->SetVisibility(true);
  attYellow->SetForceAuxEdgeVisible(true); 
  fLogicS1SuppPoronFront->SetVisAttributes(attYellow);
  fLogicPoronLat->SetVisAttributes(attYellow);
  fLogicPoronFront->SetVisAttributes(attYellow);
  fLogicPoronLatPO->SetVisAttributes(attYellow);
  fLogicPoronFrontPO->SetVisAttributes(attYellow);
  fLogicPoronLatX->SetVisAttributes(attYellow);
  fLogicPoronLatY->SetVisAttributes(attYellow);
  fLogicCFSuppPoron->SetVisAttributes(attYellow);
  fLogicRealTrapPoronLayer->SetVisAttributes(attYellow);
  fLogicPoronPlateO->SetVisAttributes(attYellow);
  fLogicPoronPlateV->SetVisAttributes(attYellow);
  fLogicExternalPoronSupport->SetVisAttributes(attYellow);


  G4VisAttributes* attBrown = new G4VisAttributes(G4Colour::Brown());
  attBrown->SetVisibility(true);
  attBrown->SetForceAuxEdgeVisible(true);
  fLogicS1SuppBack->SetVisAttributes(attBrown);
  fLogicS1SuppFront->SetVisAttributes(attBrown);
  fLogicS1SuppHoleBar->SetVisAttributes(attBrown);
  fLogicS1SuppBar->SetVisAttributes(attBrown);
  fLogicS1SuppThinBar->SetVisAttributes(attBrown);
  fLogicS1SuppThinBarBack->SetVisAttributes(attBrown);
  fLogicCFLat->SetVisAttributes(attBrown);
  fLogicCFFront->SetVisAttributes(attBrown); 
  fLogicCFLatPO->SetVisAttributes(attBrown);
  fLogicCFFrontPO->SetVisAttributes(attBrown);
  fLogicCFSupp->SetVisAttributes(attBrown);
  fLogicCFSuppTop->SetVisAttributes(attBrown);
  fLogicCFCrystalPanelDown->SetVisAttributes(attBrown);
  fLogicCFBlockContainer->SetVisAttributes(attBrown);
  fLogicCFCrystalSideX->SetVisAttributes(attBrown);
  fLogicCFCrystalSideYBig->SetVisAttributes(attBrown);
  fLogicRealTrapCFLayer->SetVisAttributes(attBrown);
  fLogicCFCrystalPanelDown2->SetVisAttributes(attBrown);

  fLogicCFVetoLatX->SetVisAttributes(attBrown);
  fLogicCFVetoLatY->SetVisAttributes(attBrown);


  G4VisAttributes* attGreen = new G4VisAttributes(G4Colour::Green());
  attGreen->SetVisibility(true);
  attGreen->SetForceAuxEdgeVisible(true);
  fLogicTeflonContainer->SetVisAttributes(attGreen);


  G4VisAttributes* attRed = new G4VisAttributes(G4Colour::Red());
  attRed->SetVisibility(true);
  attRed->SetForceAuxEdgeVisible(true);
  


  G4VisAttributes* attBlue = new G4VisAttributes(G4Colour::Blue());
  attBlue->SetVisibility(true);
  attBlue->SetForceAuxEdgeVisible(true);
  fLogicCrystalActiveBlock->SetVisAttributes(attBlue);
  fLogicVetoLatX->SetVisAttributes(attBlue);
  fLogicVetoLatX2->SetVisAttributes(attBlue);
  fLogicVetoLatY->SetVisAttributes(attBlue);
  fLogicVetoLatY2->SetVisAttributes(attBlue);
  fLogicRealTrapVetoLayer->SetVisAttributes(attBlue);

 
  
  G4VisAttributes * attInvisible = new G4VisAttributes();
  attInvisible->SetVisibility(false);
  attInvisible->SetForceAuxEdgeVisible(false);
  fLogicS1->SetVisAttributes(attInvisible);  
  fLogicS1ScintContainerP->SetVisAttributes(attInvisible);
  fLogicS1ScintContainerM->SetVisAttributes(attInvisible);
  fLogicCaloBox->SetVisAttributes(attInvisible); 
  fLogicScintBox->SetVisAttributes(attInvisible); 
  fLogicScintLayer->SetVisAttributes(attInvisible); 
  fLogicLastScintLayer->SetVisAttributes(attInvisible); 
  fLogicCrystalBlockPlaneContainer->SetVisAttributes(attInvisible); 
  fLogicCrystalBlockRawContainer->SetVisAttributes(attInvisible); 
  fLogicCrystalBlockContainer->SetVisAttributes(attInvisible); 
  fLogicCrystalBox->SetVisAttributes(attInvisible); 

}


