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

#include "TrackerConstructionConfig2.hh"

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
#include "TrackerSD.hh"


#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerConstructionConfig2::TrackerConstructionConfig2()
  :fSolidSiSensor(0),fSolidActiveLadder(0),fSolidKaptonNLadder(0),fSolidKaptonPLadder(0),fSolidLadderPlane(0),fSolidLadder(0),
   fSolidSupport(0),fSolidSupportHole(0),fSolidSupportT(0),fSolidSupportFinal(0),fSolid(0),
  fLogicSiSensor(0),
  fLogicActiveLadderL1(0),fLogicKaptonNLadderL1(0),fLogicKaptonPLadderL1(0),fLogicLadderPlaneL1(0),fLogicLadderL1(0),fLogicSupportL1(0),fLogicL1(0),
  fLogicActiveLadderL2(0),fLogicKaptonNLadderL2(0),fLogicKaptonPLadderL2(0),fLogicLadderPlaneL2(0),fLogicLadderL2(0),fLogicSupportL2(0),fLogicL2(0),
  fPhysiL1(0),fPhysiL1SupportU(0),fPhysiL1SupportD(0),fPhysiLadderPlaneL1(0),fPhysiLadderL1(0),fPhysiKaptonNLadderL1(0),
  fPhysiKaptonPLadderL1(0),fPhysiActiveLadderL1(0),fPhysiSiSensorL1(0),
  fPhysiL2(0),fPhysiL2SupportU(0),fPhysiL2SupportD(0),fPhysiLadderPlaneL2(0),fPhysiLadderL2(0),fPhysiKaptonNLadderL2(0),
  fPhysiKaptonPLadderL2(0),fPhysiActiveLadderL2(0),fPhysiSiSensorL2(0)
{
  pMaterial = new HEPDSWMaterial();

  fSiSensNumber = 5; 

  fLadderNumber = 3;

  fSiSens_X = 41.36*mm;
  fSiSens_Y = 72.045*mm;
  fSiSens_Z = 0.3*mm;

  fActiveLadder_X =fSiSens_X*fSiSensNumber;
  fActiveLadder_Y =fSiSens_Y;
  fActiveLadder_Z =fSiSens_Z;
  
  fKaptonNLadder_X = 232.428*mm;
  fKaptonNLadder_Y = 67.2*mm;
  fKaptonNLadder_Z = 0.3*mm;
  //fKaptonNLadder_Z = 0.05*mm;

  fKaptonPLadder_X = 29*mm;
  fKaptonPLadder_Y = 72*mm;
  fKaptonPLadder_Z = 0.3*mm;
  //fKaptonPLadder_Z = 0.05*mm;

  fLadder_X = fKaptonNLadder_X;
  fLadder_Y = fSiSens_Y+5.5*mm;
  fLadder_Z = fActiveLadder_Z+fKaptonNLadder_Z+fKaptonPLadder_Z;

  fLadderPlane_X = fLadder_X;
  fLadderPlane_Y = fLadder_Y*fLadderNumber;
  fLadderPlane_Z = fLadder_Z;

  fSupport_X = 226.98*mm;
  fSupport_Y = 247.135*mm;
  fSupport_Z = 2*mm;

  fSupportHole_X = 200.98*mm; 
  fSupportHole_Y = 67.795*mm;
  fSupportHole_Z = 3*mm;

  fL_X =242.428*mm;
  fL_Y =247.135*mm;
  fL_Z =5.6*mm;

  ComputeObjectsPositioning();

  siliconMaterial = "Silicon";
 
  cfiberMaterial  = "CarbonFiber";
  kaptonMaterial  = "Kapton";
  /*
  cfiberMaterial  = "Galactic";
  kaptonMaterial  = "Galactic";
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerConstructionConfig2::~TrackerConstructionConfig2()
{
  if (pMaterial) delete pMaterial;
}

void TrackerConstructionConfig2::ComputeObjectsPositioning(){

  fSupportHoleStep_Y = fSupportHole_Y+8*mm;
  fSupportHole_Ymin  = fSupport_Y/2.-13*mm-fSupportHole_Y/2.;

  fPhysiL1_X = 0;
  fPhysiL1_Y = 0;
  fPhysiL1_Z = 147.8*mm+35.78*mm+fL_Z/2.;
  
  fPhysiL1SupportU_X = -(fL_X/2.-fSupport_X/2.);
  fPhysiL1SupportU_Y = 0;
  fPhysiL1SupportU_Z = fL_Z/2-fSupport_Z/2.;

  fPhysiL1SupportD_X = -(fL_X/2.-fSupport_X/2.);
  fPhysiL1SupportD_Y = 0;
  fPhysiL1SupportD_Z = -fL_Z/2+fSupport_Z/2.;

  fPhysiLadderPlaneL1_X = fL_X/2.-fLadder_X/2.;
  fPhysiLadderPlaneL1_Y = 0;
  fPhysiLadderPlaneL1_Z = 0;

  fPhysiKaptonNLadderL1_X = 0;
  fPhysiKaptonNLadderL1_Y = 0;
  fPhysiKaptonNLadderL1_Z = -(fLadder_Z/2.-fKaptonNLadder_Z/2.);

  fPhysiKaptonPLadderL1_X = fLadder_X/2.-fKaptonPLadder_X/2.;
  fPhysiKaptonPLadderL1_Y = 0;
  fPhysiKaptonPLadderL1_Z = fLadder_Z/2.-fKaptonPLadder_Z/2.;

  fPhysiActiveLadderL1_X = -(fLadder_X/2.-fActiveLadder_X/2.);
  fPhysiActiveLadderL1_Y = 0;
  fPhysiActiveLadderL1_Z = 0;
  
  fPhysiL2_X = 0;
  fPhysiL2_Y = 0;
  fPhysiL2_Z = 147.8*mm+45.78*mm+fL_Z/2.;
  
  fPhysiL2SupportU_X = -(fL_X/2.-fSupport_X/2.);
  fPhysiL2SupportU_Y = 0;
  fPhysiL2SupportU_Z = fL_Z/2-fSupport_Z/2.;

  fPhysiL2SupportD_X = -(fL_X/2.-fSupport_X/2.);
  fPhysiL2SupportD_Y = 0;
  fPhysiL2SupportD_Z = -fL_Z/2+fSupport_Z/2.;

  fPhysiLadderPlaneL2_X = fL_X/2.-fLadder_X/2.;
  fPhysiLadderPlaneL2_Y = 0;
  fPhysiLadderPlaneL2_Z = 0;

  fPhysiKaptonNLadderL2_X = 0;
  fPhysiKaptonNLadderL2_Y = 0;
  fPhysiKaptonNLadderL2_Z = -(fLadder_Z/2.-fKaptonNLadder_Z/2.);

  fPhysiKaptonPLadderL2_X = fLadder_X/2.-fKaptonPLadder_X/2.;
  fPhysiKaptonPLadderL2_Y = 0;
  fPhysiKaptonPLadderL2_Z = fLadder_Z/2.-fKaptonPLadder_Z/2.;

  fPhysiActiveLadderL2_X = -(fLadder_X/2.-fActiveLadder_X/2.);
  fPhysiActiveLadderL2_Y = 0;
  fPhysiActiveLadderL2_Z = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerConstructionConfig2::Builder(G4VPhysicalVolume* motherVolume)
{
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String tkSDname = "/hepd/silicon";
  TrackerSD * tkSD = new TrackerSD(tkSDname);
  SDman->AddNewDetector(tkSD);

  pMaterial->DefineMaterials();
  G4Material* vacuum = pMaterial->GetMaterial("Galactic");
  G4Material* cfMat  = pMaterial->GetMaterial(cfiberMaterial);
  G4Material* kapMat  = pMaterial->GetMaterial(kaptonMaterial);
  G4Material* silMat  = pMaterial->GetMaterial(siliconMaterial);


  fSolidSiSensor = new G4Box("SiSensor",fSiSens_X/2.,fSiSens_Y/2.,fSiSens_Z/2.);
  fSolidActiveLadder = new G4Box("ActiveLadder",fActiveLadder_X/2.,fActiveLadder_Y/2.,fActiveLadder_Z/2.);
  fSolidKaptonNLadder = new G4Box("KaptonN",fKaptonNLadder_X/2.,fKaptonNLadder_Y/2.,fKaptonNLadder_Z/2.);
  fSolidKaptonPLadder = new G4Box("KaptonP",fKaptonPLadder_X/2.,fKaptonPLadder_Y/2.,fKaptonPLadder_Z/2.);
  fSolidLadderPlane = new G4Box("LadderPlane",fLadderPlane_X/2.,fLadderPlane_Y/2.,fLadderPlane_Z/2.);
  fSolidLadder = new G4Box("Ladder",fLadder_X/2.,fLadder_Y/2.,fLadder_Z/2.);
  fSolidSupport = new G4Box("Support",fSupport_X/2.,fSupport_Y/2.,fSupport_Z/2.);
  fSolidSupportHole = new G4Box("Hole",fSupportHole_X/2.,fSupportHole_Y/2.,fSupportHole_Z/2.);
  for(int i=0;i<fLadderNumber;i++){
    if(i==0){
      fSolidSupportT.push_back(new G4SubtractionSolid("Support",
						      fSolidSupport,
						      fSolidSupportHole,
						      0,
						      G4ThreeVector(0,-fSupportHole_Ymin+i*fSupportHoleStep_Y,0)));
    }else{
      fSolidSupportT.push_back(new G4SubtractionSolid("Support",
						      fSolidSupportT[i-1],
						      fSolidSupportHole,
						      0,
						      G4ThreeVector(0,-fSupportHole_Ymin+i*fSupportHoleStep_Y,0)));
    }
  }
  
  fSolidSupportFinal = fSolidSupportT.back();
  
  fSolidSupportT.clear();
  fSolidSupportT.resize(0);

  fSolid = new G4Box("Layer",fL_X/2.,fL_Y/2.,fL_Z/2.);
  
  
  fLogicSiSensor = new G4LogicalVolume(fSolidSiSensor,silMat,"SiSensor");
  fLogicActiveLadderL1 = new G4LogicalVolume(fSolidActiveLadder,vacuum,"ActiveLadderL1");
  fLogicKaptonNLadderL1= new G4LogicalVolume(fSolidKaptonNLadder,kapMat,"KaptonNL1");
  fLogicKaptonPLadderL1= new G4LogicalVolume(fSolidKaptonPLadder,kapMat,"KaptonPL1");
  fLogicLadderPlaneL1 = new G4LogicalVolume(fSolidLadderPlane,vacuum,"LadderPlaneL1");
  fLogicLadderL1 = new G4LogicalVolume(fSolidLadder,vacuum,"LadderL1");
  fLogicSupportL1 = new G4LogicalVolume(fSolidSupportFinal,cfMat,"SupportL1");
  fLogicL1 = new G4LogicalVolume(fSolid,vacuum,"Layer1");

  fLogicActiveLadderL2 = new G4LogicalVolume(fSolidActiveLadder,vacuum,"ActiveLadderL2");
  fLogicKaptonNLadderL2= new G4LogicalVolume(fSolidKaptonNLadder,kapMat,"KaptonNL2");
  fLogicKaptonPLadderL2= new G4LogicalVolume(fSolidKaptonPLadder,kapMat,"KaptonPL2");
  fLogicLadderPlaneL2 = new G4LogicalVolume(fSolidLadderPlane,vacuum,"LadderPlaneL2");
  fLogicLadderL2 = new G4LogicalVolume(fSolidLadder,vacuum,"LadderL2");
  fLogicSupportL2 = new G4LogicalVolume(fSolidSupportFinal,cfMat,"SupportL2");
  fLogicL2 = new G4LogicalVolume(fSolid,vacuum,"Layer2");


  fLogicSiSensor->SetSensitiveDetector(tkSD);

  fPhysiL1 = new G4PVPlacement(0,
			       G4ThreeVector(fPhysiL1_X,fPhysiL1_Y,fPhysiL1_Z),
			       "Layer1",
			       fLogicL1,
			       motherVolume,
			       false,
			       0,true);

  fPhysiL1SupportU = new G4PVPlacement(0,
				       G4ThreeVector(fPhysiL1SupportU_X,fPhysiL1SupportU_Y,fPhysiL1SupportU_Z),
				       "L1Support",
				       fLogicSupportL1,
				       fPhysiL1,
				       false,
				       0,true);

  fPhysiL1SupportD = new G4PVPlacement(0,
				       G4ThreeVector(fPhysiL1SupportD_X,fPhysiL1SupportD_Y,fPhysiL1SupportD_Z),
				       "L1Support",
				       fLogicSupportL1,
				       fPhysiL1,
				       false,
				       0,true);

    
  fPhysiLadderPlaneL1 = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiLadderPlaneL1_X,fPhysiLadderPlaneL1_Y,fPhysiLadderPlaneL1_Z),
					  "LadderPLaneL1",
					  fLogicLadderPlaneL1,
					  fPhysiL1,
					  false,
					  0,true);

  fPhysiLadderL1 = new G4PVReplica("LadderL1",
				   fLogicLadderL1,
				   fPhysiLadderPlaneL1,
				   kYAxis,
				   fLadderNumber,
				   fLadder_Y);
  
  fPhysiKaptonNLadderL1 = new G4PVPlacement(0,
					    G4ThreeVector(fPhysiKaptonNLadderL1_X,fPhysiKaptonNLadderL1_Y,fPhysiKaptonNLadderL1_Z),
					    "KaptonN",
					    fLogicKaptonNLadderL1,
					    fPhysiLadderL1,
					    false,
					    0,true);
  
  fPhysiKaptonPLadderL1 = new G4PVPlacement(0,
					    G4ThreeVector(fPhysiKaptonPLadderL1_X,fPhysiKaptonPLadderL1_Y,fPhysiKaptonPLadderL1_Z),
					    "KaptonP",
					    fLogicKaptonPLadderL1,
					    fPhysiLadderL1,
					    false,
					    0,true);
  
  fPhysiActiveLadderL1 = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiActiveLadderL1_X,fPhysiActiveLadderL1_Y,fPhysiActiveLadderL1_Z),
					   "ActiveLadder",
					   fLogicActiveLadderL1,
					   fPhysiLadderL1,
					   false,
					   0,true);
    
    
  fPhysiSiSensorL1 = new G4PVReplica("Sensor",
				     fLogicSiSensor,
				     fPhysiActiveLadderL1,
				     kXAxis,
				     fSiSensNumber,
				     fSiSens_X);
    

  
  fPhysiL2 = new G4PVPlacement(0,
			       G4ThreeVector(fPhysiL2_X,fPhysiL2_Y,fPhysiL2_Z),
			       "Layer2",
			       fLogicL2,
			       motherVolume,
			       false,
			       0,true);
    
  fPhysiL2SupportU = new G4PVPlacement(0,
				       G4ThreeVector(fPhysiL2SupportU_X,fPhysiL2SupportU_Y,fPhysiL2SupportU_Z),
				       "L2Support",
				       fLogicSupportL2,
				       fPhysiL2,
				       false,
				       0,true);
    
  fPhysiL2SupportD = new G4PVPlacement(0,
				       G4ThreeVector(fPhysiL2SupportD_X,fPhysiL2SupportD_Y,fPhysiL2SupportD_Z),
				       "L2Support",
				       fLogicSupportL2,
				       fPhysiL2,
				       false,
				       0,true);
    
  fPhysiLadderPlaneL2 = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiLadderPlaneL2_X,fPhysiLadderPlaneL2_Y,fPhysiLadderPlaneL2_Z),
					  "LadderPLaneL2",
					  fLogicLadderPlaneL2,
					  fPhysiL2,
					  false,
					  0,true);
    
  fPhysiLadderL2 = new G4PVReplica("LadderL2",
				   fLogicLadderL2,
				   fPhysiLadderPlaneL2,
				   kYAxis,
				   fLadderNumber,
				   fLadder_Y);
    
  fPhysiKaptonNLadderL2 = new G4PVPlacement(0,
					    G4ThreeVector(fPhysiKaptonNLadderL2_X,fPhysiKaptonNLadderL2_Y,fPhysiKaptonNLadderL2_Z),
					    "KaptonN",
					    fLogicKaptonNLadderL2,
					    fPhysiLadderL2,
					    false,
					    0,true);
    
  fPhysiKaptonPLadderL2 = new G4PVPlacement(0,
					    G4ThreeVector(fPhysiKaptonPLadderL2_X,fPhysiKaptonPLadderL2_Y,fPhysiKaptonPLadderL2_Z),
					    "KaptonP",
					    fLogicKaptonPLadderL2,
					    fPhysiLadderL2,
					    false,
					    0,true);
    
  fPhysiActiveLadderL2 = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiActiveLadderL2_X,fPhysiActiveLadderL2_Y,fPhysiActiveLadderL2_Z),
					   "ActiveLadder",
					   fLogicActiveLadderL2,
					   fPhysiLadderL2,
					   false,
					   0,true);
    
  
  fPhysiSiSensorL2 = new G4PVReplica("Sensor",
				     fLogicSiSensor,
				     fPhysiActiveLadderL2,
				     kXAxis,
				     fSiSensNumber,
				     fSiSens_X);
 
    
  //Visualization Attribute
  G4VisAttributes* attBlue = new G4VisAttributes(G4Colour::Blue());
  attBlue->SetVisibility(true);
  attBlue->SetForceAuxEdgeVisible(true);
  fLogicSiSensor->SetVisAttributes(attBlue);

  G4VisAttributes* attMagenta = new G4VisAttributes(G4Colour::Magenta());
  attMagenta->SetVisibility(true);
  attMagenta->SetForceAuxEdgeVisible(true);
  fLogicKaptonNLadderL1->SetVisAttributes(attMagenta);
  fLogicKaptonNLadderL2->SetVisAttributes(attMagenta);

  G4VisAttributes* attBrown = new G4VisAttributes(G4Colour::Brown());
  attBrown->SetVisibility(true);
  attBrown->SetForceAuxEdgeVisible(true);
  fLogicSupportL1->SetVisAttributes(attBrown);
  fLogicSupportL2->SetVisAttributes(attBrown);

  G4VisAttributes* attGreen = new G4VisAttributes(G4Colour::Green());
  attGreen->SetVisibility(true);
  attGreen->SetForceAuxEdgeVisible(true);
  fLogicKaptonPLadderL1->SetVisAttributes(attGreen);
  fLogicKaptonPLadderL2->SetVisAttributes(attGreen);

  G4VisAttributes * attInvisible = new G4VisAttributes();
  attInvisible->SetVisibility(false);
  attInvisible->SetForceAuxEdgeVisible(false);
  fLogicActiveLadderL1->SetVisAttributes(attInvisible);
  fLogicActiveLadderL2->SetVisAttributes(attInvisible);
  fLogicLadderL1->SetVisAttributes(attInvisible);
  fLogicLadderL2->SetVisAttributes(attInvisible);
  fLogicLadderPlaneL1->SetVisAttributes(attInvisible);
  fLogicLadderPlaneL2->SetVisAttributes(attInvisible);
  fLogicL1->SetVisAttributes(attInvisible);
  fLogicL2->SetVisAttributes(attInvisible);

}


