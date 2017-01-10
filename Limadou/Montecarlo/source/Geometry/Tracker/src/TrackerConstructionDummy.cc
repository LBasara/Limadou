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

#include "TrackerConstructionDummy.hh"

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

TrackerConstructionDummy::TrackerConstructionDummy()
  :fSolidSiSensor(0),fSolidActiveLadder(0),fSolidLadderPlane(0),fSolidLadder(0),fSolid(0),
   fLogicSiSensor(0),
   fLogicActiveLadderL1(0),fLogicLadderPlaneL1(0),fLogicLadderL1(0),fLogicL1(0),
   fLogicActiveLadderL2(0),fLogicLadderPlaneL2(0),fLogicLadderL2(0),fLogicL2(0),
   fPhysiL1(0),fPhysiLadderPlaneL1(0),fPhysiLadderL1(0),fPhysiActiveLadderL1(0),fPhysiSiSensorL1(0),
   fPhysiL2(0),fPhysiLadderPlaneL2(0),fPhysiLadderL2(0),fPhysiActiveLadderL2(0),fPhysiSiSensorL2(0)
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
  

  fLadder_X = fActiveLadder_X;
  fLadder_Y = fActiveLadder_Y+5.5*mm;
  fLadder_Z = fActiveLadder_Z;

  fLadderPlane_X = fLadder_X;
  fLadderPlane_Y = fLadder_Y*fLadderNumber;
  fLadderPlane_Z = fLadder_Z;

  fL_X = fLadderPlane_X;
  fL_Y = fLadderPlane_Y;
  fL_Z = fLadderPlane_Z;

  ComputeObjectsPositioning();

  siliconMaterial = "Silicon";
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerConstructionDummy::~TrackerConstructionDummy()
{
  if (pMaterial) delete pMaterial;
}

void TrackerConstructionDummy::ComputeObjectsPositioning(){

  fPhysiL1_X = 0;
  fPhysiL1_Y = 0;
  fPhysiL1_Z = 115*mm+10*mm+fL_Z/2.;
  
  fPhysiL2_X = 0;
  fPhysiL2_Y = 0;
  fPhysiL2_Z = fPhysiL1_Z+10*mm+fL_Z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerConstructionDummy::Builder(G4VPhysicalVolume* motherVolume)
{
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String tkSDname = "/hepd/silicon";
  TrackerSD * tkSD = new TrackerSD(tkSDname);
  SDman->AddNewDetector(tkSD);

  pMaterial->DefineMaterials();
  G4Material* vacuum = pMaterial->GetMaterial("Galactic");
  G4Material* silMat  = pMaterial->GetMaterial(siliconMaterial);

  fSolidSiSensor = new G4Box("SiSensor",fSiSens_X/2.,fSiSens_Y/2.,fSiSens_Z/2.);
  fSolidActiveLadder = new G4Box("ActiveLadder",fActiveLadder_X/2.,fActiveLadder_Y/2.,fActiveLadder_Z/2.);
  fSolidLadder = new G4Box("Ladder",fLadder_X/2.,fLadder_Y/2.,fLadder_Z/2.);
  fSolidLadderPlane = new G4Box("LadderPlane",fLadderPlane_X/2.,fLadderPlane_Y/2.,fLadderPlane_Z/2.);
  
  fSolid = new G4Box("Layer",fL_X/2.,fL_Y/2.,fL_Z/2.);
  
  fLogicSiSensor = new G4LogicalVolume(fSolidSiSensor,silMat,"SiSensor");

  fLogicActiveLadderL1 = new G4LogicalVolume(fSolidActiveLadder,vacuum,"ActiveLadderL1");
  fLogicLadderL1 = new G4LogicalVolume(fSolidLadder,vacuum,"LadderL1");
  fLogicLadderPlaneL1 = new G4LogicalVolume(fSolidLadderPlane,vacuum,"LadderPlaneL1");
  fLogicL1 = new G4LogicalVolume(fSolid,vacuum,"Layer1");

  fLogicActiveLadderL2 = new G4LogicalVolume(fSolidActiveLadder,vacuum,"ActiveLadderL2");
  fLogicLadderL2 = new G4LogicalVolume(fSolidLadder,vacuum,"LadderL2");
  fLogicLadderPlaneL2 = new G4LogicalVolume(fSolidLadderPlane,vacuum,"LadderPlaneL2");
  fLogicL2 = new G4LogicalVolume(fSolid,vacuum,"Layer2");

  fLogicSiSensor->SetSensitiveDetector(tkSD);

  fPhysiL1 = new G4PVPlacement(0,
			       G4ThreeVector(fPhysiL1_X,fPhysiL1_Y,fPhysiL1_Z),
			       "Layer1",
			       fLogicL1,
			       motherVolume,
			       false,
			       0,true);
    
  fPhysiLadderPlaneL1 = new G4PVPlacement(0,
					  G4ThreeVector(0,0,0),
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
  
  fPhysiActiveLadderL1 = new G4PVPlacement(0,
					   G4ThreeVector(0,0,0),
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
     
  fPhysiLadderPlaneL2 = new G4PVPlacement(0,
					  G4ThreeVector(0,0,0),
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
        
  fPhysiActiveLadderL2 = new G4PVPlacement(0,
					   G4ThreeVector(0,0,0),
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


