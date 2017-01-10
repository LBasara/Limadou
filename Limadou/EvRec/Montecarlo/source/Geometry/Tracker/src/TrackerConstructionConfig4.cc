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

#include "TrackerConstructionConfig4.hh"

#include "G4NistManager.hh"
#include "HEPDSWMaterial.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
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
#include "TrackerSD.hh"


#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerConstructionConfig4::TrackerConstructionConfig4()
  :fSolidTracker(0),fSolidLayer(0),fSolidCFFrameHole(0),fSolidCFFrame(0),
   fSolidPoronFrameHole(0),fSolidPoronFrame(0),fSolidSiliconPlate(0),fSolidSiliconSensor(0),
   fSolidLadderBox(0),fSolidLadder(0),fSolidKaptonS(0),fSolidKaptonK(0),fSolidHybrid(0),
   fSolidHeatSink(0),fSolidRingK(0),fSolidRingS(0),
   fLogicTracker(0),fLogicLayer(0),fLogicCFFrame(0),fLogicPoronFrame(0),fLogicSiliconPlateP(0),fLogicSiliconPlateM(0),
   fLogicSiliconSensor(0),fLogicLadderBox(0),fLogicLadder(0),fLogicKaptonS(0),fLogicKaptonK(0),
   fLogicHybrid(0),fLogicHeatSink(0),fLogicRingK(0),fLogicRingS(0),
   fPhysiTracker(0),fPhysiLayer(0),fPhysiCFFrame(0),fPhysiPoronFrame(0),fPhysiLadderBox(0),
   fPhysiLadder(0),fPhysiRingK(0),fPhysiKaptonK(0),fPhysiHybridK(0),fPhysiSiliconPlateP(0),
   fPhysiSiliconSensorP(0),fPhysiSiliconPlateM(0),fPhysiSiliconSensorM(0),fPhysiKaptonS(0),
  fPhysiHybridS(0),fPhysiRingS(0),fPhysiHeatSink(0)
{
  pMaterial = new HEPDSWMaterial();
  fLayerNumber = 2;
  fLadderNumber = 3;

  fTracker_X = 267.3*mm;
  fTracker_Y = 305.5*mm;
  fTracker_Z = 20.0*mm;

  fLayer_X = 267.3*mm;
  fLayer_Y = 305.5*mm;
  fLayer_Z = 10.0*mm;

  fCFFrame_X = 267.3*mm;
  fCFFrame_Y = 305.5*mm;
  fCFFrame_Z = 10.0*mm;
  
  fCFMiddleFrameHole_X = 255.3*mm;
  fCFMiddleFrameHole_Y = 305.5*mm;
  fCFMiddleFrameHole_Z = 6.3*mm;

  fCFFrameHole1_X = 247.9*mm;
  fCFFrameHole1_Y = 60.5*mm;
  fCFFrameHole1_Z = 6.0*mm;

  fCFFrameHole2_X = 72.6*mm;
  fCFFrameHole2_Y = 215.0*mm;
  fCFFrameHole2_Z = 11.0*mm;

  fPoronFrame_X = 255.3*mm;
  fPoronFrame_Y = 299.5*mm;
  fPoronFrame_Z = 6.3*mm;
	      
  fPoronMiddleFrameHole_X = 247.3*mm;
  fPoronMiddleFrameHole_Y = 299.5*mm;
  fPoronMiddleFrameHole_Z = 2.8*mm;

  fPoronFrameHole1_X = 247.9*mm;
  fPoronFrameHole1_Y = 66.5*mm;
  fPoronFrameHole1_Z = 7.0*mm;

  fPoronFrameHole2_X = 72.6*mm;
  fPoronFrameHole2_Y = 224.3*mm;
  fPoronFrameHole2_Z = 7.0*mm;

  fSiSens_X = 77.606*mm;
  fSiSens_Y = 109.65*mm;
  fSiSens_Z = 0.3*mm;

  fSiActiveSens_X = 71.606*mm;
  fSiActiveSens_Y = 105.85*mm;
  fSiActiveSens_Z = 0.3*mm;

  fLadderBox_X = 246.9*mm;
  fLadderBox_Y = 295.5*mm;
  fLadderBox_Z = 2.8*mm;

  fLadderBoxEnd_X = 246.9*mm;
  fLadderBoxEnd_Y = 60.5*mm;
  fLadderBoxEnd_Z = 4.41*mm;

  fLadder_X = 82.3*mm;
  fLadder_Y = 295.5*mm;
  fLadder_Z = 2.8*mm;
  
  fLadderEnd_X = 82.3*mm;
  fLadderEnd_Y = 60.5*mm;
  fLadderEnd_Z = 4.41*mm;
  
  fKaptonS_X = 72.0*mm;
  fKaptonS_Y = 9.01*mm;
  fKaptonS_Z = 2.4*mm;

  fKaptonK_X = 68.7*mm;
  fKaptonK_Y = 222.69*mm;
  fKaptonK_Z = 0.1*mm;

  fHybrid_X = 76*mm;
  fHybrid_Y = 60*mm;
  fHybrid_Z = 0.5*mm;

  fHeatSink_X = 78.3*mm;
  fHeatSink_Y = 60*mm;
  fHeatSink_Z = 2.32*mm;

  fRingK_X = 82.3*mm;
  fRingK_Y = 295.5*mm;
  fRingK_Z = 2.8*mm;

  fRingKHole1_X = 78.3*mm;
  fRingKHole1_Y = 295.5*mm;
  fRingKHole1_Z = 2*mm;

  fRingKHole2_X = 72.6*mm;
  fRingKHole2_Y = 215*mm;
  fRingKHole2_Z = 2.9*mm;

  fRingKHole3_X = 72.6*mm;
  fRingKHole3_Y = 63.5*mm;
  fRingKHole3_Z = 2.9*mm;

  fRingS_X = 78*mm;
  fRingS_Y = 221*mm;
  fRingS_Z = 1*mm;

  fRingSHole_X = 72.6*mm;
  fRingSHole_Y = 215*mm;
  fRingSHole_Z = 1.5*mm;

  ComputeObjectsPositioning();

  siliconMaterial = "Silicon";
  cfiberMaterial  = "CarbonFiber";
  kaptonMaterial  = "Kapton";
  poronMaterial   = "Poron";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerConstructionConfig4::~TrackerConstructionConfig4()
{
  if (pMaterial) delete pMaterial;
}

void TrackerConstructionConfig4::ComputeObjectsPositioning(){

  transCFFrameHole0_Y =-6*mm;
  transCFFrameHole0_Z =-0.15*mm;
  transCFFrameHole1_Y =-(fCFFrame_Y/2.-fCFFrameHole1_Y/2.);
  transCFFrameHole1_Z =-3*mm;
  transCFFrameHole2_X =82.5*mm;
  transCFFrameHole2_Y =fCFFrame_Y/2.-fCFFrameHole2_Y/2.-22*mm;
  transCFFrameHole3_X =0;
  transCFFrameHole3_Y =fCFFrame_Y/2.-fCFFrameHole2_Y/2.-22*mm;
  transCFFrameHole4_X =-82.5*mm;
  transCFFrameHole4_Y =fCFFrame_Y/2.-fCFFrameHole2_Y/2.-22*mm;
  transCFFrameHole5_X =82.5*mm;
  transCFFrameHole5_Y =fCFFrame_Y/2.-fCFFrameHole2_Y/2.-248*mm;
  transCFFrameHole6_X =0;
  transCFFrameHole6_Y =fCFFrame_Y/2.-fCFFrameHole2_Y/2.-248*mm;
  transCFFrameHole7_X =-82.5*mm;
  transCFFrameHole7_Y =fCFFrame_Y/2.-fCFFrameHole2_Y/2.-248*mm;

  transPoronFrameHole0_Y =-4*mm;
  transPoronFrameHole0_Z =-0.61*mm;
  transPoronFrameHole1_Y =-(fPoronFrame_Y/2.-fPoronFrameHole1_Y/2.);
  transPoronFrameHole2_X =82.5*mm;
  transPoronFrameHole2_Y =fPoronFrame_Y/2.-fPoronFrameHole2_Y/2.-6.7*mm;
  transPoronFrameHole3_X =0;
  transPoronFrameHole3_Y =fPoronFrame_Y/2.-fPoronFrameHole2_Y/2.-6.7*mm;
  transPoronFrameHole4_X =-82.5*mm;
  transPoronFrameHole4_Y =fPoronFrame_Y/2.-fPoronFrameHole2_Y/2.-6.7*mm;;

  transRingKHole1_Y = -9*mm;
  transRingKHole1_Z = fRingK_Z/2.-fRingKHole1_Z/2.-1*mm;
  transRingKHole2_Y = fRingK_Y/2.-fRingKHole2_Y/2.-12*mm;
  transRingKHole3_Y = fRingK_Y/2.-fRingKHole3_Y/2.-233*mm;
  transRingSHole_Y = fRingS_Y/2.-fRingSHole_Y/2.-3*mm;  

  transLadderBox_Y = -(fLadderBox_Y/2.-fLadderBoxEnd_Y/2.);
  transLadderBox_Z = fLadderBox_Z/2.-fLadderBoxEnd_Z/2.;

  transLadder_Y = -(fLadder_Y/2.-fLadderEnd_Y/2.);
  transLadder_Z = fLadder_Z/2.-fLadderEnd_Z/2.;

  transKaptonS_Y = 0.1*mm;
  transKaptonS_Z = -0.1*mm;

  fPhysiTracker_X = 0;
  fPhysiTracker_Y = -23.55;
  fPhysiTracker_Z = 353.2*mm;
  
  fPhysiPoronFrame_X = 0;
  fPhysiPoronFrame_Y = -(fCFFrame_Y/2.-fPoronFrame_Y/2.);
  fPhysiPoronFrame_Z = -0.15*mm;
  
  fPhysiLadderBox_X = 0;
  fPhysiLadderBox_Y = -(fCFFrame_Y/2.-fLadderBox_Y/2.);
  fPhysiLadderBox_Z = -0.76*mm;

  fPhysiRingK_X = 0;
  fPhysiRingK_Y = 0;
  fPhysiRingK_Z = 0;

  fPhysiKaptonK_X = 0;
  fPhysiKaptonK_Y = fLadder_Y/2.-fKaptonK_Y/2.-12.71*mm;
  fPhysiKaptonK_Z = fLadder_Z/2.-fKaptonK_Z/2.-1.09*mm;

  fPhysiHybridK_X = 0;
  fPhysiHybridK_Y = -(fLadder_Y/2.-fHybrid_Y/2.);
  fPhysiHybridK_Z = 0.06*mm;

  fPhysiSiliconPlateP_X = 0;
  fPhysiSiliconPlateP_Y = fLadder_Y/2.-fSiSens_Y/2.-9.5*mm;
  fPhysiSiliconPlateP_Z = fLadder_Z/2.-fSiSens_Z/2.-1.24*mm;

  fPhysiSiliconSensorP_X = 0;
  fPhysiSiliconSensorP_Y = fSiSens_Y/2.-fSiActiveSens_Y/2.-3*mm;
  fPhysiSiliconSensorP_Z = 0;

  fPhysiSiliconPlateM_X = 0;
  fPhysiSiliconPlateM_Y = fLadder_Y/2.-fSiSens_Y/2.-119.25*mm;
  fPhysiSiliconPlateM_Z = fLadder_Z/2.-fSiSens_Z/2.-1.24*mm;

  fPhysiSiliconSensorM_X = 0;
  fPhysiSiliconSensorM_Y = fSiSens_Y/2.-fSiActiveSens_Y/2.-0.8*mm;
  fPhysiSiliconSensorM_Z = 0;

  fPhysiKaptonS_X = 0;
  fPhysiKaptonS_Y = fLadder_Y/2.-fKaptonS_Y/2.-226.4*mm;
  fPhysiKaptonS_Z = fLadder_Z/2.-fKaptonS_Z/2.-1.61*mm;

  fPhysiHybridS_X = 0;
  fPhysiHybridS_Y = -(fLadder_Y/2.-fHybrid_Y/2.);
  fPhysiHybridS_Z = -2.76*mm;

  fPhysiRingS_X = 0;
  fPhysiRingS_Y = fLadder_Y/2.-fRingS_Y/2.-9*mm;
  fPhysiRingS_Z = fLadder_Z/2.-fRingS_Z/2.-1.8*mm;
  
  fPhysiHeatSink_X = 0;
  fPhysiHeatSink_Y = -(fLadder_Y/2.-fHeatSink_Y/2.);
  fPhysiHeatSink_Z = -1.35*mm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerConstructionConfig4::Builder(G4VPhysicalVolume* motherVolume)
{
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String tkSDname = "/hepd/silicon";
  TrackerSD * tkSD = new TrackerSD(tkSDname);
  SDman->AddNewDetector(tkSD);

  pMaterial->DefineMaterials();
  G4Material* vacuum    = pMaterial->GetMaterial("Galactic");
  G4Material* poronMat  = pMaterial->GetMaterial(poronMaterial);
  G4Material* cfMat     = pMaterial->GetMaterial(cfiberMaterial);
  G4Material* kapMat    = pMaterial->GetMaterial(kaptonMaterial);
  G4Material* silMat    = pMaterial->GetMaterial(siliconMaterial);
  G4Material* fr4Mat    = pMaterial->GetMaterial("FR4");
  G4Material* copperMat = pMaterial->GetMaterial("Copper");

  G4RotationMatrix* myRot = new G4RotationMatrix; 

  fSolidTracker = new G4Box("fSolidTracker",fTracker_X/2.,fTracker_Y/2.,fTracker_Z/2.);
  fSolidLayer  = new G4Box("fSolidLayer",fLayer_X/2.,fLayer_Y/2.,fLayer_Z/2.);

  G4ThreeVector transCFFrameHole0(0,transCFFrameHole0_Y,transCFFrameHole0_Z);
  G4ThreeVector transCFFrameHole1(0,transCFFrameHole1_Y,transCFFrameHole1_Z);
  G4ThreeVector transCFFrameHole2(transCFFrameHole2_X,transCFFrameHole2_Y,0);
  G4ThreeVector transCFFrameHole3(transCFFrameHole3_X,transCFFrameHole3_Y,0);
  G4ThreeVector transCFFrameHole4(transCFFrameHole4_X,transCFFrameHole4_Y,0);
  G4ThreeVector transCFFrameHole5(transCFFrameHole5_X,transCFFrameHole5_Y,0);
  G4ThreeVector transCFFrameHole6(transCFFrameHole6_X,transCFFrameHole6_Y,0);
  G4ThreeVector transCFFrameHole7(transCFFrameHole7_X,transCFFrameHole7_Y,0);

  fSolidCFFrameHole = new G4Box("fSolidCFFrameHole",fCFFrameHole2_X/2.,fCFFrameHole2_Y/2.,fCFFrameHole2_Z/2.);

  fSolidCFFrame = new G4SubtractionSolid("fSolidCFFrame",
					 new G4SubtractionSolid("CFFrame7",
								new G4SubtractionSolid("CFFrame6",
										       new G4SubtractionSolid("CFFrame5",
													      new G4SubtractionSolid("CFFrame4",
																     new G4SubtractionSolid("CFFrame3",
																			    new G4SubtractionSolid("CFFrame2",
																						   new G4SubtractionSolid("CFFrame1",
																									  new G4Box("CFFrameBox",fCFFrame_X/2.,fCFFrame_Y/2.,fCFFrame_Z/2.),
																									  new G4Box("CFMiddleFrameHole",fCFMiddleFrameHole_X/2.,fCFMiddleFrameHole_Y/2.,fCFMiddleFrameHole_Z/2.),
																									  myRot,transCFFrameHole0),
																						   new G4Box("CFBottomFrameHole",fCFFrameHole1_X/2.,fCFFrameHole1_Y/2.,fCFFrameHole1_Z/2.),
																						   myRot,transCFFrameHole1),
																			    fSolidCFFrameHole,myRot,transCFFrameHole2),
																     fSolidCFFrameHole,myRot,transCFFrameHole3),
													      fSolidCFFrameHole,myRot,transCFFrameHole4),
										       fSolidCFFrameHole,myRot,transCFFrameHole5),
								fSolidCFFrameHole,myRot,transCFFrameHole6),
					 fSolidCFFrameHole,myRot,transCFFrameHole7);

  G4ThreeVector transPoronFrameHole0(0,transPoronFrameHole0_Y,transPoronFrameHole0_Z);
  G4ThreeVector transPoronFrameHole1(0,transPoronFrameHole1_Y,0);
  G4ThreeVector transPoronFrameHole2(transPoronFrameHole2_X,transPoronFrameHole2_Y,0);
  G4ThreeVector transPoronFrameHole3(transPoronFrameHole3_X,transPoronFrameHole3_Y,0);
  G4ThreeVector transPoronFrameHole4(transPoronFrameHole4_X,transPoronFrameHole4_Y,0);
  fSolidPoronFrameHole = new G4Box("fSolidPoronFrameHole",fPoronFrameHole2_X/2.,fPoronFrameHole2_Y/2.,fPoronFrameHole2_Z/2.);

  fSolidPoronFrame = new G4SubtractionSolid("fSolidPoronFrame",
					    new G4SubtractionSolid("PoronFrame4",
								   new G4SubtractionSolid("PoronFrame3",
											  new G4SubtractionSolid("PoronFrame2",
														 new G4SubtractionSolid("PoronFrame1",
																	new G4Box("PoronFrameBox",fPoronFrame_X/2.,fPoronFrame_Y/2.,fPoronFrame_Z/2.),
																	new G4Box("PoronMiddleFrameHole",fPoronMiddleFrameHole_X/2.,fPoronMiddleFrameHole_Y/2.,fPoronMiddleFrameHole_Z/2.),
																	myRot,transPoronFrameHole0),
														 new G4Box("PoronBottomFrameHole",fPoronFrameHole1_X/2.,fPoronFrameHole1_Y/2.,fPoronFrameHole1_Z/2.),
														 myRot,transPoronFrameHole1),
											  fSolidPoronFrameHole,myRot,transPoronFrameHole2),
								   fSolidPoronFrameHole,myRot,transPoronFrameHole3),
					    fSolidPoronFrameHole,myRot,transPoronFrameHole4);
  
  
  
  fSolidSiliconPlate   = new G4Box("fSolidSiliconPlate",fSiSens_X/2.,fSiSens_Y/2.,fSiSens_Z/2.);
  fSolidSiliconSensor  = new G4Box("fSolidSiliconSensor",fSiActiveSens_X/2.,fSiActiveSens_Y/2.,fSiActiveSens_Z/2.);

  G4ThreeVector transLadderBox(0,transLadderBox_Y,transLadderBox_Z);
  fSolidLadderBox      = new G4UnionSolid("fSolidLadderBox",
					  new G4Box("LadderBox",fLadderBox_X/2.,fLadderBox_Y/2.,fLadderBox_Z/2.),
					  new G4Box("LadderBoxEnd",fLadderBoxEnd_X/2.,fLadderBoxEnd_Y/2.,fLadderBoxEnd_Z/2.),
					  myRot,transLadderBox);

  G4ThreeVector transLadder(0,transLadder_Y,transLadder_Z);
  fSolidLadder        = new G4UnionSolid("fSolidLadder",
					 new G4Box("Ladder",fLadder_X/2.,fLadder_Y/2.,fLadder_Z/2.),
					 new G4Box("LadderEnd",fLadderEnd_X/2.,fLadderEnd_Y/2.,fLadderEnd_Z/2.),
					 myRot,transLadder);
  
  G4ThreeVector transKaptonS(0,transKaptonS_Y,transKaptonS_Z);
  fSolidKaptonS        = new G4SubtractionSolid("fSolidKaptonS",
						new G4Box("fKaptonS1",fKaptonS_X/2.,fKaptonS_Y/2.,fKaptonS_Z/2.),
						new G4Box("fKaptonS2",fKaptonS_X/2.+3*mm,fKaptonS_Y/2.,fKaptonS_Z/2.),
						myRot,transKaptonS);
  
  fSolidKaptonK        = new G4Box("fSolidKaptonK",fKaptonK_X/2.,fKaptonK_Y/2.,fKaptonK_Z/2.);
  fSolidHybrid         = new G4Box("fSolidHybrid",fHybrid_X/2.,fHybrid_Y/2.,fHybrid_Z/2.);
  fSolidHeatSink       = new G4Box("fSolidHeatSink",fHeatSink_X/2.,fHeatSink_Y/2.,fHeatSink_Z/2.);
  
  G4ThreeVector transRingKHole1(0,transRingKHole1_Y,transRingKHole1_Z);
  G4ThreeVector transRingKHole2(0,transRingKHole2_Y,0);
  G4ThreeVector transRingKHole3(0,transRingKHole3_Y,0);
  
  fSolidRingK = new G4SubtractionSolid("fSolidRingK",  
				       new G4SubtractionSolid("RingKStep2",
							      new G4SubtractionSolid("RingKStep1",
										     new G4Box("RingKBlock",fRingK_X/2.,fRingK_Y/2.,fRingK_Z/2.),
										     new G4Box("RingKHole1",fRingKHole1_X/2.,fRingKHole1_Y/2.,fRingKHole1_Z/2.),
										     myRot,transRingKHole1),
							      new G4Box("RingKHole2",fRingKHole2_X/2.,fRingKHole2_Y/2.,fRingKHole2_Z/2.),
							      myRot,transRingKHole2),
				       new G4Box("RingKHole3",fRingKHole3_X/2.,fRingKHole3_Y/2.,fRingKHole3_Z/2.),
				       myRot,transRingKHole3);
  
  G4ThreeVector transRingSHole(0,transRingSHole_Y,0);
  fSolidRingS = new G4SubtractionSolid("fSolidRingS",
				       new G4Box("RingSBlock",fRingS_X/2.,fRingS_Y/2.,fRingS_Z/2.),
				       new G4Box("RingSHole",fRingSHole_X/2.,fRingSHole_Y/2.,fRingSHole_Z/2.),
				       myRot,transRingSHole);
  
  fLogicTracker       = new G4LogicalVolume(fSolidTracker,vacuum,"fLogicTracker");
  fLogicLayer         = new G4LogicalVolume(fSolidLayer,vacuum,"fLogicLayer");
  fLogicCFFrame       = new G4LogicalVolume(fSolidCFFrame,cfMat,"fLogicCFFrame");
  fLogicPoronFrame    = new G4LogicalVolume(fSolidPoronFrame,poronMat,"fLogicPoronFrame");
  fLogicSiliconPlateP = new G4LogicalVolume(fSolidSiliconPlate,silMat,"fLogicSiliconPlateP");
  fLogicSiliconPlateM = new G4LogicalVolume(fSolidSiliconPlate,silMat,"fLogicSiliconPlateM");
  fLogicSiliconSensor = new G4LogicalVolume(fSolidSiliconSensor,silMat,"fLogicSiliconSensor");
  fLogicLadderBox     = new G4LogicalVolume(fSolidLadderBox,vacuum,"fLogicLadderBox");
  fLogicLadder        = new G4LogicalVolume(fSolidLadder,vacuum,"fLogicLadder");
  fLogicKaptonS       = new G4LogicalVolume(fSolidKaptonS,kapMat,"fLogicKaptonS");
  fLogicKaptonK       = new G4LogicalVolume(fSolidKaptonK,kapMat,"fLogicKaptonK");
  fLogicHybrid        = new G4LogicalVolume(fSolidHybrid,fr4Mat,"fLogicHybrid");
  fLogicHeatSink      = new G4LogicalVolume(fSolidHeatSink,copperMat,"fLogicHeatSink");
  fLogicRingK         = new G4LogicalVolume(fSolidRingK,cfMat,"fLogicRingK");
  fLogicRingS         = new G4LogicalVolume(fSolidRingS,cfMat,"fLogicRingS");
  
  fLogicSiliconSensor->SetSensitiveDetector(tkSD);


  G4RotationMatrix * rotTracker = new G4RotationMatrix;
  rotTracker->rotateZ(180*deg);
  rotTracker->rotateX(180*deg);

  fPhysiTracker = new G4PVPlacement(rotTracker,
				    G4ThreeVector(fPhysiTracker_X,fPhysiTracker_Y,fPhysiTracker_Z),
				    "Tracker",
				    fLogicTracker,
				    motherVolume,
				    false,0,true);

  fPhysiLayer = new G4PVReplica("Layer",
				fLogicLayer,
				fPhysiTracker,
				kZAxis,
				fLayerNumber,
				fLayer_Z);

  fPhysiCFFrame = new G4PVPlacement(0,
				    G4ThreeVector(0,0,0),
				    "CFFrame",
				    fLogicCFFrame,
				    fPhysiLayer,
				    false,0,true);


  fPhysiPoronFrame = new G4PVPlacement(0,
				       G4ThreeVector(fPhysiPoronFrame_X,fPhysiPoronFrame_Y,fPhysiPoronFrame_Z),
				       "PoronFrame",
				       fLogicPoronFrame,
				       fPhysiLayer,
				       false,0,true);

  fPhysiLadderBox = new G4PVPlacement(0,
				      G4ThreeVector(fPhysiLadderBox_X,fPhysiLadderBox_Y,fPhysiLadderBox_Z),
				      "LadderBox",
				      fLogicLadderBox,
				      fPhysiLayer,
				      false,0,true);
  
  fPhysiLadder = new G4PVReplica("Ladder",
				 fLogicLadder,
				 fPhysiLadderBox,
				 kXAxis,
				 fLadderNumber,
				 fLadder_X);

  fPhysiRingK = new G4PVPlacement(0,
				  G4ThreeVector(fPhysiRingK_X,fPhysiRingK_Y,fPhysiRingK_Z),
				  "RingK",
				  fLogicRingK,
				  fPhysiLadder,
				  false,0,true);				  

  fPhysiKaptonK = new G4PVPlacement(0,
				    G4ThreeVector(fPhysiKaptonK_X,fPhysiKaptonK_Y,fPhysiKaptonK_Z),
				    "KaptonK",
				    fLogicKaptonK,
				    fPhysiLadder,
				    false,0,true);				  
  
  fPhysiHybridK = new G4PVPlacement(0,
				    G4ThreeVector(fPhysiHybridK_X,fPhysiHybridK_Y,fPhysiHybridK_Z),
				    "HybridK",
				    fLogicHybrid,
				    fPhysiLadder,
				    false,0,true);				  
  
  fPhysiSiliconPlateP = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiSiliconPlateP_X,fPhysiSiliconPlateP_Y,fPhysiSiliconPlateP_Z),
					  "SiliconPlateP",
					  fLogicSiliconPlateP,
					  fPhysiLadder,
					  false,0,true);				  

  fPhysiSiliconSensorP = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiSiliconSensorP_X,fPhysiSiliconSensorP_Y,fPhysiSiliconSensorP_Z),
					  "SiliconSensorP",
					  fLogicSiliconSensor,
					  fPhysiSiliconPlateP,
					  false,0,true);				  
  
  fPhysiSiliconPlateM = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiSiliconPlateM_X,fPhysiSiliconPlateM_Y,fPhysiSiliconPlateM_Z),
					  "SiliconPlateM",
					  fLogicSiliconPlateM,
					  fPhysiLadder,
					  false,0,true);				  

  fPhysiSiliconSensorM = new G4PVPlacement(0,
					   G4ThreeVector(fPhysiSiliconSensorM_X,fPhysiSiliconSensorM_Y,fPhysiSiliconSensorM_Z),
					  "SiliconSensorM",
					  fLogicSiliconSensor,
					  fPhysiSiliconPlateM,
					  false,0,true);				  
  
  fPhysiKaptonS = new G4PVPlacement(0,
				    G4ThreeVector(fPhysiKaptonS_X,fPhysiKaptonS_Y,fPhysiKaptonS_Z),
				    "KaptonS",
				    fLogicKaptonS,
				    fPhysiLadder,
				    false,0,true);				  
  
  fPhysiHybridS = new G4PVPlacement(0,
				    G4ThreeVector(fPhysiHybridS_X,fPhysiHybridS_Y,fPhysiHybridS_Z),
				    "HybridS",
				    fLogicHybrid,
				    fPhysiLadder,
  				    false,0,true);	
  

  fPhysiHeatSink = new G4PVPlacement(0,
				    G4ThreeVector(fPhysiHeatSink_X,fPhysiHeatSink_Y,fPhysiHeatSink_Z),
				    "HeatSink",
				    fLogicHeatSink,
				    fPhysiLadder,
				    false,0,true);	
  
  fPhysiRingS = new G4PVPlacement(0,
				  G4ThreeVector(fPhysiRingS_X,fPhysiRingS_Y,fPhysiRingS_Z),
				  "RingS",
				  fLogicRingS,
				  fPhysiLadder,
				  false,0,true);
  			  
  
  //Visualization Attribute
  G4VisAttributes* attBlue = new G4VisAttributes(G4Colour::Blue());
  attBlue->SetVisibility(true);
  attBlue->SetForceAuxEdgeVisible(true);
  fLogicSiliconPlateP->SetVisAttributes(attBlue); 
  fLogicSiliconPlateM->SetVisAttributes(attBlue); 

  G4VisAttributes* attGray = new G4VisAttributes(G4Colour::Gray());
  attGray->SetVisibility(true);
  attGray->SetForceAuxEdgeVisible(true);
  fLogicSiliconSensor->SetVisAttributes(attGray);

  G4VisAttributes* attYellow = new G4VisAttributes(G4Colour::Yellow());
  attYellow->SetVisibility(true);
  attYellow->SetForceAuxEdgeVisible(true);
  fLogicPoronFrame->SetVisAttributes(attYellow);   

  G4VisAttributes* attMagenta = new G4VisAttributes(G4Colour::Magenta());
  attMagenta->SetVisibility(true);
  attMagenta->SetForceAuxEdgeVisible(true);
  fLogicKaptonS->SetVisAttributes(attMagenta);      
  fLogicKaptonK->SetVisAttributes(attMagenta);      

  G4VisAttributes* attBrown = new G4VisAttributes(G4Colour::Brown());
  attBrown->SetVisibility(true);
  attBrown->SetForceAuxEdgeVisible(true);
  fLogicCFFrame->SetVisAttributes(attBrown);
  fLogicRingK->SetVisAttributes(attBrown);        
  fLogicRingS->SetVisAttributes(attBrown);        
      
  G4VisAttributes* attGreen = new G4VisAttributes(G4Colour::Green());
  attGreen->SetVisibility(true);
  attGreen->SetForceAuxEdgeVisible(true);
  fLogicHybrid->SetVisAttributes(attGreen);       

  G4VisAttributes* attRed = new G4VisAttributes(G4Colour::Red());
  attRed->SetVisibility(true);
  attRed->SetForceAuxEdgeVisible(true);
  fLogicHybrid->SetVisAttributes(attRed); 
      
  G4VisAttributes* attBlack = new G4VisAttributes(G4Colour::Black());
  attBlack->SetVisibility(true);
  attBlack->SetForceAuxEdgeVisible(true);
  fLogicHeatSink->SetVisAttributes(attBlack);     
  
  G4VisAttributes * attInvisible = new G4VisAttributes();
  attInvisible->SetVisibility(false);
  attInvisible->SetForceAuxEdgeVisible(false);
  fLogicTracker->SetVisAttributes(attInvisible);
  fLogicLayer->SetVisAttributes(attInvisible);        
  fLogicLadderBox->SetVisAttributes(attInvisible);    
  fLogicLadder->SetVisAttributes(attInvisible);       

}


