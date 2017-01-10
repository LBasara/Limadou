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

#include "TrackerConstructionConfig1.hh"

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

TrackerConstructionConfig1::TrackerConstructionConfig1()
  :fSolidSiSensor(0),fSolidActiveLadderL1(0),fSolidKaptonNLadderL1(0),fSolidKaptonPLadderL1(0),fSolidLadderPlaneL1(0),fSolidLadderL1(0),
   fSolidL1Support(0),fSolidL1SupportHoleU(0),fSolidL1SupportHoleD(0),fSolidL1SupportU(0),fSolidL1SupportD(0),fSolidL1SupportUFinal(0),fSolidL1SupportDFinal(0),fSolidL1(0),
   fSolidActiveLadderL2(0),fSolidKaptonNLadderL2(0),fSolidKaptonPLadderL2(0),fSolidLadderPlaneL2(0),fSolidLadderL2(0),
  fSolidL2Support(0),fSolidL2SupportHoleU(0),fSolidL2SupportHoleD(0),fSolidL2SupportU(0),fSolidL2SupportD(0),fSolidL2SupportUFinal(0),fSolidL2SupportDFinal(0),fSolidL2(0),
  fLogicSiSensor(0),
  fLogicActiveLadderL1(0),fLogicKaptonNLadderL1(0),fLogicKaptonPLadderL1(0),fLogicLadderPlaneL1(0),fLogicLadderL1(0),fLogicL1SupportU(0),fLogicL1SupportD(0),fLogicL1(0),
  fLogicActiveLadderL2(0),fLogicKaptonNLadderL2(0),fLogicKaptonPLadderL2(0),fLogicLadderPlaneL2(0),fLogicLadderL2(0),fLogicL2SupportU(0),fLogicL2SupportD(0),fLogicL2(0),
  fPhysiL1(0),fPhysiL1SupportU(0),fPhysiL1SupportD(0),fPhysiLadderPlaneL1(0),fPhysiLadderL1(0),fPhysiKaptonNLadderL1(0),
  fPhysiKaptonPLadderL1(0),fPhysiActiveLadderL1(0),fPhysiSiSensorL1(0),
  fPhysiL2(0),fPhysiL2SupportU(0),fPhysiL2SupportD(0),fPhysiLadderPlaneL2(0),fPhysiLadderL2(0),fPhysiKaptonNLadderL2(0),
  fPhysiKaptonPLadderL2(0),fPhysiActiveLadderL2(0),fPhysiSiSensorL2(0)
{
  pMaterial = new HEPDSWMaterial();

  fSiSensNumberL1 = 6; 
  fSiSensNumberL2 = 9; 

  fLadderNumberL1 = 4;
  fLadderNumberL2 = 5;

  fSiSens_X = 72.045*mm;
  fSiSens_Y = 41.36*mm;
  fSiSens_Z = 0.3*mm;

  fActiveLadderL1_X =fSiSens_X;
  fActiveLadderL1_Y =fSiSens_Y*fSiSensNumberL1;
  fActiveLadderL1_Z =fSiSens_Z;
  
  fKaptonNLadderL1_X = 67.2*mm;
  fKaptonNLadderL1_Y = 273.828*mm;
  fKaptonNLadderL1_Z = 0.3*mm;
  //fKaptonNLadderL1_Z = 0.05*mm;

  fKaptonPLadderL1_X = 72*mm;
  fKaptonPLadderL1_Y = 29*mm;
  fKaptonPLadderL1_Z = 0.3*mm;
  //fKaptonPLadderL1_Z = 0.05*mm;

  fLadderL1_X = fSiSens_X+5.5*mm;
  fLadderL1_Y = fKaptonNLadderL1_Y;
  fLadderL1_Z = fActiveLadderL1_Z+fKaptonNLadderL1_Z+fKaptonPLadderL1_Z;

  fLadderPlaneL1_X = (fSiSens_X+5.5*mm)*fLadderNumberL1;
  fLadderPlaneL1_Y = fLadderL1_Y;
  fLadderPlaneL1_Z = fLadderL1_Z;

  fL1Support_X = 314.68*mm;
  fL1Support_Y = 258.36*mm;
  fL1Support_Z = 3*mm;

  fL1SupportHoleU_X = 66.045*mm;
  fL1SupportHoleU_Y = 238.36*mm; 
  fL1SupportHoleU_Z = 4*mm;

  fL1SupportHoleD_X = 66.045*mm;
  fL1SupportHoleD_Y = 38.4*mm; 
  fL1SupportHoleD_Z = 4*mm;

  fL1_X =314.68*mm;;
  fL1_Y =278.828*mm;
  fL1_Z =8.6*mm;

  fActiveLadderL2_X =fSiSens_X;
  fActiveLadderL2_Y =fSiSens_Y*fSiSensNumberL2;
  fActiveLadderL2_Z =fSiSens_Z;

  fKaptonNLadderL2_X = 67.2*mm;
  fKaptonNLadderL2_Y = 398.028*mm;
  fKaptonNLadderL2_Z = 0.3*mm;
  //fKaptonNLadderL2_Z = 0.05*mm;

  fKaptonPLadderL2_X = 72*mm;
  fKaptonPLadderL2_Y = 29*mm;
  fKaptonPLadderL2_Z = 0.3*mm;
  //fKaptonPLadderL2_Z = 0.05*mm;

  fLadderL2_X = fSiSens_X+5.5*mm;
  fLadderL2_Y = fKaptonNLadderL2_Y;
  fLadderL2_Z = fActiveLadderL2_Z+fKaptonNLadderL2_Z+fKaptonPLadderL2_Z;

  fLadderPlaneL2_X = (fSiSens_X+5.5*mm)*fLadderNumberL2;
  fLadderPlaneL2_Y = fLadderL2_Y;
  fLadderPlaneL2_Z = fLadderL2_Z;

  fL2Support_X = 392.225*mm;
  fL2Support_Y = 382.58*mm;
  fL2Support_Z = 3*mm;

  fL2SupportHoleU_X = 66.045*mm;
  fL2SupportHoleU_Y = 362.58*mm; 
  fL2SupportHoleU_Z = 4*mm;

  fL2SupportHoleD_X = 66.045*mm;
  fL2SupportHoleD_Y = 38.4*mm; 
  fL2SupportHoleD_Z = 4*mm;

  fL2_X = 392.225*mm;
  fL2_Y = 403.028*mm;
  fL2_Z = 8.6*mm;

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

TrackerConstructionConfig1::~TrackerConstructionConfig1()
{
  if (pMaterial) delete pMaterial;
}

void TrackerConstructionConfig1::ComputeObjectsPositioning(){

  fL1SupportHoleUStep_X = fL1SupportHoleU_X+11.5*mm;
  fL1SupportHoleU_Xmin  = fL1Support_X/2.-8*mm-fL1SupportHoleU_X/2.;

  fL1SupportHoleDStep_X = fL1SupportHoleD_X+11.5*mm;
  fL1SupportHoleD_Xmin  = fL1Support_X/2.-8*mm-fL1SupportHoleD_X/2.;

  fL1SupportHoleDStep_Y = fL1SupportHoleD_Y+3*mm;
  fL1SupportHoleD_Ymin  = fL1Support_Y/2.-8*mm-fL1SupportHoleD_Y/2.;

  fPhysiL1_X = 0;
  fPhysiL1_Y = 0;
  fPhysiL1_Z = 151.8*mm+18.36*mm+fL1_Z/2.;
  
  fPhysiL1SupportU_X = 0;
  fPhysiL1SupportU_Y = -(fL1_Y/2.-fL1Support_Y/2.);
  fPhysiL1SupportU_Z = fL1_Z/2-fL1Support_Z/2.;

  fPhysiL1SupportD_X = 0;
  fPhysiL1SupportD_Y = -(fL1_Y/2.-fL1Support_Y/2.);
  fPhysiL1SupportD_Z = -fL1_Z/2+fL1Support_Z/2.;

  fPhysiLadderPlaneL1_X = 0;
  fPhysiLadderPlaneL1_Y = fL1_Y/2.-fLadderL1_Y/2.;
  fPhysiLadderPlaneL1_Z = 0;

  fPhysiKaptonNLadderL1_X = 0;
  fPhysiKaptonNLadderL1_Y = 0;
  fPhysiKaptonNLadderL1_Z = -(fLadderL1_Z/2.-fKaptonNLadderL1_Z/2.);

  fPhysiKaptonPLadderL1_X = 0;
  fPhysiKaptonPLadderL1_Y = fLadderL1_Y/2.-fKaptonPLadderL1_Y/2.;
  fPhysiKaptonPLadderL1_Z = fLadderL1_Z/2.-fKaptonPLadderL1_Z/2.;

  fPhysiActiveLadderL1_X = 0;
  fPhysiActiveLadderL1_Y = -(fLadderL1_Y/2.-fActiveLadderL1_Y/2.);
  fPhysiActiveLadderL1_Z = 0;
  
  fL2SupportHoleUStep_X = fL2SupportHoleU_X+11.5*mm;
  fL2SupportHoleU_Xmin  = fL2Support_X/2.-8*mm-fL2SupportHoleU_X/2.;

  fL2SupportHoleDStep_X = fL2SupportHoleD_X+11.5*mm;
  fL2SupportHoleD_Xmin  = fL2Support_X/2.-8*mm-fL2SupportHoleD_X/2.;

  fL2SupportHoleDStep_Y = fL2SupportHoleD_Y+3*mm;
  fL2SupportHoleD_Ymin  = fL2Support_Y/2.-8*mm-fL2SupportHoleD_Y/2.;

  fPhysiL2_X = 0;
  fPhysiL2_Y = 0;
  fPhysiL2_Z = 151.8*mm+38.36*mm+fL2_Z/2.;
  
  fPhysiL2SupportU_X = 0;
  fPhysiL2SupportU_Y = -(fL2_Y/2.-fL2Support_Y/2.);
  fPhysiL2SupportU_Z = fL2_Z/2-fL2Support_Z/2.;

  fPhysiL2SupportD_X = 0;
  fPhysiL2SupportD_Y = -(fL2_Y/2.-fL2Support_Y/2.);
  fPhysiL2SupportD_Z = -fL2_Z/2+fL2Support_Z/2.;

  fPhysiLadderPlaneL2_X = 0;
  fPhysiLadderPlaneL2_Y = fL2_Y/2.-fLadderL2_Y/2.;
  fPhysiLadderPlaneL2_Z = 0;

  fPhysiKaptonNLadderL2_X = 0;
  fPhysiKaptonNLadderL2_Y = 0;
  fPhysiKaptonNLadderL2_Z = -(fLadderL2_Z/2.-fKaptonNLadderL2_Z/2.);

  fPhysiKaptonPLadderL2_X = 0;
  fPhysiKaptonPLadderL2_Y = fLadderL2_Y/2.-fKaptonPLadderL2_Y/2.;
  fPhysiKaptonPLadderL2_Z = fLadderL2_Z/2.-fKaptonPLadderL2_Z/2.;

  fPhysiActiveLadderL2_X = 0;
  fPhysiActiveLadderL2_Y = -(fLadderL2_Y/2.-fActiveLadderL2_Y/2.);
  fPhysiActiveLadderL2_Z = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerConstructionConfig1::Builder(G4VPhysicalVolume* motherVolume)
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
  fSolidActiveLadderL1 = new G4Box("ActiveLadder",fActiveLadderL1_X/2.,fActiveLadderL1_Y/2.,fActiveLadderL1_Z/2.);
  fSolidKaptonNLadderL1 = new G4Box("KaptonN",fKaptonNLadderL1_X/2.,fKaptonNLadderL1_Y/2.,fKaptonNLadderL1_Z/2.);
  fSolidKaptonPLadderL1 = new G4Box("KaptonP",fKaptonPLadderL1_X/2.,fKaptonPLadderL1_Y/2.,fKaptonPLadderL1_Z/2.);
  fSolidLadderPlaneL1 = new G4Box("LadderPlane",fLadderPlaneL1_X/2.,fLadderPlaneL1_Y/2.,fLadderPlaneL1_Z/2.);
  fSolidLadderL1 = new G4Box("Ladder",fLadderL1_X/2.,fLadderL1_Y/2.,fLadderL1_Z/2.);
  fSolidL1Support = new G4Box("Support",fL1Support_X/2.,fL1Support_Y/2.,fL1Support_Z/2.);
  fSolidL1SupportHoleU = new G4Box("HoleUp",fL1SupportHoleU_X/2.,fL1SupportHoleU_Y/2.,fL1SupportHoleU_Z/2.);
  fSolidL1SupportHoleD = new G4Box("HoleDown",fL1SupportHoleD_X/2.,fL1SupportHoleD_Y/2.,fL1SupportHoleD_Z/2.);
  for(int i=0;i<fLadderNumberL1;i++){
    if(i==0){
      fSolidL1SupportU.push_back(new G4SubtractionSolid("L1SupportUp",
							fSolidL1Support,
							fSolidL1SupportHoleU,
							0,
							G4ThreeVector(-fL1SupportHoleU_Xmin+i*fL1SupportHoleUStep_X,0,0)));
    }else{
      fSolidL1SupportU.push_back(new G4SubtractionSolid("L1SupportUp",
							fSolidL1SupportU[i-1],
							fSolidL1SupportHoleU,
							0,
							G4ThreeVector(-fL1SupportHoleU_Xmin+i*fL1SupportHoleUStep_X,0,0)));
    }
  }
  for(int i=0;i<fLadderNumberL1;i++){
    for(int j=0;j<fSiSensNumberL1;j++){
      if(j==0&&i==0){
	fSolidL1SupportD.push_back(new G4SubtractionSolid("L1SupportDown",
							  fSolidL1Support,
							  fSolidL1SupportHoleD,
							  0,
							  G4ThreeVector(-fL1SupportHoleD_Xmin+i*fL1SupportHoleDStep_X,
									-fL1SupportHoleD_Ymin+j*fL1SupportHoleDStep_Y,0)));
      }else{
	fSolidL1SupportD.push_back(new G4SubtractionSolid("L1SupportDown",
							  fSolidL1SupportD[(j+(fSiSensNumberL1*i))-1],
							  fSolidL1SupportHoleD,
							  0,
							  G4ThreeVector(-fL1SupportHoleD_Xmin+i*fL1SupportHoleDStep_X,
									-fL1SupportHoleD_Ymin+j*fL1SupportHoleDStep_Y,0)));
	
      }
    }
  }

  fSolidL1SupportUFinal = fSolidL1SupportU.back();
  fSolidL1SupportDFinal = fSolidL1SupportD.back();
  
  fSolidL1SupportU.clear();
  fSolidL1SupportU.resize(0);

  fSolidL1SupportD.clear();
  fSolidL1SupportD.resize(0);

  fSolidL1 = new G4Box("Layer1",fL1_X/2.,fL1_Y/2.,fL1_Z/2.);
  

  fSolidActiveLadderL2 = new G4Box("ActiveLadder",fActiveLadderL2_X/2.,fActiveLadderL2_Y/2.,fActiveLadderL2_Z/2.);
  fSolidKaptonNLadderL2 = new G4Box("KaptonN",fKaptonNLadderL2_X/2.,fKaptonNLadderL2_Y/2.,fKaptonNLadderL2_Z/2.);
  fSolidKaptonPLadderL2 = new G4Box("KaptonP",fKaptonPLadderL2_X/2.,fKaptonPLadderL2_Y/2.,fKaptonPLadderL2_Z/2.);
  fSolidLadderPlaneL2 = new G4Box("LadderPlane",fLadderPlaneL2_X/2.,fLadderPlaneL2_Y/2.,fLadderPlaneL2_Z/2.);
  fSolidLadderL2 = new G4Box("Ladder",fLadderL2_X/2.,fLadderL2_Y/2.,fLadderL2_Z/2.);
  fSolidL2Support = new G4Box("Support",fL2Support_X/2.,fL2Support_Y/2.,fL2Support_Z/2.);
  fSolidL2SupportHoleU = new G4Box("HoleUp",fL2SupportHoleU_X/2.,fL2SupportHoleU_Y/2.,fL2SupportHoleU_Z/2.);
  fSolidL2SupportHoleD = new G4Box("HoleDown",fL2SupportHoleD_X/2.,fL2SupportHoleD_Y/2.,fL2SupportHoleD_Z/2.);
  for(int i=0;i<fLadderNumberL2;i++){
    if(i==0){
      fSolidL2SupportU.push_back(new G4SubtractionSolid("L2SupportUp",
							fSolidL2Support,
							fSolidL2SupportHoleU,
							0,
							G4ThreeVector(-fL2SupportHoleU_Xmin+i*fL2SupportHoleUStep_X,0,0)));
    }else{
      fSolidL2SupportU.push_back(new G4SubtractionSolid("L2SupportUp",
							fSolidL2SupportU[i-1],
							fSolidL2SupportHoleU,
							0,
							G4ThreeVector(-fL2SupportHoleU_Xmin+i*fL2SupportHoleUStep_X,0,0)));
    }
  }
  for(int i=0;i<fLadderNumberL2;i++){
    for(int j=0;j<fSiSensNumberL2;j++){
      if(i==0&&j==0){
	fSolidL2SupportD.push_back(new G4SubtractionSolid("L2SupportDown",
							  fSolidL2Support,
							  fSolidL2SupportHoleD,
							  0,
							  G4ThreeVector(-fL2SupportHoleD_Xmin+i*fL2SupportHoleDStep_X,
									-fL2SupportHoleD_Ymin+j*fL2SupportHoleDStep_Y,0)));
      }else{
	fSolidL2SupportD.push_back(new G4SubtractionSolid("L2SupportDown",
							  fSolidL2SupportD[(j+(fSiSensNumberL2*i))-1],
							  fSolidL2SupportHoleD,
							  0,
							  G4ThreeVector(-fL2SupportHoleD_Xmin+i*fL2SupportHoleDStep_X,
									-fL2SupportHoleD_Ymin+j*fL2SupportHoleDStep_Y,0)));

      }
    }
  }

  fSolidL2SupportUFinal = fSolidL2SupportU.back();
  fSolidL2SupportDFinal = fSolidL2SupportD.back();
  
  fSolidL2SupportU.clear();
  fSolidL2SupportU.resize(0);

  fSolidL2SupportD.clear();
  fSolidL2SupportD.resize(0);

  fSolidL2 = new G4Box("Layer2",fL2_X/2.,fL2_Y/2.,fL2_Z/2.);
  
  fLogicSiSensor = new G4LogicalVolume(fSolidSiSensor,silMat,"SiSensor");
  fLogicActiveLadderL1 = new G4LogicalVolume(fSolidActiveLadderL1,vacuum,"ActiveLadder");
  fLogicKaptonNLadderL1= new G4LogicalVolume(fSolidKaptonNLadderL1,kapMat,"KaptonN");
  fLogicKaptonPLadderL1= new G4LogicalVolume(fSolidKaptonPLadderL1,kapMat,"KaptonN");
  fLogicLadderPlaneL1 = new G4LogicalVolume(fSolidLadderPlaneL1,vacuum,"LadderPlane");
  fLogicLadderL1 = new G4LogicalVolume(fSolidLadderL1,vacuum,"Ladder");
  fLogicL1SupportU = new G4LogicalVolume(fSolidL1SupportUFinal,cfMat,"L1Support");
  fLogicL1SupportD = new G4LogicalVolume(fSolidL1SupportDFinal,cfMat,"L1Support");
  fLogicL1 = new G4LogicalVolume(fSolidL1,vacuum,"Layer1");

  fLogicActiveLadderL2 = new G4LogicalVolume(fSolidActiveLadderL2,vacuum,"ActiveLadder");
  fLogicKaptonNLadderL2= new G4LogicalVolume(fSolidKaptonNLadderL2,kapMat,"KaptonN");
  fLogicKaptonPLadderL2= new G4LogicalVolume(fSolidKaptonPLadderL2,kapMat,"KaptonN");
  fLogicLadderPlaneL2 = new G4LogicalVolume(fSolidLadderPlaneL2,vacuum,"LadderPLane");
  fLogicLadderL2 = new G4LogicalVolume(fSolidLadderL2,vacuum,"Ladder");
  fLogicL2SupportU = new G4LogicalVolume(fSolidL2SupportUFinal,cfMat,"L2Support");
  fLogicL2SupportD = new G4LogicalVolume(fSolidL2SupportDFinal,cfMat,"L2Support");
  fLogicL2 = new G4LogicalVolume(fSolidL2,vacuum,"Layer1");

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
				       fLogicL1SupportU,
				       fPhysiL1,
				       false,
				       0,true);
  
  fPhysiL1SupportD = new G4PVPlacement(0,
				       G4ThreeVector(fPhysiL1SupportD_X,fPhysiL1SupportD_Y,fPhysiL1SupportD_Z),
				       "L1Support",
				       fLogicL1SupportD,
				       fPhysiL1,
				       false,
				       0,true);
  
  fPhysiLadderPlaneL1 = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiLadderPlaneL1_X,fPhysiLadderPlaneL1_Y,fPhysiLadderPlaneL1_Z),
					  "LadderPLane",
					  fLogicLadderPlaneL1,
					  fPhysiL1,
					  false,
					  0,true);

  fPhysiLadderL1 = new G4PVReplica("Ladder",
				   fLogicLadderL1,
				   fPhysiLadderPlaneL1,
				   kXAxis,
				   fLadderNumberL1,
				   fLadderL1_X);
  
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
				     kYAxis,
				     fSiSensNumberL1,
				     fSiSens_Y);
    

  
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
				       fLogicL2SupportU,
				       fPhysiL2,
				       false,
				       0,true);
    
  fPhysiL2SupportD = new G4PVPlacement(0,
				       G4ThreeVector(fPhysiL2SupportD_X,fPhysiL2SupportD_Y,fPhysiL2SupportD_Z),
				       "L2Support",
				       fLogicL2SupportD,
				       fPhysiL2,
				       false,
				       0,true);
    
  fPhysiLadderPlaneL2 = new G4PVPlacement(0,
					  G4ThreeVector(fPhysiLadderPlaneL2_X,fPhysiLadderPlaneL2_Y,fPhysiLadderPlaneL2_Z),
					  "LadderPLane",
					  fLogicLadderPlaneL2,
					  fPhysiL2,
					  false,
					  0,true);
    
  fPhysiLadderL2 = new G4PVReplica("Ladder",
				   fLogicLadderL2,
				   fPhysiLadderPlaneL2,
				   kXAxis,
				   fLadderNumberL2,
				   fLadderL2_X);
    
    
    
    
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
				     kYAxis,
				     fSiSensNumberL2,
				     fSiSens_Y);
 
    
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
  fLogicL1SupportU->SetVisAttributes(attBrown);
  fLogicL2SupportU->SetVisAttributes(attBrown);
  fLogicL1SupportD->SetVisAttributes(attBrown);
  fLogicL2SupportD->SetVisAttributes(attBrown);

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
  fLogicLadderPlaneL1->SetVisAttributes(attInvisible);
  fLogicLadderL2->SetVisAttributes(attInvisible);
  fLogicLadderPlaneL2->SetVisAttributes(attInvisible);
  fLogicL1->SetVisAttributes(attInvisible);
  fLogicL2->SetVisAttributes(attInvisible);

}


