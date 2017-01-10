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
/// \file electromagnetic/TestEm3/include/TrackerConstructionConfig4.hh
/// \brief Definition of the TrackerConstructionConfig4 class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef TrackerConstructionConfig4_h
#define TrackerConstructionConfig4_h 1

#include "globals.hh"
#include <vector>

class G4Box;
class G4SubtractionSolid;
class G4UnionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class HEPDSWMaterial;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackerConstructionConfig4
{
public:
  
   TrackerConstructionConfig4();
  ~TrackerConstructionConfig4();

  inline void SetSiliconMaterial(G4String aMat){siliconMaterial=aMat;}
  inline void SetKaptonMaterial(G4String aMat){kaptonMaterial=aMat;}
  inline void SetCarbonFiberMaterial(G4String aMat){cfiberMaterial=aMat;}
  inline void SetPoronMaterial(G4String aMat){poronMaterial=aMat;}

  void Builder(G4VPhysicalVolume* motherVolume);
  
private:

  void ComputeObjectsPositioning();

  G4String cfiberMaterial;
  G4String kaptonMaterial;
  G4String siliconMaterial;
  G4String poronMaterial;

  G4int fLayerNumber;
  G4int fLadderNumber;

  G4double  fTracker_X;
  G4double  fTracker_Y;
  G4double  fTracker_Z;

  G4double  fLayer_X;
  G4double  fLayer_Y;
  G4double  fLayer_Z;

  G4double  fCFFrame_X;
  G4double  fCFFrame_Y;
  G4double  fCFFrame_Z;
  
  G4double  fCFMiddleFrameHole_X;
  G4double  fCFMiddleFrameHole_Y;
  G4double  fCFMiddleFrameHole_Z;

  G4double  fCFFrameHole1_X;
  G4double  fCFFrameHole1_Y;
  G4double  fCFFrameHole1_Z;

  G4double  fCFFrameHole2_X;
  G4double  fCFFrameHole2_Y;
  G4double  fCFFrameHole2_Z;

  G4double  fPoronFrame_X;
  G4double  fPoronFrame_Y;
  G4double  fPoronFrame_Z;
	      
  G4double  fPoronMiddleFrameHole_X;
  G4double  fPoronMiddleFrameHole_Y;
  G4double  fPoronMiddleFrameHole_Z;

  G4double  fPoronFrameHole1_X;
  G4double  fPoronFrameHole1_Y;
  G4double  fPoronFrameHole1_Z;

  G4double  fPoronFrameHole2_X;
  G4double  fPoronFrameHole2_Y;
  G4double  fPoronFrameHole2_Z;

  G4double  fSiSens_X;
  G4double  fSiSens_Y;
  G4double  fSiSens_Z;

  G4double  fSiActiveSens_X;
  G4double  fSiActiveSens_Y;
  G4double  fSiActiveSens_Z;

  G4double  fLadderBox_X;
  G4double  fLadderBox_Y;
  G4double  fLadderBox_Z;

  G4double  fLadderBoxEnd_X;
  G4double  fLadderBoxEnd_Y;
  G4double  fLadderBoxEnd_Z;

  G4double  fLadder_X;
  G4double  fLadder_Y;
  G4double  fLadder_Z;
  
  G4double  fLadderEnd_X;
  G4double  fLadderEnd_Y;
  G4double  fLadderEnd_Z;
  
  G4double  fKaptonS_X;
  G4double  fKaptonS_Y;
  G4double  fKaptonS_Z;

  G4double  fKaptonK_X;
  G4double  fKaptonK_Y;
  G4double  fKaptonK_Z;

  G4double  fHybrid_X;
  G4double  fHybrid_Y;
  G4double  fHybrid_Z;

  G4double  fHeatSink_X;
  G4double  fHeatSink_Y;
  G4double  fHeatSink_Z;

  G4double  fRingK_X;
  G4double  fRingK_Y;
  G4double  fRingK_Z;

  G4double  fRingKHole1_X;
  G4double  fRingKHole1_Y;
  G4double  fRingKHole1_Z;

  G4double  fRingKHole2_X;
  G4double  fRingKHole2_Y;
  G4double  fRingKHole2_Z;

  G4double  fRingKHole3_X;
  G4double  fRingKHole3_Y;
  G4double  fRingKHole3_Z;

  G4double  fRingS_X;
  G4double  fRingS_Y;
  G4double  fRingS_Z;

  G4double  fRingSHole_X;
  G4double  fRingSHole_Y;
  G4double  fRingSHole_Z;

  G4double  transCFFrameHole0_Y;
  G4double  transCFFrameHole0_Z;
  G4double  transCFFrameHole1_Y;
  G4double  transCFFrameHole1_Z;
  G4double  transCFFrameHole2_X;
  G4double  transCFFrameHole2_Y;
  G4double  transCFFrameHole3_X;
  G4double  transCFFrameHole3_Y;
  G4double  transCFFrameHole4_X;
  G4double  transCFFrameHole4_Y;
  G4double  transCFFrameHole5_X;
  G4double  transCFFrameHole5_Y;
  G4double  transCFFrameHole6_X;
  G4double  transCFFrameHole6_Y;
  G4double  transCFFrameHole7_X;
  G4double  transCFFrameHole7_Y;

  G4double  transPoronFrameHole0_Y;
  G4double  transPoronFrameHole0_Z;
  G4double  transPoronFrameHole1_Y;
  G4double  transPoronFrameHole2_X;
  G4double  transPoronFrameHole2_Y;
  G4double  transPoronFrameHole3_X;
  G4double  transPoronFrameHole3_Y;
  G4double  transPoronFrameHole4_X;
  G4double  transPoronFrameHole4_Y;

  G4double  transRingKHole1_Y;
  G4double  transRingKHole1_Z;
  G4double  transRingKHole2_Y;
  G4double  transRingKHole3_Y;
  G4double  transRingSHole_Y;

  G4double transLadderBox_Y;
  G4double transLadderBox_Z;

  G4double transLadder_Y;
  G4double transLadder_Z;

  G4double transKaptonS_Y;
  G4double transKaptonS_Z;

  G4double  fPhysiTracker_X;
  G4double  fPhysiTracker_Y;
  G4double  fPhysiTracker_Z;
  
  G4double  fPhysiPoronFrame_X;
  G4double  fPhysiPoronFrame_Y;
  G4double  fPhysiPoronFrame_Z;
  
  G4double  fPhysiLadderBox_X;
  G4double  fPhysiLadderBox_Y;
  G4double  fPhysiLadderBox_Z;

  G4double  fPhysiRingK_X;
  G4double  fPhysiRingK_Y;
  G4double  fPhysiRingK_Z;

  G4double  fPhysiKaptonK_X;
  G4double  fPhysiKaptonK_Y;
  G4double  fPhysiKaptonK_Z;

  G4double  fPhysiHybridK_X;
  G4double  fPhysiHybridK_Y;
  G4double  fPhysiHybridK_Z;

  G4double  fPhysiSiliconPlateP_X;
  G4double  fPhysiSiliconPlateP_Y;
  G4double  fPhysiSiliconPlateP_Z;

  G4double  fPhysiSiliconSensorP_X;
  G4double  fPhysiSiliconSensorP_Y;
  G4double  fPhysiSiliconSensorP_Z;

  G4double  fPhysiSiliconPlateM_X;
  G4double  fPhysiSiliconPlateM_Y;
  G4double  fPhysiSiliconPlateM_Z;

  G4double  fPhysiSiliconSensorM_X;
  G4double  fPhysiSiliconSensorM_Y;
  G4double  fPhysiSiliconSensorM_Z;

  G4double  fPhysiKaptonS_X;
  G4double  fPhysiKaptonS_Y;
  G4double  fPhysiKaptonS_Z;

  G4double  fPhysiHybridS_X;
  G4double  fPhysiHybridS_Y;
  G4double  fPhysiHybridS_Z;

  G4double  fPhysiRingS_X;
  G4double  fPhysiRingS_Y;
  G4double  fPhysiRingS_Z;

  G4double  fPhysiHeatSink_X;
  G4double  fPhysiHeatSink_Y;
  G4double  fPhysiHeatSink_Z;

  HEPDSWMaterial* pMaterial;

  G4Box* fSolidTracker;
  G4Box* fSolidLayer;
  G4Box* fSolidCFFrameHole;
  G4SubtractionSolid* fSolidCFFrame;
  G4Box* fSolidPoronFrameHole;
  G4SubtractionSolid* fSolidPoronFrame;
  G4Box* fSolidSiliconPlate;
  G4Box* fSolidSiliconSensor;
  G4UnionSolid* fSolidLadderBox;
  G4UnionSolid* fSolidLadder;
  G4SubtractionSolid* fSolidKaptonS;
  G4Box* fSolidKaptonK;
  G4Box* fSolidHybrid;
  G4Box* fSolidHeatSink;
  G4SubtractionSolid* fSolidRingK;
  G4SubtractionSolid* fSolidRingS;
  
  G4LogicalVolume* fLogicTracker;
  G4LogicalVolume* fLogicLayer;
  G4LogicalVolume* fLogicCFFrame;
  G4LogicalVolume* fLogicPoronFrame;
  G4LogicalVolume* fLogicSiliconPlateP;
  G4LogicalVolume* fLogicSiliconPlateM;
  G4LogicalVolume* fLogicSiliconSensor;
  G4LogicalVolume* fLogicLadderBox;
  G4LogicalVolume* fLogicLadder;
  G4LogicalVolume* fLogicKaptonS;
  G4LogicalVolume* fLogicKaptonK;
  G4LogicalVolume* fLogicHybrid;
  G4LogicalVolume* fLogicHeatSink;
  G4LogicalVolume* fLogicRingK;
  G4LogicalVolume* fLogicRingS;

  G4VPhysicalVolume* fPhysiTracker;
  G4VPhysicalVolume* fPhysiLayer;
  G4VPhysicalVolume* fPhysiCFFrame;
  G4VPhysicalVolume* fPhysiPoronFrame;
  G4VPhysicalVolume* fPhysiLadderBox;
  G4VPhysicalVolume* fPhysiLadder;
  G4VPhysicalVolume* fPhysiRingK;
  G4VPhysicalVolume* fPhysiKaptonK;
  G4VPhysicalVolume* fPhysiHybridK;
  G4VPhysicalVolume* fPhysiSiliconPlateP;
  G4VPhysicalVolume* fPhysiSiliconSensorP;
  G4VPhysicalVolume* fPhysiSiliconPlateM;
  G4VPhysicalVolume* fPhysiSiliconSensorM;
  G4VPhysicalVolume* fPhysiKaptonS;
  G4VPhysicalVolume* fPhysiHybridS;
  G4VPhysicalVolume* fPhysiRingS;
  G4VPhysicalVolume* fPhysiHeatSink;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

