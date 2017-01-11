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
/// \file electromagnetic/TestEm3/include/CalorimeterConstructionDummyOptical.hh
/// \brief Definition of the CalorimeterConstructionDummyOptical class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CalorimeterConstructionDummyOptical_h
#define CalorimeterConstructionDummyOptical_h 1

#include "globals.hh"
#include "G4AssemblyVolume.hh"

class G4IntersectionSolid;
class G4Box;
class G4UnionSolid;
class G4SubtractionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class HEPDSWMaterial;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CalorimeterConstructionDummyOptical 
{
public:
  
   CalorimeterConstructionDummyOptical();
  ~CalorimeterConstructionDummyOptical();

  inline void SetVetoMaterial(G4String aMat){vetoMaterial=aMat;}
  inline void SetCaloMaterial(G4String aMat){scintMaterial=aMat;}
  inline void SetCaloMaterial(G4String aMat1,G4String aMat2){scintMaterial=aMat1;crystalMaterial=aMat2;}

  void SetNumberOfCrystalLayer(G4int aVal);

  void Builder(G4VPhysicalVolume* motherVolume);
  
private:
  void ComputeObjectsPositioning();

  G4int              fNbOfCaloLayers;
  G4int              fNbOfScintLayers;
  G4int              fNbOfCrystalLayers;

  G4String           vetoMaterial;
  G4String           scintMaterial;
  G4String           crystalMaterial;

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

  G4double fVetoLatX_X;
  G4double fVetoLatX_Y;
  G4double fVetoLatX_Z;

  G4double fVetoLatY_X;
  G4double fVetoLatY_Y;
  G4double fVetoLatY_Z;

  G4double fVetoBottom_X;
  G4double fVetoBottom_Y;
  G4double fVetoBottom_Z;

  G4double fVetoUp_X;
  G4double fVetoUp_Y;
  G4double fVetoUp_Z;


  G4double fPhysiVetoLatX_X;
  G4double fPhysiVetoLatX_Y;
  G4double fPhysiVetoLatX_Z;

  G4double fPhysiVetoLatY_X;
  G4double fPhysiVetoLatY_Y;
  G4double fPhysiVetoLatY_Z;

  G4double fPhysiVetoBottom_X;
  G4double fPhysiVetoBottom_Y;
  G4double fPhysiVetoBottom_Z;

  G4double fPhysiVetoUp_X;
  G4double fPhysiVetoUp_Y;
  G4double fPhysiVetoUp_Z;

  G4double  fPhysiScintBox_X;
  G4double  fPhysiScintBox_Y;
  G4double  fPhysiScintBox_Z;

  G4double  fPhysiCrystalBox_X;
  G4double  fPhysiCrystalBox_Y;
  G4double  fPhysiCrystalBox_Z;

  HEPDSWMaterial*      pMaterial;

  G4Box* fSolidCaloBox;
  G4Box* fSolidScintBox;
  G4Box* fSolidCrystalBox;
  G4Box* fSolidLayer;
  G4Box* fSolidVetoLatX;	   
  G4Box* fSolidVetoLatY;	   
  G4Box* fSolidVetoBottom;	   
  G4SubtractionSolid* fSolidVetoUp;	   

  G4LogicalVolume* fLogicCaloBox;
  G4LogicalVolume* fLogicScintBox;
  G4LogicalVolume* fLogicCrystalBox;
  G4LogicalVolume* fLogicLayer;
  G4LogicalVolume* fLogicScintLayer;
  G4LogicalVolume* fLogicCrystalLayer;
  G4LogicalVolume* fLogicVetoLatX;	   
  G4LogicalVolume* fLogicVetoLatY;	   
  G4LogicalVolume* fLogicVetoBottom;	   
  G4LogicalVolume* fLogicVetoUp;	   

  G4VPhysicalVolume* fPhysiCaloBox;
  G4VPhysicalVolume* fPhysiScintBox;
  G4VPhysicalVolume* fPhysiCrystalBox;
  G4VPhysicalVolume* fPhysiLayer;
//  G4VPhysicalVolume* fPhysiScintLayer;
  G4AssemblyVolume* fPhysiScintLayer;
  G4VPhysicalVolume* fPhysiCrystalLayer;
  G4VPhysicalVolume* fPhysiVetoLatXP;
  G4VPhysicalVolume* fPhysiVetoLatXM;
  G4VPhysicalVolume* fPhysiVetoLatYP;
  G4VPhysicalVolume* fPhysiVetoLatYM;
  G4VPhysicalVolume* fPhysiVetoBottom;
  G4VPhysicalVolume* fPhysiVetoUp;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif














				
				
				
