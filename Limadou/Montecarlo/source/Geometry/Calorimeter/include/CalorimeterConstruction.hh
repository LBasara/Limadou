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
/// \file electromagnetic/TestEm3/include/CalorimeterConstruction.hh
/// \brief Definition of the CalorimeterConstruction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CalorimeterConstruction_h
#define CalorimeterConstruction_h 1

#include "globals.hh"
#include "CalorimeterConstructionConfig1.hh"
#include "CalorimeterConstructionConfig2.hh"
#include "CalorimeterConstructionConfig3.hh"
#include "CalorimeterConstructionConfig4.hh"
#include "CalorimeterConstructionConfig5.hh"
#include "CalorimeterConstructionConfig6.hh"
#include "CalorimeterConstructionDummy.hh"
#include "CalorimeterConstructionDummyOptical.hh"

class G4IntersectionSolid;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CalorimeterConstruction 
{
public:
  
  CalorimeterConstruction();
  ~CalorimeterConstruction();
  
  void SetVetoMaterial(G4String aMat);
  //  void SetCaloMaterial(G4String aMat);
  void SetCaloMaterial(G4String aMat1,G4String aMat2);
  void SetPoronMaterial(G4String aMat);
  void SetCarbonFiberMaterial(G4String aMat);
  void SetHoneyCombMaterial(G4String aMat);
  void SetNumberOfCrystalLayer(G4int aVal);
  void Builder(G4String config,G4VPhysicalVolume* motherVolume);
  
private:
  CalorimeterConstructionConfig1* theConfig1;
  CalorimeterConstructionConfig2* theConfig2;
  CalorimeterConstructionConfig3* theConfig3;
  CalorimeterConstructionConfig4* theConfig4;
  CalorimeterConstructionConfig5* theConfig5;
  CalorimeterConstructionConfig6* theConfig6;
  CalorimeterConstructionDummy* theDummy;
  CalorimeterConstructionDummyOptical* theDummyOptical;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif














				
				
				
