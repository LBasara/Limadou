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
/// \file electromagnetic/TestEm3/include/TrackerConstruction.hh
/// \brief Definition of the TrackerConstruction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef TrackerConstruction_h
#define TrackerConstruction_h 1

#include "globals.hh"
#include "TrackerConstructionConfig1.hh"
#include "TrackerConstructionConfig2.hh"
#include "TrackerConstructionConfig3.hh"
#include "TrackerConstructionConfig4.hh"
#include "TrackerConstructionDummy.hh"

class G4IntersectionSolid;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackerConstruction 
{
public:
  
  TrackerConstruction();
  ~TrackerConstruction();

  void SetSiliconMaterial(G4String aMat);
  void SetKaptonMaterial(G4String aMat);
  void SetCarbonFiberMaterial(G4String aMat);
  void SetPoronMaterial(G4String aMat);

 
  void Builder(G4String config,G4VPhysicalVolume* motherVolume);
  
private:
  TrackerConstructionConfig1* theConfig1;
  TrackerConstructionConfig2* theConfig2;
  TrackerConstructionConfig3* theConfig3;
  TrackerConstructionConfig4* theConfig4;
  TrackerConstructionDummy* theDummy;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif














				
				
				
