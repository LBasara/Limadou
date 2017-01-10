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

#include "TrackerConstruction.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerConstruction::TrackerConstruction()
{
  theConfig1 = new TrackerConstructionConfig1();
  theConfig2 = new TrackerConstructionConfig2();
  theConfig3 = new TrackerConstructionConfig3();
  theConfig4 = new TrackerConstructionConfig4();
  theDummy   = new TrackerConstructionDummy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerConstruction::~TrackerConstruction()
{
  delete theConfig1;
  delete theConfig2;
  delete theConfig3;
  delete theConfig4;
  delete theDummy;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerConstruction::Builder(G4String config,G4VPhysicalVolume* motherVolume)
{
  if(config=="Config1")
    theConfig1->Builder(motherVolume);
  else if(config=="Config2")
    theConfig2->Builder(motherVolume);
  else if(config=="Config3")
    theConfig3->Builder(motherVolume);
  else if(config=="Config4")
    theConfig4->Builder(motherVolume);
  else if(config=="Dummy")
    theDummy->Builder(motherVolume);
  else
    G4cout<<"You have select an unavailable configuration"<<std::endl;
}


void TrackerConstruction::SetSiliconMaterial(G4String aMat){
  theConfig1->SetSiliconMaterial(aMat);
  theConfig2->SetSiliconMaterial(aMat);
  theConfig3->SetSiliconMaterial(aMat);
  theConfig4->SetSiliconMaterial(aMat);
  theDummy->SetSiliconMaterial(aMat);
}

void TrackerConstruction::SetKaptonMaterial(G4String aMat){
  theConfig1->SetKaptonMaterial(aMat);
  theConfig2->SetKaptonMaterial(aMat);
  theConfig3->SetKaptonMaterial(aMat);
  theConfig4->SetKaptonMaterial(aMat);

}
void TrackerConstruction::SetCarbonFiberMaterial(G4String aMat){
  theConfig1->SetCarbonFiberMaterial(aMat);
  theConfig2->SetCarbonFiberMaterial(aMat);
  theConfig3->SetCarbonFiberMaterial(aMat);
  theConfig4->SetCarbonFiberMaterial(aMat);

}
void TrackerConstruction::SetPoronMaterial(G4String aMat){
  theConfig4->SetPoronMaterial(aMat);
}
