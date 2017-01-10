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

#include "HEPDBoxConstruction.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HEPDBoxConstruction::HEPDBoxConstruction()
{
  theConfig1 = new HEPDBoxConstructionConfig1();
  theConfig2 = new HEPDBoxConstructionConfig2();
  theConfig3 = new HEPDBoxConstructionConfig3();
  theConfig4 = new HEPDBoxConstructionConfig4();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HEPDBoxConstruction::~HEPDBoxConstruction()
{
  delete theConfig1;
  delete theConfig2;
  delete theConfig3;
  delete theConfig4;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HEPDBoxConstruction::Builder(G4String config,G4VPhysicalVolume* motherVolume)
{
  if(config=="Config1")
    theConfig1->Builder(motherVolume);
  else if(config=="Config2")
    theConfig2->Builder(motherVolume);
  else if(config=="Config3")
    theConfig3->Builder(motherVolume);
  else if(config=="Config4")
    theConfig4->Builder(motherVolume);
  else
    G4cout<<"You have select an unavailable configuration"<<std::endl;
}


void HEPDBoxConstruction::SetBlanketMaterial(G4String aMat){
  theConfig1->SetBlanketMaterial(aMat);
  theConfig2->SetBlanketMaterial(aMat);
  theConfig4->SetBlanketMaterial(aMat);
}

void HEPDBoxConstruction::SetBlanketMaterial(G4String aMat1, G4String aMat2, G4String aMat3,G4String aMat4){
  theConfig3->SetBlanketMaterial(aMat1,aMat2,aMat3,aMat4);
}

void HEPDBoxConstruction::SetWallMaterial(G4String aMat){
  theConfig1->SetWallMaterial(aMat);
}
void HEPDBoxConstruction::SetWallMaterial(G4String aMat1,G4String aMat2){
  theConfig4->SetWallMaterial(aMat1,aMat2);
}
