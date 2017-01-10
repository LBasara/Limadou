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

#include "CalorimeterConstruction.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterConstruction::CalorimeterConstruction()
{
  theConfig1 = new CalorimeterConstructionConfig1();
  theConfig2 = new CalorimeterConstructionConfig2();
  theConfig3 = new CalorimeterConstructionConfig3();
  theConfig4 = new CalorimeterConstructionConfig4();
  theConfig5 = new CalorimeterConstructionConfig5();
  theConfig6 = new CalorimeterConstructionConfig6();
  theDummy   = new CalorimeterConstructionDummy();
  theDummyOptical   = new CalorimeterConstructionDummyOptical();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterConstruction::~CalorimeterConstruction()
{
  delete theConfig1;
  delete theConfig2;
  delete theConfig3;
  delete theConfig4;
  delete theConfig5; 
  delete theConfig6;
  delete theDummy;
  delete theDummyOptical;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterConstruction::Builder(G4String config,G4VPhysicalVolume* motherVolume)
{
  if(config=="Config1")
    theConfig1->Builder(motherVolume);
  else if(config=="Config2")
    theConfig2->Builder(motherVolume);
  else if(config=="Config3")
    theConfig3->Builder(motherVolume);
  else if(config=="Config4")
    theConfig4->Builder(motherVolume);
  else if(config=="Config5")
    theConfig5->Builder(motherVolume); 
  else if(config=="Config6")
    theConfig6->Builder(motherVolume);
  else if(config=="Dummy")
    theDummy->Builder(motherVolume);
  else if(config=="DummyOptical")
    theDummyOptical->Builder(motherVolume);
  else
    G4cout<<"You have select an unavailable configuration"<<std::endl;
}

void CalorimeterConstruction::SetVetoMaterial(G4String aMat){
  theConfig1->SetVetoMaterial(aMat);
  theConfig2->SetVetoMaterial(aMat);
  theConfig3->SetVetoMaterial(aMat);
  theConfig4->SetVetoMaterial(aMat);
  theConfig5->SetVetoMaterial(aMat);
  theConfig6->SetVetoMaterial(aMat);

  theDummy->SetVetoMaterial(aMat);
}

// void CalorimeterConstruction::SetCaloMaterial(G4String aMat){
//   theConfig1->SetCaloMaterial(aMat);
//   theConfig2->SetCaloMaterial(aMat);
//   theConfig3->SetCaloMaterial(aMat);
//   theConfig4->SetCaloMaterial(aMat);
//   theConfig5->SetCaloMaterial(aMat); 
//   theConfig6->SetCaloMaterial(aMat);
//   theDummy->SetCaloMaterial(aMat);
//   theDummyOptical->SetCaloMaterial(aMat);
// }

void CalorimeterConstruction::SetCaloMaterial(G4String aMat1,G4String aMat2){
  theConfig3->SetCaloMaterial(aMat1,aMat2);
  theConfig4->SetCaloMaterial(aMat1,aMat2);
  theConfig5->SetCaloMaterial(aMat1,aMat2);
  theConfig6->SetCaloMaterial(aMat1,aMat2);
  theDummy->SetCaloMaterial(aMat1,aMat2);
}

void CalorimeterConstruction::SetNumberOfCrystalLayer(G4int aVal){
  theConfig3->SetNumberOfCrystalLayer(aVal);
  theDummy->SetNumberOfCrystalLayer(aVal);
}

void CalorimeterConstruction::SetPoronMaterial(G4String aMat){
  theConfig1->SetPoronMaterial(aMat);
  theConfig2->SetPoronMaterial(aMat);
  theConfig3->SetPoronMaterial(aMat);
  theConfig4->SetPoronMaterial(aMat);
  theConfig5->SetPoronMaterial(aMat);
  theConfig6->SetPoronMaterial(aMat);

}
void CalorimeterConstruction::SetCarbonFiberMaterial(G4String aMat){
  theConfig1->SetCarbonFiberMaterial(aMat);
  theConfig2->SetCarbonFiberMaterial(aMat);
  theConfig3->SetCarbonFiberMaterial(aMat);
  theConfig4->SetCarbonFiberMaterial(aMat);
  theConfig5->SetCarbonFiberMaterial(aMat);
  theConfig6->SetCarbonFiberMaterial(aMat);

}

void CalorimeterConstruction::SetHoneyCombMaterial(G4String aMat){
  theConfig1->SetHoneyCombMaterial(aMat);
  theConfig2->SetHoneyCombMaterial(aMat);
  theConfig3->SetHoneyCombMaterial(aMat);
  theConfig4->SetHoneyCombMaterial(aMat);
  theConfig5->SetHoneyCombMaterial(aMat);
  theConfig6->SetHoneyCombMaterial(aMat);

}


