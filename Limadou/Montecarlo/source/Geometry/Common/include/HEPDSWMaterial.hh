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
//
// $Id$
//
//    **********************************
//    *                                *
//    *      BrachyMaterial.hh          *
//    *                                *
//    **********************************
//
//Code developed by: Susanna Guatelli, guatelli@ge.infn.it
//
//This class manages the elements and materials needed by the simulation
// set-up ...
//
#ifndef HEPDSWMaterial_H
#define HEPDSWMaterial_H 1
#include "globals.hh"
class G4Material;
class G4MaterialPropertiesTable;
class HEPDSWMaterial
{ 
public:
  HEPDSWMaterial();
  ~ HEPDSWMaterial();

public:
  void  DefineMaterials();
  G4Material* GetMaterial(G4String); //returns the material

private:
  G4Material* matAir;
  G4Material* matAl; 
  G4Material* matTi; 
  G4Material* matGraphite; 
  G4Material* matCarbonFiber; 
  G4Material* matH2O;
  G4Material* matMgB2;
  G4Material* matGlass;
  G4Material* matCu;
  G4Material* nylon;  
  G4Material* teflon;
  G4Material* mylar;
  G4Material* beta;
  G4Material* nextel;
  G4Material* kevlar;
  G4Material* vacuum;
  G4Material* betaCloth; 
  G4Material* eterogeneousNextel;
  G4Material* kevlarVacuum;
  G4Material* polyethylene;
  G4Material* polyacrylate;
  G4Material* evoh;
  G4Material* nomex; 
  G4Material* matPoly;
  G4Material* matSci;
  G4Material* matSciOpt;
  G4Material* matLct;
  G4Material* matSili;
  G4Material* matSiO2;
  G4Material* matG10;
  G4Material* matCsI;
  G4Material* matBGO;
  G4Material* matSiNx;
  G4Material* matlAr;
  G4Material* matlArEm3;
  G4Material* matLead;
  G4Material* matLeadSb;
  G4Material* matLeadGlass;
  G4Material* matLYSO;
  G4Material* matLYSOOpt;
  G4Material* matAerog;
  G4Material* matCO2;
  G4Material* matPoron;
  G4Material* matKapton;
  G4Material* matDacron;
  G4Material* matEpoxy;
  G4Material* matFR4;
  G4Material* matHCAlu;
  G4Material* matEJ200;      
  //  G4Material* matTeflon;

  //  G4MaterialPropertiesTable* MPTMylar;  
  //  G4MaterialPropertiesTable* MPTTeflon;
  G4MaterialPropertiesTable* MPTEJ200;
  G4MaterialPropertiesTable* MPTScintillator;
  G4MaterialPropertiesTable* MPTLyso;

};
#endif
