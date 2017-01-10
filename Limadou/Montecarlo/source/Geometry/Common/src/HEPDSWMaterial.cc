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
// Code developed by:
//  F.Ambroglini
//
//    *******************************
//    *                             *
//    *    HEPDSWMaterial.cc        *
//    *                             *
//    *******************************
//
// $Id$
//

#include "HEPDSWMaterial.hh"

#include "globals.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4MaterialTable.hh"
#include "G4RunManager.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"

HEPDSWMaterial::HEPDSWMaterial():
  matAir(0), matAl(0), matTi(0), matGraphite(0),matCarbonFiber(0),matH2O(0),matMgB2(0),matGlass(0),matCu(0),
  nylon(0), teflon(0), mylar(0), 
  beta(0), nextel(0), kevlar(0),
  vacuum(0), betaCloth(0), eterogeneousNextel(0), kevlarVacuum(0),
  polyethylene(0), polyacrylate(0), evoh(0), nomex(0), 
  matPoly(0),matSci(0),matSciOpt(0),matLct(0),matSili(0),matSiO2(0),matG10(0),
  matCsI(0),matBGO(0),matSiNx(0),matlAr(0),matlArEm3(0),matLead(0),matLeadSb(0),matLeadGlass(0),matLYSO(0),
  matLYSOOpt(0),matAerog(0),matCO2(0),matPoron(0),matKapton(0),matDacron(0),matEpoxy(0),matFR4(0),matHCAlu(0),MPTScintillator(0),MPTLyso(0)
{;}

HEPDSWMaterial::~HEPDSWMaterial()
{
  delete nomex;
  delete evoh;
  delete polyacrylate;
  delete polyethylene;
  delete kevlarVacuum;
  delete eterogeneousNextel;
  delete betaCloth;
  delete vacuum;
  delete kevlar;
  delete nextel;
  delete beta;
  delete mylar;
  delete teflon;
  delete nylon;
  delete matH2O;
  delete matAir;
  delete matAl;
  delete matTi;
  delete matGraphite;
  delete matGlass;
  delete matCu;
  delete matPoly;
  delete matSci;
  delete matSciOpt;
  delete matLct;
  delete matSili;
  delete matSiO2;
  delete matG10;
  delete matCsI;
  delete matBGO;
  delete matSiNx;
  delete matlAr;
  delete matlArEm3;
  delete matLead;
  delete matLeadSb;
  delete matLeadGlass;
  delete matLYSO;
  delete matLYSOOpt;
  delete matAerog;
  delete matCO2;
  delete matPoron;
  delete matKapton;
  delete matDacron;
  delete matEpoxy;
  delete matFR4;
  delete matHCAlu;

}

void HEPDSWMaterial::DefineMaterials()
{

  // This function illustrates the possible ways to define materials using 
  // G4 database on G4Elements
  G4NistManager* manager = G4NistManager::Instance();
  manager->SetVerbose(0);
  //
  // define Elements
  //
  G4Element* elH  = manager->FindOrBuildElement(1);
  G4Element* elB  = manager->FindOrBuildElement(5);
  G4Element* elC  = manager->FindOrBuildElement(6);
  G4Element* elN  = manager->FindOrBuildElement(7);
  G4Element* elO  = manager->FindOrBuildElement(8);
  G4Element* elF  = manager->FindOrBuildElement(9);
  //G4Element* elNa = manager->FindOrBuildElement(11);
  G4Element* elMg = manager->FindOrBuildElement(12);
  G4Element* elAl = manager->FindOrBuildElement(13);
  G4Element* elSi = manager->FindOrBuildElement(14);
  //G4Element* elP  = manager->FindOrBuildElement(15);  
  //G4Element* elS  = manager->FindOrBuildElement(16);  
  G4Element* elCl = manager->FindOrBuildElement(17);  
  //G4Element* elK  = manager->FindOrBuildElement(19);  
  //G4Element* elCa = manager->FindOrBuildElement(20);  
  G4Element* elTi = manager->FindOrBuildElement(22);  
  //G4Element* elFe = manager->FindOrBuildElement(26);  
  //G4Element* elCu = manager->FindOrBuildElement(29);  
  //G4Element* elZn = manager->FindOrBuildElement(30);  
  G4Element* elGe = manager->FindOrBuildElement(32);
  G4Element* elAs = manager->FindOrBuildElement(33);
  //G4Element* elRb = manager->FindOrBuildElement(37);
  //G4Element* elSr = manager->FindOrBuildElement(38);
  G4Element* elY = manager->FindOrBuildElement(39);
  //G4Element* elZr = manager->FindOrBuildElement(40);
  G4Element* elSb = manager->FindOrBuildElement(51);
  G4Element* elI  = manager->FindOrBuildElement(53);
  G4Element* elCs = manager->FindOrBuildElement(55);
  G4Element* elCe = manager->FindOrBuildElement(58);
  G4Element* elLu = manager->FindOrBuildElement(71);
  G4Element* elPb = manager->FindOrBuildElement(82);
  G4Element* elBi = manager->FindOrBuildElement(83);

  //
  // define an Element from isotopes, by relative abundance
  //
  G4int iz, n;                       //iz=number of protons  in an isotope;
                                     // n=number of nucleons in an isotope;
  G4int   ncomponents;                                     
  G4double abundance;                                     

  G4Isotope* U5 = new G4Isotope("U235", iz=92, n=235,235.01*g/mole);
  G4Isotope* U8 = new G4Isotope("U238", iz=92, n=238,238.03*g/mole);

  G4Element* U  = new G4Element("enriched Uranium", "U", ncomponents=2);
  U->AddIsotope(U5, abundance= 90.*perCent);
  U->AddIsotope(U8, abundance= 10.*perCent);

  // Define required materials

  //  G4double A;  // atomic mass
  //G4double Z;  // atomic number
  G4double d;  // density
  
  // Air material
  matAir = manager->FindOrBuildMaterial("G4_AIR");
  // Aluminium
  matAl =  new G4Material("Aluminium",13.,26.98*g/mole,2.700*g/cm3);  
  //Titanium
  matTi = new G4Material("Titanium",22.,47.867*g/mole,4.54*g/cm3);
  //MgB2
  d = 2.500*g/cm3;
  matMgB2 = new G4Material("MgB2",d,2);
  matMgB2->AddElement(elMg,1);
  matMgB2->AddElement(elB,2);
  //Graphite
  d = 2.210*g/cm3;
  matGraphite = new G4Material("Graphite",d,1);
  matGraphite->AddElement(elC,1);
  //CarbonFiber
  d = 1.55*g/cm3;
  matCarbonFiber = new G4Material("CarbonFiber",d,1);
  matCarbonFiber->AddElement(elC,1);
  //Glass
  d = 1.220*g/cm3;
  matGlass = new G4Material("Glass",d,2);
  matGlass->AddElement(elSi,1);
  matGlass->AddElement(elO,2);
  //Copper
  matCu = new G4Material("Copper",29.,63.55*g/mole,8.960*g/cm3);
  // Water
  d = 1.000*g/cm3;
  matH2O = new G4Material("Water",d,2);
  matH2O->AddElement(elH,2);
  matH2O->AddElement(elO,1);
  matH2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  matH2O->SetChemicalFormula("H_2O");
  
   //Teflon
  d = 2.2 *g/cm3;
  teflon = new G4Material("Teflon",d,2);
  teflon -> AddElement(elF,0.759817);
  teflon -> AddElement(elC,0.240183);

  //nylon (alenia spazio)
  d = 1.14 *g/cm3;
  nylon = new G4Material("nylon",d,4);
  nylon -> AddElement(elH,0.108);
  nylon -> AddElement(elC,0.68);
  nylon -> AddElement(elN,0.099);
  nylon -> AddElement(elO,0.113);

  //Mylar (alenia spazio)
  d= 1.4 *g/cm3;
  mylar = new G4Material("mylar",d,3);
  mylar -> AddElement(elH,0.042);
  mylar -> AddElement(elC,0.625);
  mylar -> AddElement(elO,0.333);

  //beta cloth
  G4double betaDensity = 2.3 *g/cm3;
  beta = new G4Material("beta",betaDensity,2);
  beta -> AddElement(elO,0.53);
  beta -> AddElement(elSi,0.47);
 
  G4double nextelDensity = 2.7 * g/cm3;
  nextel = new G4Material("nextel",nextelDensity,4);
  nextel -> AddElement(elB,0.04);
  nextel -> AddElement(elO,0.52);
  nextel -> AddElement(elAl,0.33);
  nextel -> AddElement(elSi,0.11);

  //kevlar
  G4double kevlarDensity = 1.44 *g/cm3;
  kevlar = new G4Material("kevlar",d,4);
  kevlar -> AddElement(elH,0.04);
  kevlar -> AddElement(elC,0.71);
  kevlar -> AddElement(elO,0.12);
  kevlar -> AddElement(elN,0.13);
  
  //vacuum
  G4double vacuumDensity = 1.e-25 *g/cm3;
  G4double pressure = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  vacuum = new G4Material("Galactic",1.,1.008*g/mole,vacuumDensity,kStateGas,temperature,pressure);
  d = (vacuumDensity*0.44)+ (betaDensity*0.56);
  betaCloth = new G4Material("betacloth",d,2);
  betaCloth -> AddMaterial (beta, 0.56);
  betaCloth -> AddMaterial (vacuum,0.44);

  d = (vacuumDensity*0.73)+ (nextelDensity*0.27); 
  eterogeneousNextel = new G4Material("Nextel312AF62",d,2);
  eterogeneousNextel -> AddMaterial (vacuum, 0.73);
  eterogeneousNextel -> AddMaterial (nextel, 0.27);

  d = (vacuumDensity*0.44)+ (kevlarDensity*0.56);
  kevlarVacuum = new G4Material("kevlarVacuum",d,2);
  kevlarVacuum -> AddMaterial (vacuum, 0.44);
  kevlarVacuum -> AddMaterial (kevlar, 0.56);

  d = 0.94 * g/cm3;
  polyethylene = new G4Material("polyethylene",d,2);
  polyethylene -> AddElement(elH,0.14);
  polyethylene -> AddElement(elC,0.86);
  
  d = 1.19 * g/cm3;
  polyacrylate = new G4Material("polyacrylate",d,3);
  polyacrylate -> AddElement(elH,0.08);
  polyacrylate -> AddElement(elC,0.60);
  polyacrylate -> AddElement(elO,0.32);
 
  d = 1.17 * g/cm3;
  evoh = new G4Material("evoh",d,3);
  evoh -> AddElement(elH,0.11);
  evoh -> AddElement(elC,0.67);
  evoh -> AddElement(elO,0.22);
  
  G4double nomexDensity = 0.98 * g/cm3;
  nomex = new G4Material("nomex",nomexDensity,5);
  nomex -> AddElement(elH,0.04);
  nomex -> AddElement(elC,0.54);
  nomex -> AddElement(elN,0.09);
  nomex -> AddElement(elO,0.10);
  nomex -> AddElement(elCl,0.23);

  matPoly = new G4Material("Polystyrene",1.032*g/cm3,2);
  matPoly->AddElement(elC,1);
  matPoly->AddElement(elH,1);

  matSci = new G4Material("Scintillator",1.032*g/cm3,2);
  matSci->AddElement(elC,19);
  matSci->AddElement(elH,21);

  matSciOpt = new G4Material("ScintillatorOpt",1.032*g/cm3,2);
  matSciOpt->AddElement(elC,19);
  matSciOpt->AddElement(elH,21);
  const G4int scinum = 14;
  G4double sci_Energy[scinum] ={2.25*eV,2.38*eV,2.48*eV,2.58*eV,2.70*eV,2.79*eV,2.82*eV,2.92*eV,2.95*eV,3.10*eV,3.18*eV,3.26*eV,3.44*eV,3.54*eV};
  G4double sci_RIND[scinum]={1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58};
  G4double sci_ABSL[scinum]={210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm,210*cm};
  G4double sci_SCINT[scinum]={0,0.02,0.08,0.2,0.4,0.6,0.8,1,0.8,0.2,0.1,0.07,0.02,0};
  MPTScintillator = new G4MaterialPropertiesTable();
  MPTScintillator->AddProperty("RINDEX",sci_Energy,sci_RIND,scinum);
  MPTScintillator->AddProperty("ABSLENGTH",sci_Energy,sci_ABSL,scinum);
  MPTScintillator->AddProperty("FASTCOMPONENT",sci_Energy, sci_SCINT,scinum);
  MPTScintillator->AddConstProperty("SCINTILLATIONYIELD",13./keV);
  MPTScintillator->AddConstProperty("RESOLUTIONSCALE",1.0);
  MPTScintillator->AddConstProperty("FASTTIMECONSTANT", 2.1*ns);
  matSciOpt->SetMaterialPropertiesTable(MPTScintillator);



// EJ-200   calorimeter scintillation material  added 01/12/14
  
  //matEJ200 = manager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  matEJ200 = new G4Material("ScintillatorEJ200",1.032*g/cm3,2);
  matEJ200->AddElement(elC,19);
  matEJ200->AddElement(elH,21);

  matEJ200->GetIonisation()->SetBirksConstant(0.16*mm/MeV);
  
  const G4int ej200num = 55;
  G4double ej200_Energy[ej200num]={2.481*eV, 2.487*eV, 2.498*eV, 2.515*eV, 2.527*eV, 2.544*eV, 2.558*eV, 2.572*eV, 2.585*eV, 2.594*eV, 2.608*eV, 2.622*eV, 2.638*eV, 2.655*eV, 2.674*eV, 2.687*eV, 2.708*eV, 2.743*eV, 2.762*eV, 2.77*eV, 2.787*eV, 2.798*eV, 2.807*eV, 2.822*eV, 2.831*eV, 2.843*eV, 2.859*eV, 2.88*eV, 2.894*eV, 2.902*eV, 2.916*eV, 2.926*eV, 2.935*eV, 2.94*eV, 2.947*eV, 2.952*eV, 2.955*eV, 2.961*eV, 2.966*eV, 2.971*eV, 2.977*eV, 2.982*eV, 2.989*eV, 2.996*eV, 3.002*eV, 3.007*eV, 3.018*eV, 3.023*eV, 3.034*eV, 3.044*eV, 3.055*eV, 3.07*eV, 3.083*eV, 3.102*eV, 3.124*eV};
  G4double ej200_SCINT[ej200num] = {0.05845, 0.06153, 0.07076, 0.08307, 0.09845, 0.11691, 0.1323, 0.15383, 0.17537, 0.19691, 0.22767, 0.26152, 0.30767, 0.3569, 0.39997, 0.42151, 0.45227, 0.52303, 0.56918, 0.5938, 0.64918, 0.68302, 0.71686, 0.76917, 0.80609, 0.85224, 0.90147, 0.95069, 0.98147, 0.99685, 0.99992, 0.99377, 0.98147, 0.94762, 0.90762, 0.8707, 0.84301, 0.80301, 0.73841, 0.68918, 0.63687, 0.58765, 0.53226, 0.46458, 0.41843, 0.36305, 0.30459, 0.26152, 0.20306, 0.16614, 0.11076, 0.04615, 0.01846, 0.00923, 0.00308};
  G4double ej200_RIND[ej200num]  = {1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,
  1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58};
  
  // Attenuation length found between 300 and 400 cm; for the moment we take 380*cm as from NIM A 552 (2005) 449
  
  G4double ej200_ABSL[ej200num]  = { 380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm,380*cm};

  MPTEJ200 = new G4MaterialPropertiesTable();
  MPTEJ200->AddProperty("FASTCOMPONENT", ej200_Energy, ej200_SCINT, ej200num)->SetSpline(true);
  MPTEJ200->AddProperty("RINDEX",        ej200_Energy, ej200_RIND,  ej200num)->SetSpline(true);
  MPTEJ200->AddProperty("ABSLENGTH",     ej200_Energy, ej200_ABSL,  ej200num)->SetSpline(true);
  MPTEJ200->AddConstProperty("SCINTILLATIONYIELD",10574/MeV); // CHECK THE CORRECT NUMBER
  MPTEJ200->AddConstProperty("RESOLUTIONSCALE",1.0);
  MPTEJ200->AddConstProperty("FASTTIMECONSTANT",2.1*ns);
  MPTEJ200->AddConstProperty("SLOWTIMECONSTANT",2.1*ns);
  MPTEJ200->AddConstProperty("SLOWSCINTILLATIONRISETIME",0.9*ns);
  MPTEJ200->AddConstProperty("FASTSCINTILLATIONRISETIME",0.9*ns);
  MPTEJ200->AddConstProperty("YIELDRATIO",1.0);
  matEJ200->SetMaterialPropertiesTable(MPTEJ200);

///






  matLct = new G4Material("Lucite",1.185*g/cm3,3);
  matLct->AddElement(elC, 59.97*perCent);
  matLct->AddElement(elH, 8.07*perCent);
  matLct->AddElement(elO, 31.96*perCent);

  matSili = new G4Material("Silicon",2.330*g/cm3,1);
  matSili->AddElement(elSi,1);

  matSiO2 = new G4Material("quartz",2.200*g/cm3,2);
  matSiO2->AddElement(elSi,1);
  matSiO2->AddElement(elO,2);

  matG10 = new G4Material("NemaG10",1.700*g/cm3,4);
  matG10->AddElement(elSi, 1);
  matG10->AddElement(elO , 2);
  matG10->AddElement(elC , 3);
  matG10->AddElement(elH , 3);

  matCsI = new G4Material("CsI",4.534*g/cm3,2);
  matCsI->AddElement(elCs, 1);
  matCsI->AddElement(elI , 1);
  matCsI->GetIonisation()->SetMeanExcitationEnergy(553.1*eV);

  matBGO = new G4Material("BGO",7.10*g/cm3, 3);
  matBGO->AddElement(elO , 12);
  matBGO->AddElement(elGe,  3);
  matBGO->AddElement(elBi,  4);

  matSiNx= new G4Material("SiNx", 3.1 *g/cm3,3);
  matSiNx-> AddElement(elSi, 300);
  matSiNx-> AddElement(elN, 310);
  matSiNx-> AddElement(elH, 6);

  matlAr = manager->FindOrBuildMaterial("G4_lAr");

  matlArEm3 = new G4Material("liquidArgon", 1.390*g/cm3,1);
  matlArEm3->AddMaterial(matlAr,1.0);
  
  matLead = new G4Material("Lead",11.35*g/cm3,1);
  matLead->AddElement(elPb,1.0);

  matLeadSb = new G4Material("LeadSb",11.35*g/cm3,2);
  matLeadSb->AddElement(elSb,4.*perCent);
  matLeadSb->AddElement(elPb,96.*perCent);

  matLYSO = new G4Material("LYSO", 7.1*g/cm3, 5, kStateSolid);
  matLYSO->AddElement(elLu, 71.43*perCent);
  matLYSO->AddElement(elY, 4.03*perCent);
  matLYSO->AddElement(elSi, 6.37*perCent);
  matLYSO->AddElement(elO, 18.14*perCent);
  matLYSO->AddElement(elCe, 0.02*perCent);

  matLYSOOpt = new G4Material("LYSOOpt", 7.1*g/cm3, 5, kStateSolid);
  matLYSOOpt->AddElement(elLu, 71.43*perCent);
  matLYSOOpt->AddElement(elY, 4.03*perCent);
  matLYSOOpt->AddElement(elSi, 6.37*perCent);
  matLYSOOpt->AddElement(elO, 18.14*perCent);
  matLYSOOpt->AddElement(elCe, 0.02*perCent);


  const G4int lysonum = 19;
  G4double lyso_Energy[lysonum]={1.91*eV,2.07*eV,2.16*eV,2.25*eV,2.37*eV,2.48*eV,2.68*eV,2.76*eV,2.86*eV,2.95*eV,2.97*eV,3.02*eV,3.05*eV,3.08*eV,3.16*eV,3.21*eV,3.28*eV,3.35*eV,3.54};
  G4double lyso_SCINT[lysonum] = {0,0.04,0.08,0.16,0.32,0.48,0.84,0.94,0.96,1,0.92,0.84,0.72,0.6,0.36,0.2,0.08,0.04,0.00};
  G4double lyso_RIND[lysonum]  = {1.81,1.81,1.81,1.81,1.81,1.81,1.81,1.81,1.81,1.81,1.81,1.81,1.81,1.81,1.81,1.81,1.81,1.81,1.81};
  //Attenuation length come from this paper "Optimization of the effective light attenuation length of YAP:Ce and LYSO:Ce crystals for a novel geometrical PET concept"
  G4double lyso_ABSL[lysonum]  = { 20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm,20*cm};

  MPTLyso = new G4MaterialPropertiesTable();
  MPTLyso->AddProperty("FASTCOMPONENT", lyso_Energy, lyso_SCINT, lysonum);
  MPTLyso->AddProperty("RINDEX",        lyso_Energy, lyso_RIND,  lysonum);
  MPTLyso->AddProperty("ABSLENGTH",     lyso_Energy, lyso_ABSL,  lysonum);
  MPTLyso->AddConstProperty("SCINTILLATIONYIELD",32/keV);
  MPTLyso->AddConstProperty("RESOLUTIONSCALE",1.0);
  MPTLyso->AddConstProperty("FASTTIMECONSTANT",42.*ns);
  matLYSOOpt->SetMaterialPropertiesTable(MPTLyso);

  matLeadGlass = new G4Material("LeadGlass",6.22*g/cm3,5);
  matLeadGlass->AddElement(elPb,75.1938*perCent);
  matLeadGlass->AddElement(elO,15.6453*perCent);
  matLeadGlass->AddElement(elSi,8.0866*perCent);
  matLeadGlass->AddElement(elTi,0.8092*perCent);
  matLeadGlass->AddElement(elAs,0.2651*perCent);
  matLeadGlass->GetIonisation()->SetMeanExcitationEnergy(526.4*eV);

  matAerog = new G4Material("Aerogel",0.200*g/cm3,3);
  matAerog->AddMaterial(matSiO2,62.5*perCent);
  matAerog->AddMaterial(matH2O ,37.4*perCent);
  matAerog->AddElement (elC   , 0.1*perCent);

  matPoron = new G4Material("Poron",0.320*g/cm3,4);
  matPoron->AddElement(elC,20*perCent);
  matPoron->AddElement(elO,40*perCent);
  matPoron->AddElement(elH,20*perCent);
  matPoron->AddElement(elN,20*perCent);

  matCO2 =  new G4Material("CarbonicGas",27.*mg/cm3,2,kStateGas,325.*kelvin,50.*atmosphere);
  matCO2->AddElement(elC,1);
  matCO2->AddElement(elO,2);

  matKapton = new G4Material("Kapton", 1.39*g/cm3,3);
  matKapton->AddElement(elO,2);
  matKapton->AddElement(elC,5);
  matKapton->AddElement(elH,4);

  matDacron = new G4Material("Dacron",1.37*g/cm3,3);
  matDacron->AddElement(elC,10);
  matDacron->AddElement(elH,8);
  matDacron->AddElement(elO,4);

  matEpoxy = new G4Material("Epoxy",1.2*g/cm3,2);
  matEpoxy->AddElement(elH,2);
  matEpoxy->AddElement(elC,2);

  matFR4 = new G4Material("FR4"  , 1.86*g/cm3, 2);
  matFR4->AddMaterial(matSiO2,0.528);
  matFR4->AddMaterial(matEpoxy,0.472);

  matHCAlu = new G4Material("HCAluminium", 1.35*g/cm3, 2);
  matHCAlu->AddMaterial(matAl,0.5);
  matHCAlu->AddMaterial(matAir,0.5);


}

G4Material* HEPDSWMaterial::GetMaterial(G4String material)
{
  G4Material* pttoMaterial = G4Material::GetMaterial(material); 
  return pttoMaterial; 
}
