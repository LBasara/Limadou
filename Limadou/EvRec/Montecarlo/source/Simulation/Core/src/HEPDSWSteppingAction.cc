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
//
//

#include "HEPDSWSteppingAction.hh"
#include "HEPDSWPrimaryGeneratorAction.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TrackStatus.hh"

HEPDSWSteppingAction::HEPDSWSteppingAction(HEPDSWPrimaryGeneratorAction*)
{

  Xpos = -999.999;
  Ypos = -999.999;
  Zpos = -999.999;
  TrackID = -1;
  stepNumber=0;

}

HEPDSWSteppingAction::~HEPDSWSteppingAction()
{
}

void HEPDSWSteppingAction::UserSteppingAction(const G4Step* )
{
  /*
void HEPDSWSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  if(fabs((aStep->GetPreStepPoint()->GetPosition()).getX()-Xpos)<0.01*mm && 
     fabs((aStep->GetPreStepPoint()->GetPosition()).getY()-Ypos)<0.01*mm && 
     fabs((aStep->GetPreStepPoint()->GetPosition()).getZ()-Zpos)<0.01*mm && 
     aStep->GetTrack()->GetTrackID()==TrackID){
      stepNumber++;
      if(stepNumber>50){
	Xpos = -999.999;
	Ypos = -999.999;
	Zpos = -999.999;
	TrackID = -1;
	stepNumber=0;
	aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
	G4Cout<<"WARNING EVENT WITH TRACK ABORTED BY OUR XCHECK"<<G4Endl;
      }
    }else{
      Xpos=(aStep->GetPreStepPoint()->GetPosition()).getX();
      Ypos=(aStep->GetPreStepPoint()->GetPosition()).getY();
      Zpos=(aStep->GetPreStepPoint()->GetPosition()).getZ();
      TrackID=aStep->GetTrack()->GetTrackID();
      stepNumber=0;
    }
  */
}
 
