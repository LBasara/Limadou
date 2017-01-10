#include "HEPDSWTrackingAction.hh"
#include "G4ThreeVector.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "Track.hh"
#include "Vertex.hh"
#include "MCTruthSD.hh"

////////////////////////////////////////////////////////////////////////////////
//
HEPDSWTrackingAction::HEPDSWTrackingAction ()
{
  fTrackingActionMessenger = new HEPDSWTrackingActionMessenger(this);
  saveTrack=true;
  savePrimaryParticleVertex=true;
}
////////////////////////////////////////////////////////////////////////////////
//
HEPDSWTrackingAction::~HEPDSWTrackingAction ()
{delete fTrackingActionMessenger;}
////////////////////////////////////////////////////////////////////////////////
//
void HEPDSWTrackingAction::PostUserTrackingAction (const G4Track* track)
{
  if (savePrimaryParticleVertex && _tracks.size() > 0) {
    // See if the track is in the tracking list
    std::map<const G4Track*, int>::iterator trackIter = _tracks.find(track);
    if (trackIter != _tracks.end()) {
      // Check process type
      const G4Step *step = track->GetStep();
      int primaryPDG = track->GetDefinition()->GetPDGEncoding();
      // Check if interaction is quasi-elastic
      float maxSecondaryEnergy = step->GetPreStepPoint()->GetKineticEnergy() * 0.9;
      const G4TrackVector *secondaries = step->GetSecondary();
      unsigned int nSecondaries = secondaries->size();
      bool inelastic = true;
      bool quasielastic = false;
      G4ThreeVector vertexPosition = step->GetPreStepPoint()->GetPosition();
      G4String vertexVolume = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
      for (unsigned int is = 0; is < nSecondaries; is++) {
	G4Track *currentSecondaryTrack = (*secondaries)[is];
	if (currentSecondaryTrack->GetDefinition()->GetPDGEncoding() == primaryPDG
	    && currentSecondaryTrack->GetKineticEnergy() > maxSecondaryEnergy) {
	  _tracks[currentSecondaryTrack] = trackIter->second;
	  inelastic = false;
	  quasielastic = true;
	  break;
	}
      }
      Vertex* storeVertex = new Vertex(quasielastic,inelastic,vertexVolume,vertexPosition);
      MCTruthSD* sd = dynamic_cast<MCTruthSD*>(G4SDManager::GetSDMpointer()->FindSensitiveDetector("/hepd/mctruth"));
      sd->GetVertexCollection()->insert(storeVertex);
      _tracks.erase(trackIter);
      
    }
  } 
}

////////////////////////////////////////////////////////////////////////////////
//
void HEPDSWTrackingAction::PreUserTrackingAction (const G4Track* aTrack)
{ 
  // Check whether it is a primary particle
  G4PrimaryParticle *primary = aTrack->GetDynamicParticle()->GetPrimaryParticle();
  if (primary != NULL)
    _tracks[aTrack] = aTrack->GetTrackID();
  
  if(saveTrack){
    G4int         TrackId = aTrack->GetTrackID();
    G4int         PDG = aTrack->GetDefinition()->GetPDGEncoding();
    G4ThreeVector position = aTrack->GetVertexPosition();
    G4ThreeVector dir = aTrack->GetVertexMomentumDirection();
    G4double      kinE = aTrack->GetVertexKineticEnergy();
    G4int         motherTrackId = aTrack->GetParentID();
    G4String      vertexVolumeName =  aTrack->GetVolume()->GetName();
    G4String      creatorProcessName = "unknown";
    if (aTrack->GetCreatorProcess() != 0){
      creatorProcessName = aTrack->GetCreatorProcess()->GetProcessName();
    }
    Track* storeTrack = new Track(TrackId,PDG,position,dir,kinE,motherTrackId,vertexVolumeName,creatorProcessName);
    
    MCTruthSD* sd = dynamic_cast<MCTruthSD*>(G4SDManager::GetSDMpointer()->FindSensitiveDetector("/hepd/mctruth"));
    sd->GetTrackCollection()->insert(storeTrack);
  }
  
}
