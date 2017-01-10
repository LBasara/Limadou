#ifndef HEPDSWTrackingAction_h
#define HEPDSWTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "HEPDSWTrackingActionMessenger.hh"
#include "globals.hh"

#include <vector>
#include <map>

#include "G4PrimaryParticle.hh"
#include "G4Track.hh"

class HEPDSWTrackingActionMessenger;

////////////////////////////////////////////////////////////////////////////////
//
class HEPDSWTrackingAction : public G4UserTrackingAction
{
public:
  HEPDSWTrackingAction ();
  ~HEPDSWTrackingAction ();

  void PreUserTrackingAction (const G4Track* theTrack);
  void PostUserTrackingAction (const G4Track* theTrack);

  inline void SetSaveTrack(G4bool aVal){saveTrack=aVal;}
  inline void SetSavePrimaryParticleVertex(G4bool aVal){savePrimaryParticleVertex=aVal;}

private:
  HEPDSWTrackingActionMessenger* fTrackingActionMessenger;
  G4bool saveTrack;
  G4bool savePrimaryParticleVertex;

  std::map<const G4Track*, int> _tracks; 

};
////////////////////////////////////////////////////////////////////////////////
#endif
