////////////////////////////////////////////////////////////////////////////////
//
#include "Vertex.hh"

G4Allocator<Vertex> VertexAllocator;
////////////////////////////////////////////////////////////////////////////////
//
Vertex::Vertex ():thePosition(0)
{
  theVolumeName="";
  theInelastic=false;
  theQuasielastic=false;
}

////////////////////////////////////////////////////////////////////////////////
//
Vertex::Vertex (G4bool aQuasielastic,G4bool aInelastic, G4String aVolumeName,
		G4ThreeVector aPosition)
{ 
  theVolumeName          = aVolumeName;
  thePosition            = aPosition;	 
  theQuasielastic        = aQuasielastic;
  theInelastic           = aInelastic;
}
////////////////////////////////////////////////////////////////////////////////
//
Vertex::~Vertex ()
{}
////////////////////////////////////////////////////////////////////////////////
//
Vertex::Vertex (const Vertex& right) : G4VHit()
{
  theVolumeName          = right.theVolumeName;
  thePosition            = right.thePosition;	 
  theQuasielastic        = right.theQuasielastic;
  theInelastic           = right.theInelastic;
}
////////////////////////////////////////////////////////////////////////////////
//
const Vertex& Vertex::operator= (const Vertex& right)
{ 
  theVolumeName          = right.theVolumeName;
  thePosition            = right.thePosition;	 
  theQuasielastic        = right.theQuasielastic;
  theInelastic           = right.theInelastic;
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
//
int Vertex::operator== (const Vertex& ) const
{
  return 0;
}

