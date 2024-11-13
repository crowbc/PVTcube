// User-Defined Header
#include "PVTcubePhysics.hh"
// Constructor
PVTcubePhysicsList::PVTcubePhysicsList()
{
	RegisterPhysics(new G4EmStandardPhysics());
	RegisterPhysics(new G4OpticalPhysics());
	RegisterPhysics(new G4DecayPhysics());
	RegisterPhysics(new G4RadioactiveDecayPhysics());
}
// Destructor
PVTcubePhysicsList::~PVTcubePhysicsList()
{}
