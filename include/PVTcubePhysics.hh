#ifndef PVTCUBEPHYSICS_HH
#define PVTCUBEPHYSICS_HH
// include Geant4 Physics packages
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
// write the class
class PVTcubePhysicsList : public G4VModularPhysicsList
{
public:
	// Constructor
	PVTcubePhysicsList();
	// Destructor
	~PVTcubePhysicsList();
};
#endif
