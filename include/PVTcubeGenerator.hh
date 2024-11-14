#ifndef PVTCUBEGENERATOR_HH
#define PVTCUBEGENERATOR_HH
// included Geant4 header files
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4RunManager.hh"
// class definition
class PVTcubePrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
	PVTcubePrimaryGenerator();
	~PVTcubePrimaryGenerator();
	
	virtual void GeneratePrimaries(G4Event* anEvent);
	
private:
	G4GeneralParticleSource* fParticleSource;
	// useful constants
};
#endif
