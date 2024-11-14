// included user defined header file
#include "PVTcubeGenerator.hh"
// Constructor
PVTcubePrimaryGenerator::PVTcubePrimaryGenerator()
{
	fParticleSource = new G4GeneralParticleSource();
}
// Destructor
PVTcubePrimaryGenerator::~PVTcubePrimaryGenerator()
{
	delete fParticleSource;
}
// GeneratePrimaries()
void PVTcubePrimaryGenerator::GeneratePrimaries(G4Event *PVTcubeEvent)
{
	// generate vertex
	fParticleSource->GeneratePrimaryVertex(PVTcubeEvent);
}
