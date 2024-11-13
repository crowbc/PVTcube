// included user defined header file
#include "PVTcubeGenerator.hh"
// Constructor
PVTcubePrimaryGenerator::PVTcubePrimaryGenerator()
{
	//fParticleGun = new G4ParticleGun(1);
	fParticleSource = new G4GeneralParticleSource();
	// define the particle
	/*G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName = "geantino";//"ion";"opticalphoton";"e-";"e+";"gamma";
	G4ParticleDefinition *particle = particleTable->FindParticle(particleName);
	// set ion properties for Na-22 decay source
	Z = 11;
	A = 22;
	charge = 0.*eplus;
	energy = 0.*MeV;
	// Note: Activity is assumed to be 2.7897E05 Bq, so runs should be multiples of that to simulate 1s. 
	// Place source on -z face of NuLat detector.
	G4ThreeVector pos(0., 0.25*m + 1.5*in, -15.0*cm);
	// Set initial momentum on random cone into the detector from approximately 20 cm away, centered on the -z face
	// note: to be truly uniform in theta, must use inverse cosine of random number between 1 and cos(cone_angle) -- however, this is close enough for sufficiently small values of cone_angle
	theta = G4UniformRand()*cone_angle;
	phi = G4UniformRand()*360*deg;
	G4ThreeVector mom(1., 0., 0.);//mom(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
	// set particle gun properties
	G4double eCo60Lo = 1.173*MeV;// deprecated (Z = 27, A = 60 for this source)
	G4double eCo60Hi = 1.332*MeV;// deprecated (Z = 27, A = 60 for this source)
	G4double eNa22 = 1.275*MeV;// Assume annihilation gammas trigger coincidence detector - see Knoll ch 1
	G4double eminusEn = 2.0*MeV;
	G4double eplusEn = 2.0*MeV;
	fParticleGun->SetParticlePosition(pos);
	fParticleGun->SetParticleMomentumDirection(mom);
	fParticleGun->SetParticleDefinition(particle);
	if(particle==G4Geantino::Geantino())
	{
		fParticleGun->SetParticleMomentum(energy);
	}
	if(particleName=="gamma")
	{
		fParticleGun->SetParticleMomentum(eNa22);
	}
	if(particleName=="ion")
	{
		particle = G4IonTable::GetIonTable()->GetIon(Z, A, energy);
		fParticleGun->SetParticleCharge(charge);
	}
	if(particleName=="opticalphoton")
	{
		fParticleGun->SetParticleMomentum(2.95*eV);
	}
	if(particleName == "e-")
	{
		fParticleGun->SetParticleMomentum(eminusEn);
	}
	if(particleName == "e+")
	{
		fParticleGun->SetParticleMomentum(eplusEn);
	}/**/
}
// Destructor
PVTcubePrimaryGenerator::~PVTcubePrimaryGenerator()
{
	//delete fParticleGun;
	delete fParticleSource;
}
// GeneratePrimaries()
void PVTcubePrimaryGenerator::GeneratePrimaries(G4Event *PVTcubeEvent)
{
	// generate vertex
	//fParticleGun->GeneratePrimaryVertex(PVTcubeEvent);
	fParticleSource->GeneratePrimaryVertex(PVTcubeEvent);
}
