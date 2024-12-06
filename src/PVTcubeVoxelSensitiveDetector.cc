// included user defined header
#include "PVTcubeVoxelSensitiveDetector.hh"
// Constructor
PVTcubeVoxelSensitiveDetector::PVTcubeVoxelSensitiveDetector(G4String name, G4int nCubes) : G4VSensitiveDetector(name)
{
	nVox = nCubes;
}
// Destructor
PVTcubeVoxelSensitiveDetector::~PVTcubeVoxelSensitiveDetector()
{}
// Process Hits in the Sensitive Detector
G4bool PVTcubeVoxelSensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROHist)
{
	G4int evt, pID, trkID;
	G4double eDepTot = aStep->GetTotalEnergyDeposit();
	G4double xPos, yPos, zPos, time, pX0, pY0, pZ0, pMag;
	G4Track *track = aStep->GetTrack();
	trkID = track->GetTrackID();
	G4String pName = track->GetDefinition()->GetParticleName();
	pID = ParticleNameToIDNumber(pName);
	//G4ThreeVector posHit, momHit;
	if (eDepTot == 0.)
	{
		return true;
	}
	// don't count hits from optical photons
	if (pID == 100)
	{
		return true;
	}
	const G4VTouchable *touch = aStep->GetPreStepPoint()->GetTouchable();
	G4VPhysicalVolume *physVol = touch->GetVolume();
	G4int copyNo = physVol->GetCopyNo();
	// populate variables
	evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	//posHit = aStep->GetTrack()->GetPosition();
	xPos = aStep->GetTrack()->GetPosition().x();
	yPos = aStep->GetTrack()->GetPosition().y();
	zPos = aStep->GetTrack()->GetPosition().z();
	time = aStep->GetTrack()->GetGlobalTime();
	pX0 = aStep->GetTrack()->GetVertexMomentumDirection().x();
	pY0 = aStep->GetTrack()->GetVertexMomentumDirection().y();
	pZ0 = aStep->GetTrack()->GetVertexMomentumDirection().z();
	pMag = aStep->GetTrack()->GetMomentum().mag();
	// initialize analysis manager and fill Ntuples
	G4AnalysisManager *aMan = G4AnalysisManager::Instance();
	// Fill Ntuple columns with Voxel energy deposition information
	aMan->FillNtupleIColumn(2, 0, evt);
	aMan->FillNtupleIColumn(2, 1, pID);
	aMan->FillNtupleIColumn(2, 2, trkID);
	aMan->FillNtupleIColumn(2, 3, copyNo);
	aMan->FillNtupleDColumn(2, 4, eDepTot);
	aMan->FillNtupleDColumn(2, 5, xPos);
	aMan->FillNtupleDColumn(2, 6, yPos);
	aMan->FillNtupleDColumn(2, 7, zPos);
	aMan->FillNtupleDColumn(2, 8, time);
	aMan->FillNtupleDColumn(2, 9, pX0);
	aMan->FillNtupleDColumn(2, 10, pY0);
	aMan->FillNtupleDColumn(2, 11, pZ0);
	aMan->FillNtupleDColumn(2, 12, pMag);
	aMan->AddNtupleRow(2);
	// Return value
	return true;
}
/* ---------------------------------------------- */
/*  Convert a steps particle name to an ID number */
/*  specific to the PVTcube analysis              */
/* ---------------------------------------------- */
G4int PVTcubeVoxelSensitiveDetector::ParticleNameToIDNumber(G4String name)
{
	G4int num;
	if(name == "gamma"){
		num = 1;
	}
	else if(name == "e"){
		num = 2;
	}
	else if(name == "e+"){
		num = 3;
	}
	else if(name == "neutron"){
		num = 4;
	}
	else if(name == "proton"){
		num = 5;
	}
	else if(name == "mu+"){
		num = 6;
	}
	else if(name == "mu-"){
		num = 7;
	}
	else if(name == "alpha"){
		num = 8;
	}
	else if(name == "Li7"){
		num = 9;
	}
	else if(name == "opticalphoton"){
		num = 100;
	}
	else{
		num = 0;
	}
	return num;
}
