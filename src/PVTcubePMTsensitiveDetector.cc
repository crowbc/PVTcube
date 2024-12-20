// Header for user defined library
#include "PVTcubePMTsensitiveDetector.hh"
// Constructor
PVTcubePMTSensitiveDetector::PVTcubePMTSensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{
	queff = 0.25;// use uniform value of 25%
}
// Destructor
PVTcubePMTSensitiveDetector::~PVTcubePMTSensitiveDetector()
{}
// Process Hits Function
G4bool PVTcubePMTSensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROHist)
{
	//debugMsg = false;// set to true for debug messages
	G4Track *track = aStep->GetTrack();
	pName = track->GetDefinition()->GetParticleName();
	pID = ParticleNameToIDNumber(pName);
	// only do this for optical photons to simulate the regime where the PMT is actually sensitve
	if (pID == 100)
	{
		track->SetTrackStatus(fStopAndKill);
		// set up step points. Get true photon position of PMT hit
		G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
		G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
		// Get Photon Position, Momentum and Wavelength (MC truth)
		posPhoton = preStepPoint->GetPosition();
		xpos = posPhoton[0];
		ypos = posPhoton[1];
		zpos = posPhoton[2];
		momPhoton = preStepPoint->GetMomentum();
		pX0 = momPhoton[0];
		pY0 = momPhoton[1];
		pZ0 = momPhoton[2];
		// Calculates wavelength in nm (change constant to change unit order of magnitude)
		wlen = HCNM/momPhoton.mag();
		// Debug Print if enabled -- TODO: find a more efficient way to do this rather than using conditionals
		/*if(debugMsg){
			G4cout << "Photon position: " << posPhoton << "; Photon wavelength: " << wlen << " nm" << G4endl;
		}/**/
		// Get index of hit PMT
		const G4VTouchable *touch = aStep->GetPreStepPoint()->GetTouchable();
		fID = touch->GetCopyNumber();
		// Debug Print if enabled
		/*if(debugMsg)
		{
			G4cout << "Copy number: " << fID << G4endl;
		}/**/
		// Get position of hit PMT (not photon position)
		G4VPhysicalVolume *physVol = touch->GetVolume();
		posDet = physVol->GetTranslation();
		fX = posDet[0];
		fY = posDet[1];
		fZ = posDet[2];
		fT = aStep->GetTrack()->GetGlobalTime();
		// Debug Print if enabled
		/*if(debugMsg){
			G4cout << "Detector position: " << posDet << G4endl;
		}/**/
		// initialize analysis manager and fill Ntuples
		G4AnalysisManager* Aman = G4AnalysisManager::Instance();
		// Get event number
		fEvt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
		// Fill Ntuple columns with photon MC Truth information
		Aman->FillNtupleIColumn(0, 0, fEvt);
		Aman->FillNtupleDColumn(0, 1, xpos);
		Aman->FillNtupleDColumn(0, 2, ypos);
		Aman->FillNtupleDColumn(0, 3, zpos);
		Aman->FillNtupleDColumn(0, 4, wlen);
		Aman->AddNtupleRow(0);
		// Fill Ntuple columns with sensitive detector hit information - ONLY if RNG generates a number below the quantum efficiency
		if(G4UniformRand() < queff)//quEff->Value(wlen))
		{
			Aman->FillNtupleIColumn(1, 0, fEvt);
			Aman->FillNtupleIColumn(1, 1, fID);
			Aman->FillNtupleDColumn(1, 2, fX);
			Aman->FillNtupleDColumn(1, 3, fY);
			Aman->FillNtupleDColumn(1, 4, fZ);
			Aman->FillNtupleDColumn(1, 5, fT);
			Aman->AddNtupleRow(1);
			/*if(debugMsg){
				G4cout << "PMT Hit registered for PMT # " << fID << G4endl;
			}/**/
		}
	}
	// return value
	return true;
}
// Particle Name to ID number converter
G4int PVTcubePMTSensitiveDetector::ParticleNameToIDNumber(G4String name)
{
	G4int num;
	if(name == "gamma"){
		num=1;
	}
	else if(name == "e-"){
		num=2;
	}
	else if(name == "e+"){
		num=3;
	}
	else if(name == "neutron"){
		num=4;
	}
	else if(name == "proton"){
		num=5;
	}
	else if(name == "mu+"){
		num=6;
	}
	else if(name == "mu-"){
		num=7;
	}
	else if(name == "alpha"){
		num=8;
	}
	else if(name == "Li7"){
		num=9;
	}
	else if(name == "opticalphoton"){
		num=100;
	}
	else num=0;
	return num;
}
