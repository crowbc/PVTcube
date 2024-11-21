// user defined header file for class
#include "PVTcubeRun.hh"
// Constructor
PVTcubeRunAction::PVTcubeRunAction()
{
	// initialize analysis manager
	G4AnalysisManager *Aman = G4AnalysisManager::Instance();
	// create N tuple for photon truth information in NuLat PMT hits
	Aman->CreateNtuple("PVTcube_PMT_Truth", "PVTcube_PMT_Truth");// Ntuple 0
	Aman->CreateNtupleIColumn("fEvent");// col 0
	Aman->CreateNtupleDColumn("fX");// col 1
	Aman->CreateNtupleDColumn("fY");// col 2
	Aman->CreateNtupleDColumn("fZ");// col 3
	Aman->CreateNtupleDColumn("fWlen");// col 4
	Aman->FinishNtuple(0);
	// create N tuple for PMT hit information
	Aman->CreateNtuple("PVTcube_PMT_Hits", "PVTcube_PMT_Hits");// Ntuple 1
	Aman->CreateNtupleIColumn("fEvent");// col 0
	// PMT ID
	Aman->CreateNtupleIColumn("fID");// col 1
	// PMT coords
	Aman->CreateNtupleDColumn("fX");// col 2
	Aman->CreateNtupleDColumn("fY");// col 3
	Aman->CreateNtupleDColumn("fZ");// col 4
	// time of hit
	Aman->CreateNtupleDColumn("fT");// col 5
	Aman->FinishNtuple(1);
	// create N tuple for energy deposition tracking (MC truth)
	Aman->CreateNtuple("PVTcube_Scoring", "PVTcube_Scoring");// Ntuple 2
	Aman->CreateNtupleIColumn("fEvent");// col 0
	// columns for particle ID and track ID
	Aman->CreateNtupleIColumn("fPID");// col 1
	Aman->CreateNtupleIColumn("fTrkID");// col 2
	// Voxel ID - can use to validate whether or not coordinates are true
	Aman->CreateNtupleIColumn("fID");// col 3
	// scoring value for energy deposition
	Aman->CreateNtupleDColumn("fTotEdep");// col 4
	// Voxel coords
	Aman->CreateNtupleDColumn("fX");// col 5
	Aman->CreateNtupleDColumn("fY");// col 6
	Aman->CreateNtupleDColumn("fZ");// col 7
	// time of hit
	Aman->CreateNtupleDColumn("fT");// col 8
	// columns for momentum components
	Aman->CreateNtupleDColumn("fPX0");// col 9
	Aman->CreateNtupleDColumn("fPY0");// col 10
	Aman->CreateNtupleDColumn("fPZ0");// col 11
	Aman->CreateNtupleDColumn("fPmag");// col 12
	Aman->FinishNtuple(2);
	// TODO: get trajectory info stored in Ntuple
}
// Destructor
PVTcubeRunAction::~PVTcubeRunAction()
{}
// Beginning of Run
void PVTcubeRunAction::BeginOfRunAction(const G4Run* PVTcubeRun)
{
	scoringVolume = G4LogicalVolumeStore::GetInstance()->GetVolume("logicVoxel");
	G4AnalysisManager *Aman = G4AnalysisManager::Instance();
	G4int rNum = PVTcubeRun->GetRunID();
	std::stringstream sRunID;
	sRunID << rNum;
	G4String name = "PVTcube_run";
	G4String ext = "_output.root";
	G4String fName = name + sRunID.str() + ext;
	Aman->OpenFile(fName);
	Aman->SetVerboseLevel(1);
}
// End of Run
void PVTcubeRunAction::EndOfRunAction(const G4Run* PVTcubeRun)
{
	G4AnalysisManager *Aman = G4AnalysisManager::Instance();
	// write and close the ROOT file !IMPORTANT!
	Aman->Write();
	Aman->CloseFile();
}
