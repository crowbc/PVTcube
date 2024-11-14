// user defined header file for class
#include "PVTcubeEvent.hh"
// Constructor
PVTcubeEventAction::PVTcubeEventAction(PVTcubeRunAction* aRun) : G4UserEventAction()
{
	// initialize with 0 energy
	fEdepPVTcube = 0.0;
}
// Destructor
PVTcubeEventAction::~PVTcubeEventAction()
{}
// Begin of Event Action
void PVTcubeEventAction::BeginOfEventAction(const G4Event* anEvent)
{
	// Set energy to 0 at beginning of every event
	fEdepPVTcube = 0.0;
	// Get event number, print event number for every 1000th event
	fEvent = anEvent->GetEventID();
	//G4int rNum = G4RunManager::GetRunManager()->GetCurrentEvent()->GetRunID();
	if(fEvent%1000 == 0)
	{
		if (fEvent == 0) G4cout << "Beginning of run ..." << G4endl << G4endl;
		G4cout << "Beginning of event # " << fEvent << "..." << G4endl;
	}
}
// End of Event Action
void PVTcubeEventAction::EndOfEventAction(const G4Event* anEvent)
{
	// Initialize analysis manager and fill N tuple with energy depositions
	G4AnalysisManager *Aman = G4AnalysisManager::Instance();
	// Fill scoring N-tuple
	Aman->FillNtupleIColumn(2, 0, fEvent);
	// finish Scoring
	Aman->AddNtupleRow(2);
	// TODO: add energy scoring
}
