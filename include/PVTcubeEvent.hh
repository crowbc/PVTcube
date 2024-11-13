// conditional to define class only once
#ifndef PVTCUBEEVENT_HH
#define PVTCUBEEVENT_HH
// Geant4 header files
#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4AnalysisManager.hh"
#include "G4TrackingManager.hh"
#include "G4VTrajectory.hh"
#include "G4RunManager.hh"
// user defined header files
#include "PVTcubeRun.hh"
// C++ header files
#include <vector>
using namespace std;
// define the class
class PVTcubeEventAction : public G4UserEventAction
{
public:
	// Constructor and Destructor
	PVTcubeEventAction(PVTcubeRunAction* aRun);
	~PVTcubeEventAction();
	// User Event Actions
	virtual void BeginOfEventAction(const G4Event* anEvent);
	virtual void EndOfEventAction(const G4Event* anEvent);
	// functions to add energy depositions
	void AddEdepPVTcube(G4double edep){ fEdepPVTcube+=edep; }
private:
	// variable to store energy depositions, hit voxel coordinates, time and momentum components
	G4double fEdepPVTcube, fX, fY, fZ, fT, fPX0, fPY0, fPZ0;
	/*
		variables to store event ID, index of hit voxel, particle ID, & track ID
	*/
	G4int fEvent, fID, fPID, fTrkID;
};
// end of conditional to define class only once
#endif
