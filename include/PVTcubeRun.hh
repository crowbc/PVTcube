// start conditional to define class only once
#ifndef PVTCUBERUN_HH
#define PVTCUBERUN_HH
// Geant4 header files
#include "G4UserRunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
// C++ header files
#include <vector>
using namespace std;
// define the class
class PVTcubeRunAction : public G4UserRunAction
{
public:
	// Constructor and Destructor
	PVTcubeRunAction();
	~PVTcubeRunAction();
	// User Run actions
	virtual void BeginOfRunAction(const G4Run*);
	virtual void EndOfRunAction(const G4Run*);
private:
	G4LogicalVolume* scoringVolume;
};
// end of conditional to define class only once
#endif
