// conditional to define the class only once
#ifndef PVTCUBEACTION_HH
#define PVTCUBEACTION_HH
// Geant4 header file
#include "G4VUserActionInitialization.hh"
// user defined header files
#include "PVTcubeGenerator.hh"
#include "PVTcubeRun.hh"
#include "PVTcubeEvent.hh"
//#include "PVTcubeStepping.hh"
// define the class
class PVTcubeActionInitialization : public G4VUserActionInitialization
{
public:
	// Constructor and Destructor
	PVTcubeActionInitialization();
	~PVTcubeActionInitialization();
	// Build User Action Methods
	virtual void Build() const;
	virtual void BuildForMaster() const;
};
// end of conditional to define the class only once
#endif
