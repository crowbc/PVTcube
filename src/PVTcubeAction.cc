#include "PVTcubeAction.hh"
// constructor
PVTcubeActionInitialization::PVTcubeActionInitialization()
{}
// destructor
PVTcubeActionInitialization::~PVTcubeActionInitialization()
{}
// build for master thread
void PVTcubeActionInitialization::BuildForMaster() const
{
	// Do only Run Action for Master Thread
	PVTcubeRunAction *runAction = new PVTcubeRunAction();
	SetUserAction(runAction);/**/
}
// build function
void PVTcubeActionInitialization::Build() const
{
	// Generator Action (uncomment this first)
	PVTcubePrimaryGenerator *generator = new PVTcubePrimaryGenerator();
	SetUserAction(generator);
	// Run Action
	PVTcubeRunAction *runAction = new PVTcubeRunAction();
	SetUserAction(runAction);
	// Event Action
	PVTcubeEventAction *eventAction = new PVTcubeEventAction(runAction);
	SetUserAction(eventAction);
	// Stepping Action
	/*PVTcubeSteppingAction *steppingAction = new PVTcubeSteppingAction(eventAction);
	SetUserAction(steppingAction);/**/
}
