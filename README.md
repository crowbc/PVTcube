# PVTcube

This is a Geant4 simulation of a single polyvinyltoulene (PVT) cube in a simple environment where backgrounds may be present. This simulation is written for Geant4 v11 and is not backwards compatible with older Geant4 builds than v11.0.

# build
Initial build for testing geometry and materials properties, physics of detector medium and supporting components, and PMT's. Build uses v 11.1.0 of Geant4.

Change log (in reverse chronological order starting from most recent version:)

v 1.0.0 -- 12NOV2024
	Initial Commit, created main function for initializing UI or initializing simulation in batch mode. Created particle generator, physics manager, construction manager, action manager, run manager, and event manager classes. Created CMakelists file for compiling. Created a macro for batch mode runs and a macro for interactive mode runs with a radioactive source. 