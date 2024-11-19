# PVTcube

This is a Geant4 simulation of a single polyvinyltoulene (PVT) cube in a simple environment where backgrounds may be present. This simulation is written for Geant4 v11 and is not backwards compatible with older Geant4 builds than v 11.0.

# build
Initial build for testing geometry and materials properties, physics of detector medium and supporting components, and PMT's. Build uses v 11.1.0 of Geant4.

Change log (in reverse chronological order starting from most recent version:)

18NOV2024 -- v 1.2.1: Foil Covering Fix
	Removed the +y surface of the foil covering to unblock the PMT and light guide assembly from optical contact with the PVT cube. Adjusted "PMT" radius to correct value.

15NOV2024 -- v 1.2.0: Two-fold Coincidence
	Added a second detector voxel, photomultiplier tube, lens and mu metal shield to simulate two-fold coincidence. Added aluminum foil wrapping to PVT cubes to simulate actual set-up and to eliminate optical cross-talk.

14NOV2024 -- v 1.1.1: Physical Volume Overlaps
	Corrected physical volume overlaps in detector construction class.

14NOV2024 -- v 1.1.0: Sensitive Detector Manager
	Created sensitive detector manager for recording phototube hits using a simple random number generator to simulate PMT quantum efficiency. Added mu metal shield and PMT housing to geometry in detector construction manager class.

12NOV2024 -- v 1.0.0: Initial Commit
	Created main function for initializing UI or initializing simulation in batch mode. Created particle generator, physics manager, construction manager, action manager, run manager, and event manager classes. Created CMakelists file for compiling. Created a macro for batch mode runs and a macro for interactive mode runs with a radioactive source. 
