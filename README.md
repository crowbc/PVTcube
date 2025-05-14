# PVTcube

This is a Geant4 simulation of a single polyvinyltoulene (PVT) cube in a simple environment where backgrounds may be present. This simulation is written for Geant4 v11 and is not backwards compatible with older Geant4 builds than v 11.0.

# build
Initial build for testing geometry and materials properties, physics of detector medium and supporting components, and PMT's. Build uses v 11.1.0 of Geant4.

Change log (in reverse chronological order starting from most recent version:)

14MAY2025 -- v 1.5.0 Fourfold Coincidence
	Added a fourth detector voxel, foil wrapping, photomultiplier tube, lens and mu metal shield to simulate fourfold coincidence with a secondary scattering target placed near the third cube. Added analysis scripts and ntuple merging scripts which are designed to be used on the output csv files in the directory where these files are located. Corrected geometry for cube dimensions.

12DEC2024 -- v 1.4.1 Updated Macros
	Updated batch mode macros with event count that more accurately reflects the age of the source being simulated. Changed the Event Manager to update every 100,000th event instead of every 1,000th event.
	
05DEC2024 -- v 1.4.0 Threefold Coincidence
	Added a third detector voxel, foil wrapping, photomultiplier tube, lens and mu metal shield as well as lead shielding blocks to simulate threefold coincidence with a collumation channel between the source and the third voxel. Changed default output ntuple format to csv.

22NOV2024 -- v 1.3.1 Unused Volume Fix
	Removed unused volume declarations in detector construction class. Fixed CMakeLists.txt to copy macro files to build directory.
	
20NOV2024 -- v 1.3.0 Scoring Volume Update
	Added voxel sensitive detector class and registered voxels as scoring volumes in SD manager. Implemented scoring N-Tuples in sensitive detector class.

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
