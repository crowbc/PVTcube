// conditional statement for defining class only once
#ifndef PVTCUBEDETECTORCONSTRUCTION_HH
#define PVTCUBEDETECTORCONSTRUCTION_HH
// Header files for defining volumes
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
// Nist manager and units for material properties
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
// Header files for placements
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
// Header files for geometry types
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
// Header files for sensitive detector
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
// Header files for surfaces
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
// Header file for messenger control
#include "G4GenericMessenger.hh"
// Header files for visual attribute manager
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
// Header file for user defined libraries
#include "PVTcubePMTsensitiveDetector.hh"
#include "PVTcubeVoxelSensitiveDetector.hh"
// Write the class
class PVTcubeDetectorConstruction : public G4VUserDetectorConstruction
{
public:
	PVTcubeDetectorConstruction();
	~PVTcubeDetectorConstruction();
	// method for looking up scoring volume
	G4LogicalVolume *GetPVTcubeScoringVolume() const { return fPVTcubeScoringVolume; }
	// construct function for detector factory
	virtual G4VPhysicalVolume* Construct();
private:
	// Volume declarations - naming convention: solidName for geometry volume definitions, logicName for logical volume definitions and physName for physical volume definitions
	G4Box *solidWorld, *solidVoxel, *solidTableTop, *solidTable, *solidNotePad, *solidDarkBoxOuter, *solidDarkBoxInner, *solidFoilWrapInner, *solidFoilWrapOuter, *solidLeadBlock;
	G4Trd *solidLGTrd;
	G4Tubs *solidPMT, *solidPMTGlass, *solidSource, *solidPMTshield;
	G4Cons *solidLGCone;
	G4Sphere *solidPMTconvex;
	G4SubtractionSolid *solidPMTLens, *solidDarkBox, *solidPMTshieldwall, *solidFoilWrap;
	G4IntersectionSolid *solidLG;
	G4LogicalVolume *logicWorld, *logicVoxel, *fPVTcubeScoringVolume, *logicLG, *logicPMT, *logicPMTLens, *logicTableTop, *logicTable, *logicSource, *logicNotePad, *logicDarkBox, *logicPMTshieldwall, *logicFoilWrap, *logicLeadBlock;
	G4VPhysicalVolume *physWorld, *physVoxel, *physLG, *physPMT, *physPMTLens, *physTableTop, *physTable, *physSource, *physNotePad, *physDarkBox, *physPMTshieldwall, *physFoilWrap, *physLeadBlock;
	// Declare optical surfaces
	G4OpticalSurface *surface_Al, *surface_SS;
	// Material declarations
	G4Element *H, *Be, *C, *O, *Si, *Cr, *Fe, *Ni, *Cu, *Mo, *Pb;
	G4Material *air, *EJ200, *acrylic, *aluminum, *concrete, *stainless, *lead, *borosilicateGlass, *mumetal, *paper, *plywood;
	G4MaterialPropertiesTable *mpt_air, *mpt_EJ200, *mpt_acrylic, *mpt_Al, *mpt_SS, *mpt_BorosilicateGlass;
	// useful constants:
	// physical constant for computing photon energies or wavelengths: (note - divide by wavelength in nm to get energy in eV, or divide by energy in eV to get wavelength in nm)
	const G4double HCNM = 1239.841939*eV*nm;
	// conversion factor inches to mm
	const G4double in = 25.4*mm;
	// amu in SI:
	const G4double uu = 1.66054E-27*kg;
	// mass of H and C atoms in amu
	const G4double Cmass = 12.011*uu;
	const G4double Hmass = 1.0079*uu;
	// Constants for material properties:
	// air
	const G4double rindexConst_air = 1.0;
	// EJ200 - scintillation yield, decay time, rise time, peak emission wavelength (in nm) and full width half maximum:
	const G4double scyld_EJ200 = 10000./MeV;
	const G4double dt_EJ200 = 2.1*ns;
	const G4double rt_EJ200 = 0.9*ns;
	const G4double wlmax_EJ200 = 425*nm;
	const G4double fwhm_EJ200 = 2.5*ns;
	// EJ200 light attenuation length, refractive index, C atom number density, H atom number density and electon number density:
	const G4double aLenConst_EJ200 = 100.*cm;// note in Yokley (2106) -- best fit attenuation length is 93 cm (p. 131) or 113 cm (p. 134); EJ-200 spec sheet lists 380 cm
	const G4double rindexConst_EJ200 = 1.58;
	const G4double nH_EJ200 = 5.17E22/cm3;
	const G4double nC_EJ200 = 4.69E22/cm3;
	const G4double ne_EJ200 = 3.33E23/cm3;
	const G4double rho_EJ200 = 1.023*g/cm3;
	// acrylic properties
	const G4double aLenConst_acrylic = 10.0*m;
	const G4double rindexConst_acrylic = 1.492;
	const G4double rho_acrylic = 1.180*g/cm3;
	// Reflective surface properties
	const G4double refConst_Al = 0.95;
	const G4double refConst_SS = 1.0;
	// Borosilicate Glass properties
	const G4double rho_BorosilicateGlass = 2.65*g/cm3;
	const G4double rindexConst_BorosilicateGlass = 1.55;
	const G4double aLenConst_BorosilicateGlass = 100.0*cm;
	// mu metal properties
	const G4double rho_mumetal = 8.7*g/cm3;
	// variable declarations
	// World volume size in x, y and z dimensions
	G4double xWorld, yWorld, zWorld;
	// Table dimensions
	G4double w_table, l_table, h_table, t_table;
	// Detector Parameters
	G4int nCubes;
	G4double xVoxelSize, yVoxelSize, zVoxelSize;
	G4double xVoxelSep, yVoxelSep, zVoxelSep;
	G4double lenPMT, lenLGBase, lenLGTaper, lenLGSqu, lenLGwPMT, r1_LG, r2_LG;
	G4double rPMT, tGlass, rSpherePMTsurf, tGlass_min;
	G4double x_block, y_block, z_block;
	// Other geometry parameters
	G4double w_notepad, l_notepad, h_notepad, r_source_inner, r_source_outer, t_source, t_foil;
	// Mass fraction of C and H in EJ200
	G4double Hmass_EJ200 = Hmass*nH_EJ200;
	G4double Cmass_EJ200 = Cmass*nC_EJ200;
	// mass per cc
	G4double totalMass_EJ200 = Cmass_EJ200 + Hmass_EJ200;
	G4double Cfrac_EJ200 = nC_EJ200*Cmass/totalMass_EJ200;
	G4double Hfrac_EJ200 = nH_EJ200*Hmass/totalMass_EJ200;
	// function declarations
	void DefineMaterials();
	virtual void ConstructSDandField();
	// Pointers for SD's:
	PVTcubePMTSensitiveDetector *detPMT;
	PVTcubeVoxelSensitiveDetector *detVoxel;
	// Pointer for primitive scorer
	G4VPrimitiveScorer *primVoxel;
	// Pointer to Generic Messenger Object
	G4GenericMessenger *fMessenger;
	// Pointer to Visual Attribute manager object
	G4VisAttributes *attr;
};
// end of conditional for defining class only once
#endif
