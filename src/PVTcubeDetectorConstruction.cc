// Included header files written for simulation
#include "PVTcubeDetectorConstruction.hh"
// Constructor
PVTcubeDetectorConstruction::PVTcubeDetectorConstruction()
{
	// Set Parameter Default Values - TODO: set default values in the messenger
	nCubes = 3;
	// default geometry dimensions
	xVoxelSize = 2.25*in;
	yVoxelSize = 2.25*in;
	zVoxelSize = 2.25*in;
	xVoxelSep = 0.000*in;
	yVoxelSep = 0.000*in;
	zVoxelSep = 0.000*in;
	// PMT and Light Guide dimensions
	lenPMT = 200.*mm;
	lenLGBase = 46.0*mm;
	lenLGTaper = 35.0*mm;
	lenLGSqu = 5.0*mm;
	lenLGwPMT = lenPMT+lenLGTaper+lenLGSqu;
	r1_LG = 3.465*in/2;
	r2_LG = 1.811*in/2;
	tGlass = 4*mm;
	rPMT = 23*mm;
	rSpherePMTsurf = 529.25*mm;
	tGlass_min = 0.8*mm;
	// Set world volume to needed size
	xWorld = 0.6*m;
	yWorld = 1.0*m;
	zWorld = 1.0*m;
	// set table dimensions
	w_table = xWorld;
	h_table = yWorld/2;
	l_table = zWorld;
	t_table = 0.25*in;
	// foil thickness
	t_foil = 0.02*mm;
	// set note pad dimensions
	w_notepad = 1.5*in;
	l_notepad = 1.5*in;
	h_notepad = yVoxelSize/4;
	// set source "puck" dimensions
	r_source_inner = 0.25*cm;
	r_source_outer = 0.5*cm;
	t_source = 0.05*cm;
	// set shielding block dimensions
	x_block = 10.0*cm;
	y_block = 10.0*cm;
	z_block = 20.0*cm;
	// Define the messenger and declare properties - for now only number of cubes can be varied.
	fMessenger = new G4GenericMessenger(this, "/detector/", "Detector Construction");
	fMessenger->DeclareProperty("nCubes", nCubes, "Number of PVT cubes in the set-up");
	// Define Detector Construction Materials
	DefineMaterials();
}
// Destructor
PVTcubeDetectorConstruction::~PVTcubeDetectorConstruction()
{}
// DefineMaterials()
void PVTcubeDetectorConstruction::DefineMaterials()
{
	G4NistManager *nist = G4NistManager::Instance();
	// size_t variable for property array lengths
	const size_t nI = 182;
	// Photon energy range for energy dependent material responses - corresponds to a range of 220 nm to 609 nm
	G4double photonEnergy[nI] = {
		2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV, 2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV, 2.341*eV, 2.386*eV, //10
		2.433*eV, 2.481*eV, 2.487*eV, 2.496*eV, 2.506*eV, 2.516*eV, 2.524*eV, 2.531*eV, 2.539*eV, 2.547*eV, //20
		2.554*eV, 2.561*eV, 2.569*eV, 2.577*eV, 2.586*eV, 2.595*eV, 2.605*eV, 2.614*eV, 2.622*eV, 2.630*eV, //30
		2.638*eV, 2.646*eV, 2.653*eV, 2.660*eV, 2.669*eV, 2.676*eV, 2.681*eV, 2.688*eV, 2.693*eV, 2.698*eV, //40
		2.703*eV, 2.706*eV, 2.711*eV, 2.718*eV, 2.723*eV, 2.731*eV, 2.742*eV, 2.755*eV, 2.768*eV, 2.782*eV, //50
		2.793*eV, 2.803*eV, 2.811*eV, 2.819*eV, 2.829*eV, 2.837*eV, 2.845*eV, 2.853*eV, 2.860*eV, 2.867*eV, //60
		2.875*eV, 2.882*eV, 2.888*eV, 2.894*eV, 2.900*eV, 2.907*eV, 2.913*eV, 2.919*eV, 2.924*eV, 2.930*eV, //70
		2.937*eV, 2.942*eV, 2.948*eV, 2.954*eV, 2.960*eV, 2.968*eV, 2.976*eV, 2.983*eV, 2.991*eV, 3.001*eV, //80
		3.008*eV, 3.017*eV, 3.028*eV, 3.038*eV, 3.048*eV, 3.055*eV, 3.070*eV, 3.087*eV, 3.103*eV, 3.121*eV, //90
		3.138*eV, 3.155*eV, 3.173*eV, 3.191*eV, 3.220*eV, 3.250*eV, 3.281*eV, 3.313*eV, 3.344*eV, 3.375*eV, //100
		3.403*eV, 3.439*eV, 3.479*eV, 3.522*eV, 3.566*eV, 3.611*eV, 3.644*eV, 3.684*eV, 3.731*eV, 3.780*eV, //110
		3.831*eV, 3.868*eV, 3.892*eV, 3.910*eV, 3.921*eV, 3.934*eV, 3.946*eV, 3.957*eV, 3.970*eV, 3.994*eV, //120
		4.044*eV, 4.102*eV, 4.160*eV, 4.202*eV, 4.236*eV, 4.267*eV, 4.298*eV, 4.328*eV, 4.357*eV, 4.387*eV, //130
		4.422*eV, 4.455*eV, 4.494*eV, 4.563*eV, 4.607*eV, 4.616*eV, 4.624*eV, 4.627*eV, 4.628*eV, 4.633*eV, //140
		4.640*eV, 4.642*eV, 4.649*eV, 4.656*eV, 4.661*eV, 4.666*eV, 4.678*eV, 4.685*eV, 4.692*eV, 4.699*eV, //150
		4.706*eV, 4.713*eV, 4.720*eV, 4.727*eV, 4.740*eV, 4.751*eV, 4.763*eV, 4.775*eV, 4.788*eV, 4.798*eV, //160
		4.813*eV, 4.828*eV, 4.840*eV, 4.853*eV, 4.869*eV, 4.886*eV, 4.905*eV, 4.928*eV, 4.953*eV, 5.015*eV, //170
		5.099*eV, 5.143*eV, 5.174*eV, 5.202*eV, 5.235*eV, 5.265*eV, 5.294*eV, 5.330*eV, 5.413*eV, 5.493*eV, //180
		5.556*eV, 5.611*eV}; //182
	// Wavelengths in nm
	G4double wlenNM[nI];
	// Arrays for material properties
	G4double rindex_air[nI];
	G4double rindex_acrylic[nI];
	G4double aLen_acrylic[nI];
	G4double rindex_EJ200[nI];
	G4double aLen_EJ200[nI];
	G4double sc_EJ200[nI] = {
		0.000,  0.000,  0.000,  0.000,  0.000,  0.010,  0.020,  0.035,  0.050,  0.060, //10
		0.070,  0.085,  0.090,  0.095,  0.098,  0.100,  0.110,  0.120,  0.130,  0.140, //20
		0.150,  0.160,  0.170,  0.180,  0.200,  0.220,  0.240,  0.250,  0.270,  0.290, //30  
		0.300,  0.320,  0.340,  0.350,  0.360,  0.390,  0.400,  0.420,  0.430,  0.440, //40  
		0.445,  0.450,  0.460,  0.470,  0.480,  0.500,  0.550,  0.600,  0.630,  0.700, //50  
		0.730,  0.750,  0.800,  0.830,  0.850,  0.870,  0.900,  0.920,  0.940,  0.950, //60  
		0.960,  0.970,  0.980,  0.985,  0.990,  0.995,  1.000,  1.000,  1.000,  0.995, //70  
		0.990,  0.985,  0.980,  0.970,  0.960,  0.930,  0.900,  0.870,  0.850,  0.800, //80 
		0.700,  0.600,  0.500,  0.400,  0.300,  0.220,  0.130,  0.070,  0.010,  0.000, //90  
		0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //100
		0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //110  
		0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //120  
		0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //130  
		0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //140  
		0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //150  
		0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //160  
		0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //170  
		0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, //180  
		0.000,  0.000}; //182
	G4double reflectivity_SS[nI];
	G4double reflectivity_Al[nI];
	// PMT material properties variables (including Mu-metal surface properties)
	G4double refl_MuMetal[nI], eff_MuMetal[nI], specularLope_MuMetal[nI], specularSpike_MuMetal[nI], backscatter_MuMetal[nI];
	G4double rindex_BeCuPhotoCath[nI], aLen_BeCuPhotoCath[nI];
	G4double rindex_BorosilicateGlass[nI], aLen_BorosilicateGlass[nI];
	// Define Elements
	H = nist->FindOrBuildElement("H");
	Be = nist->FindOrBuildElement("Be");
	C = nist->FindOrBuildElement("C");
	O = nist->FindOrBuildElement("O");
	Si = nist->FindOrBuildElement("Si");
	Cr = nist->FindOrBuildElement("Cr");
	Fe = nist->FindOrBuildElement("Fe");
	Ni = nist->FindOrBuildElement("Ni");
	Cu = nist->FindOrBuildElement("Cu");
	Mo = nist->FindOrBuildElement("Mo");
	Pb = nist->FindOrBuildElement("Pb");
	// Define Materials
	air = nist->FindOrBuildMaterial("G4_AIR");
	EJ200 = new G4Material("EJ200", rho_EJ200, 2);
	EJ200->AddElement(C, Cfrac_EJ200*100*perCent);
	EJ200->AddElement(H, Hfrac_EJ200*100*perCent);
	acrylic = new G4Material("acrylic", rho_acrylic, 3);
	acrylic->AddElement(O, 2);
	acrylic->AddElement(C, 5);
	acrylic->AddElement(H, 8);
	aluminum = nist->FindOrBuildMaterial("G4_Al");
	concrete = nist->FindOrBuildMaterial("G4_CONCRETE");
	stainless = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
	lead = nist->FindOrBuildMaterial("G4_Pb");
	mumetal = new G4Material("mumetal", rho_mumetal, 4);
	mumetal->AddElement(Ni, 0.80);
	mumetal->AddElement(Fe, 0.15);
	mumetal->AddElement(Cu, 0.04);
	mumetal->AddElement(Mo, 0.01);
	borosilicateGlass = new G4Material("borosilicateGlass", rho_BorosilicateGlass, 2);
	borosilicateGlass->AddElement(Si, 1);
	borosilicateGlass->AddElement(O, 2);
	paper = new G4Material("paper", 0.8*g/cm3, 3);
	paper->AddElement(O, 5);
	paper->AddElement(H, 10);
	paper->AddElement(C, 6);
	plywood = new G4Material("plywood", 0.9*g/cm3, 3);
	plywood->AddElement(O, 5);
	plywood->AddElement(H, 10);
	plywood->AddElement(C, 6);
	// Set Material Properties
	for (size_t i = 0; i<nI; i++){
		// Calculate wavelength (in nm) from photon energy
		wlenNM[i] = HCNM/photonEnergy[i];// Energy is input directly above. Note: energy is displayed in units of MeV in the G4cout statement!
		//G4cout << "Photon Energy (eV:) " << photonEnergy[i] << "; wavelength (nm:) " << wlenNM[i] << G4endl;// uncomment if debug feedback is needed
		rindex_air[i] = rindexConst_air;
		rindex_acrylic[i] = rindexConst_acrylic;
		aLen_acrylic[i] = aLenConst_acrylic;
		rindex_EJ200[i] = rindexConst_EJ200;
		aLen_EJ200[i] = aLenConst_EJ200;
		reflectivity_Al[i] = refConst_Al;
		reflectivity_SS[i] = refConst_SS;
		rindex_BorosilicateGlass[i] = rindexConst_BorosilicateGlass;
		aLen_BorosilicateGlass[i] = aLenConst_BorosilicateGlass;
		//G4cout << "physical parameters: EJ-200 specular component: " << sc_EJ200[i] << G4endl;// uncomment if debug feedback is needed
	}
	//G4cout << "refractive indices, absorption lengths and reflectivities are all assumed constant for this simulation." << G4endl;
	// Declare material properties tables and populate with values. Assign tables to materials
	mpt_air = new G4MaterialPropertiesTable();
	mpt_air->AddProperty("RINDEX", photonEnergy, rindex_air, nI);
	air->SetMaterialPropertiesTable(mpt_air);
	mpt_EJ200 = new G4MaterialPropertiesTable();
	mpt_EJ200->AddProperty("RINDEX", photonEnergy, rindex_EJ200, nI);
	mpt_EJ200->AddProperty("ABSLENGTH", photonEnergy, aLen_EJ200, nI);
	mpt_EJ200->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy, sc_EJ200, nI, true);
	mpt_EJ200->AddConstProperty("SCINTILLATIONYIELD", scyld_EJ200);
	mpt_EJ200->AddConstProperty("SCINTILLATIONTIMECONSTANT1", dt_EJ200);
	mpt_EJ200->AddConstProperty("SCINTILLATIONRISETIME1", rt_EJ200);
	mpt_EJ200->AddConstProperty("RESOLUTIONSCALE", 1.0);
	mpt_EJ200->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
	EJ200->SetMaterialPropertiesTable(mpt_EJ200);
	mpt_acrylic = new G4MaterialPropertiesTable();
	mpt_acrylic->AddProperty("RINDEX", photonEnergy, rindex_acrylic, nI);
	mpt_acrylic->AddProperty("ABSLENGTH", photonEnergy, aLen_acrylic, nI);
	acrylic->SetMaterialPropertiesTable(mpt_acrylic);
	mpt_Al = new G4MaterialPropertiesTable();
	mpt_Al->AddProperty("REFLECTIVITY", photonEnergy, reflectivity_Al, nI);
	aluminum->SetMaterialPropertiesTable(mpt_Al);
	mpt_SS = new G4MaterialPropertiesTable();
	mpt_SS->AddProperty("REFLECTIVITY", photonEnergy, reflectivity_SS, nI);
	stainless->SetMaterialPropertiesTable(mpt_SS);
	mpt_BorosilicateGlass = new G4MaterialPropertiesTable();
	mpt_BorosilicateGlass->AddProperty("RINDEX", photonEnergy, rindex_BorosilicateGlass, nI);
	mpt_BorosilicateGlass->AddProperty("ABSLENGTH", photonEnergy, aLen_BorosilicateGlass, nI);
	borosilicateGlass->SetMaterialPropertiesTable(mpt_BorosilicateGlass);
	// optical surface properties
	surface_Al = new G4OpticalSurface("surface_Al");
	surface_Al->SetType(dielectric_metal);
	surface_Al->SetFinish(ground);
	surface_Al->SetModel(unified);
	surface_Al->SetMaterialPropertiesTable(mpt_Al);
	surface_SS = new G4OpticalSurface("surface_SS");
	surface_SS->SetType(dielectric_metal);
	surface_SS->SetFinish(ground);
	surface_SS->SetModel(unified);
	surface_SS->SetMaterialPropertiesTable(mpt_SS);
	// Print Materials Table
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}
// Detector Construct() function
G4VPhysicalVolume* PVTcubeDetectorConstruction::Construct()
{
	// Define Visual Attributes object for adjusting coloring and visibility of various components - set to false to make components invisible
	attr = new G4VisAttributes(false);
	// Variables for setting positions
	G4double xPos = 0., yPos = 0., zPos = 0.;
	// Declare rotation angle
	G4double phi = 90*deg;
	// Declare vectors to construct rotation matrices
	G4ThreeVector u = G4ThreeVector(0, 0, 1);
	G4ThreeVector v = G4ThreeVector(-std::sin(phi), std::cos(phi), 0);
	G4ThreeVector w = G4ThreeVector(std::cos(phi), std::sin(phi), 0);
	// Declare rotation matrices to use for source, LG and PMT, including shield wall
	G4RotationMatrix *xRot = new G4RotationMatrix(u, v, w);
	G4RotationMatrix *yRot = new G4RotationMatrix(v, u, w);
	// World volume
	solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);
	logicWorld =  new G4LogicalVolume(solidWorld, air, "logicWorld");
	// make world volume invisiblelenPMT/2
	logicWorld->SetVisAttributes(attr);
	physWorld = new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), logicWorld, "physWorld", 0, false, 0, true);
	// make optical table filler
	solidTable = new G4Box("solidTable", w_table, h_table, l_table);
	logicTable = new G4LogicalVolume(solidTable, concrete, "logicTable");
	// off-white color for concrete base
	attr = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.2));
	logicTable->SetVisAttributes(attr);
	yPos = -h_table/2;
	physTable = new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), logicTable, "physTable", logicWorld, false, 0, true);
	// make optical table top
	solidTableTop = new G4Box("solidTableTop", w_table, t_table/2, l_table);
	logicTableTop = new G4LogicalVolume(solidTableTop, stainless, "logicTableTop");
	// gray color for optical table top
	attr = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.5));
	logicTableTop->SetVisAttributes(attr);
	yPos = h_table/2 + t_table/2;
	physTableTop = new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), logicTableTop, "physTableTop", logicWorld, false, 0, true);
	// wood dark box
	solidDarkBoxOuter = new G4Box("solidDarkBoxOuter", 22*cm, 22*cm, 60*cm);
	solidDarkBoxInner = new G4Box("solidDarkBoxInner", 21*cm, 21*cm, 59*cm);
	solidDarkBox = new G4SubtractionSolid("solidDarkBox", solidDarkBoxOuter, solidDarkBoxInner, 0, G4ThreeVector(0.,0.,0.));
	logicDarkBox = new G4LogicalVolume(solidDarkBox, plywood, "logicDarkBox");
	// magenta dark box
	attr = new G4VisAttributes(G4Colour(0.9, 0.0, 0.9, 0.05));
	logicDarkBox->SetVisAttributes(attr);
	yPos = h_table/2 + t_table + 22*cm;
	physDarkBox = new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), logicDarkBox, "physDarkBox", logicWorld, false, 0, true);
	// Note pad used as platform for source
	solidNotePad = new G4Box("solidNotePad", w_notepad, h_notepad, l_notepad);
	logicNotePad = new G4LogicalVolume(solidNotePad, paper, "logicNotePad");
	// make yellow "note pad" stack
	attr = new G4VisAttributes(G4Colour(0.69,0.69,0.0,0.5));
	logicNotePad->SetVisAttributes(attr);
	yPos = h_table/2 + t_table + 1*cm + yVoxelSize/4;
	zPos = -15.0*cm;
	physNotePad = new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), logicNotePad, "physNotePad", logicWorld, false, 0, true);
	// source "puck"
	solidSource = new G4Tubs("solidSource", r_source_inner, r_source_outer, t_source, 0, 360*deg);
	logicSource = new G4LogicalVolume(solidSource, acrylic, "logicSource");
	// red source "puck"
	attr = new G4VisAttributes(G4Colour(0.9, 0.0, 0.0, 0.5));
	logicSource->SetVisAttributes(attr);
	yPos = h_table/2 + t_table + 1*cm + yVoxelSize/2+ 0.05*cm;
	physSource = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicSource, "physSource", logicWorld, false, 0, true);
	// PVT Cube
	solidVoxel = new G4Box("solidVoxel", xVoxelSize/2, yVoxelSize/2, zVoxelSize/2);
	logicVoxel = new G4LogicalVolume(solidVoxel, EJ200, "logicVoxel");
	// make yellow colored PVT cube
	attr = new G4VisAttributes(G4Colour(0.9,0.9,0.0,0.4));
	logicVoxel->SetVisAttributes(attr);
	// Set cube as scoring volume
	fPVTcubeScoringVolume = logicVoxel;
	xPos = 10*cm + xVoxelSize/2;
	yPos = h_table/2 + t_table + 1.0*cm + yVoxelSize/2;
	physVoxel = new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), logicVoxel, "physVoxel", logicWorld, false, 0, true);
	xPos = -10*cm - xVoxelSize/2;
	physVoxel = new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), logicVoxel, "physVoxel", logicWorld, false, 1, true);
	xPos = 0*cm;
	zPos = 15.0*cm;
	physVoxel = new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), logicVoxel, "physVoxel", logicWorld, false, 2, true);
	xPos = 2.5*xVoxelSize;
	physVoxel = new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), logicVoxel, "physVoxel", logicWorld, false, 3, true);
	// Foil "wrap" for PVT cubes
	solidFoilWrapInner = new G4Box("solidFoilWrapInner", xVoxelSize/2, yVoxelSize/2, zVoxelSize/2);
	solidFoilWrapOuter = new G4Box("solidFoilWrapOuter", xVoxelSize/2 + t_foil, yVoxelSize/2, zVoxelSize/2 + t_foil);
	solidFoilWrap = new G4SubtractionSolid("solidFoilWrap", solidFoilWrapOuter, solidFoilWrapInner, 0, G4ThreeVector(0.,0.,0.));
	logicFoilWrap = new G4LogicalVolume(solidFoilWrap, aluminum, "logicFoilWrap");
	// gray color for foil "wrap"
	attr = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.5));
	logicFoilWrap->SetVisAttributes(attr);
	xPos = 10*cm + xVoxelSize/2;
	yPos = h_table/2 + t_table + 1.0*cm + yVoxelSize/2;
	zPos = -15.0*cm;
	physFoilWrap = new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), logicFoilWrap, "physFoilWrap", logicWorld, false, 0, true);
	xPos = -10*cm - xVoxelSize/2;
	physFoilWrap = new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), logicFoilWrap, "physFoilWrap", logicWorld, false, 1, true);
	xPos = 0*cm;
	zPos = 15.0*cm;
	physFoilWrap = new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), logicFoilWrap, "physFoilWrap", logicWorld, false, 2, true);
	xPos = 2.5*xVoxelSize;
	physFoilWrap = new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), logicFoilWrap, "physFoilWrap", logicWorld, false, 3, true);
	// LG for PMT
	solidLGTrd = new G4Trd("solidLGTrd", xVoxelSize/2, lenLGBase/2, yVoxelSize/2, lenLGBase/2, lenLGTaper/2);
	solidLGCone = new G4Cons("solidLGCone", 0.*cm, r1_LG, 0.*cm, r2_LG, lenLGTaper/2, 0, 360*deg);
	// Define the light guide from the intersection of the two solids defined above
	solidLG = new G4IntersectionSolid("solidLG", solidLGTrd, solidLGCone);
	logicLG = new G4LogicalVolume(solidLG, acrylic, "logicLG");
	// make the LG orange
	attr = new G4VisAttributes(G4Colour(0.7,0.3,0.0,0.4));
	logicLG->SetVisAttributes(attr);
	xPos = 10*cm + xVoxelSize/2;
	yPos = h_table/2 + t_table + 1.0*cm + yVoxelSize + lenLGTaper/2;
	zPos = -15.0*cm;
	physLG = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicLG, "physLG", logicWorld, false, 0, true);
	xPos = -10*cm - xVoxelSize/2;
	physLG = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicLG, "physLG", logicWorld, false, 1, true);
	xPos = 0*cm;
	zPos = 15.0*cm;
	physLG = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicLG, "physLG", logicWorld, false, 2, true);
	xPos = 2.5*xVoxelSize;
	physLG = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicLG, "physLG", logicWorld, false, 3, true);
	// Define Borosilicate Glass boundaries
	solidPMTGlass = new G4Tubs("solidPMTGlass", 0, rPMT, tGlass/2, 0, 360.*deg);
	solidPMTconvex = new G4Sphere("solidPMTconvex" , 0, rSpherePMTsurf, 0, 360.*deg, 0, 360.*deg);
	// Define PMT Lens (to simulate reflections at PMT/LG boundary)
	solidPMTLens = new G4SubtractionSolid("solidPMTLens", solidPMTGlass, solidPMTconvex, 0, G4ThreeVector(0,0,rSpherePMTsurf-tGlass/2+tGlass_min));
	logicPMTLens = new G4LogicalVolume(solidPMTLens, borosilicateGlass, "logicPMTLens");
	// make lenses magenta
	attr = new G4VisAttributes(G4Colour(0.3,0.0,0.2,0.4));
	logicPMTLens->SetVisAttributes(attr);
	xPos = 10*cm + xVoxelSize/2;
	yPos = h_table/2 + t_table + 1.0*cm + yVoxelSize + lenLGTaper + tGlass/2;
	zPos = -15.0*cm;
	physPMTLens = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicPMTLens, "physPMTLens", logicWorld, false, 0, true);
	xPos = -10*cm - xVoxelSize/2;
	physPMTLens = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicPMTLens, "physPMTLens", logicWorld, false, 1, true);
	xPos = 0*cm;
	zPos = 15.0*cm;
	physPMTLens = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicPMTLens, "physPMTLens", logicWorld, false, 2, true);
	xPos = 2.5*xVoxelSize;
	physPMTLens = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicPMTLens, "physPMTLens", logicWorld, false, 3, true);
	// Define "PMT(s)"
	solidPMT = new G4Tubs("solidPMT", 0, rPMT, lenPMT/2, 0, 360*deg);
	logicPMT = new G4LogicalVolume(solidPMT, air, "logicPMT");
	// make PMT volumes purple cylinders
	attr = new G4VisAttributes(G4Colour(0.5,0.0,0.5,0.5));
	logicPMT->SetVisAttributes(attr);
	xPos = 10*cm + xVoxelSize/2;
	yPos = h_table/2 + t_table + 1.0*cm + yVoxelSize + lenLGTaper + tGlass + lenPMT/2;
	zPos = -15.0*cm;
	physPMT = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicPMT, "physPMT", logicWorld, false, 0, true);
	xPos = -10*cm - xVoxelSize/2;
	physPMT = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicPMT, "physPMT", logicWorld, false, 1, true);
	xPos = 0*cm;
	zPos = 15.0*cm;
	physPMT = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicPMT, "physPMT", logicWorld, false, 2, true);
	xPos = 2.5*xVoxelSize;
	physPMT = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicPMT, "physPMT", logicWorld, false, 3, true);
	// Define "PMT" housing
	solidPMTshield = new G4Tubs("solidPMTshield", 0, rPMT, lenPMT/2 + (rPMT - r2_LG)/2, 0, 360*deg);
	solidPMTshieldwall = new G4SubtractionSolid("solidPMTshieldwall", solidPMTshield, solidPMT);
	logicPMTshieldwall = new G4LogicalVolume(solidPMTshieldwall, mumetal, "logicPMTshieldwall");
	// make mumetal shield gray and almost transparent
	attr = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.05));
	logicPMTshieldwall->SetVisAttributes(attr);
	xPos = 10*cm + xVoxelSize/2;
	zPos = -15.0*cm;
	physPMTshieldwall = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicPMTshieldwall, "physPMTshieldwall", logicWorld, false, 0, true);
	xPos = -10*cm - xVoxelSize/2;
	physPMTshieldwall = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicPMTshieldwall, "physPMTshieldwall", logicWorld, false, 1, true);
	xPos = 0*cm;
	zPos = 15.0*cm;
	physPMTshieldwall = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicPMTshieldwall, "physPMTshieldwall", logicWorld, false, 2, true);
	xPos = 2.5*xVoxelSize;
	physPMTshieldwall = new G4PVPlacement(yRot, G4ThreeVector(xPos, yPos, zPos), logicPMTshieldwall, "physPMTshieldwall", logicWorld, false, 3, true);
	// lead shielding blocks
	solidLeadBlock = new G4Box("solidLeadBlock", x_block/2, y_block/2, z_block/2);
	logicLeadBlock = new G4LogicalVolume(solidLeadBlock, lead, "logicLeadBlock");
	// make lead shielding block dark gray and almost opaque
	attr = new G4VisAttributes(G4Colour(0.25,0.25,0.25,0.95));
	logicLeadBlock->SetVisAttributes(attr);
	xPos = -3.0*cm - x_block;
	yPos = h_table/2 + t_table + 1.0*cm + y_block/2;
	zPos = 0;
	physLeadBlock = new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), logicLeadBlock, "physLeadBlock", logicWorld, false, 0, true);
	xPos = 3.0*cm + x_block;
	physLeadBlock = new G4PVPlacement(0, G4ThreeVector(xPos, yPos, zPos), logicLeadBlock, "physLeadBlock", logicWorld, false, 1, true);
	// return value
	return physWorld;
}
// Sensitive Detector Function
void PVTcubeDetectorConstruction::ConstructSDandField()
{
	// Get SDM pointer for creating new detectors
	G4SDManager *SDman = G4SDManager::GetSDMpointer();
	// Define PMT's as sensitive volumes
	detPMT = new PVTcubePMTSensitiveDetector("/PVTcubePMT");
	SDman->AddNewDetector(detPMT);
	logicPMT->SetSensitiveDetector(detPMT);
	// PVT cube voxel scoring volume
	detVoxel = new PVTcubeVoxelSensitiveDetector("/PVTcubeVoxel", nCubes);
	SDman->AddNewDetector(detVoxel);
	logicVoxel->SetSensitiveDetector(detVoxel);
}
