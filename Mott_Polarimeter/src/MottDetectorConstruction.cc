//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: MottDetectorConstruction.cc,v 3.16 2014/01/16 23:40:03 mjmchugh Exp mjmchugh $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MottDetectorConstruction.hh"
#include "MottDetectorMessenger.hh"
#include "MottDumpParameterisation.hh" 
#include "MottTrackerSD.hh"
//#include "MottPMTSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh" 
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
MottDetectorConstruction::MottDetectorConstruction() {
 
  // Initialize member variables
  solidWorld = NULL;
  logicWorld = NULL;
  physiWorld = NULL;
  
  WorldMater = NULL;
  fWorldLength = 8.0*m;
  
  solidTarget = NULL;
  logicTarget = NULL;
  physiTarget = NULL;
  
  TargetMater = NULL;
  fTargetLength = 1.0e-6*m;	// One micron
  TargetXpos = 0.0*mm;		//
  TargetYpos = 0.0*mm;		// Target home position
  TargetZpos = -3.9878*mm;	//
  
  stepLimit = NULL;

  detectorMessenger = new MottDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
MottDetectorConstruction::~MottDetectorConstruction()
{
  delete stepLimit;
  delete detectorMessenger;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* MottDetectorConstruction::Construct()
{
  std::cout << "\tEntering MottDetectorConstruction::Construct()" << std::endl;

//--------- Material definition ---------

  G4String name, symbol;
  G4double a, z;
  G4double density, temperature, pressure;
  G4double fractionmass;
  G4int nel;

  // ---- Elemental Materials ---- //

  G4Material* Be = new G4Material("Beryllium", z=4, a=9.0121831*g/mole, density=1.85*g/cm3);
  G4Material* Al = new G4Material("Aluminum", z=13, a=26.981*g/mole, density=2.699*g/cm3);
  G4Material* Cu = new G4Material("Copper", z=29, a=63.564*g/mole, density=8.96*g/cm3);
  G4Material* Au = new G4Material("Gold", z=79, a=169.97*g/mole, density=19.32*g/cm3);
  G4Material* Pb = new G4Material("Lead", z=82, a=207.2*g/mole, density=11.34*g/cm3);
  G4Material* Cs = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cs");

  // ---- Composite Materials ---- //
  
  // Stainless Steel Components
  G4Element* Fe = new G4Element(name="Iron", symbol="Fe", z=26, a=55.895*g/mole);
  G4Element* Cr = new G4Element(name="Chromium", symbol="Cr", z=24, a=51.996*g/mole);
  G4Element* Ni = new G4Element(name="Nickel", symbol="Ni", z=28, a=58.693*g/mole);
  G4Element* Mn = new G4Element(name="Manganese", symbol="Mn", z=25, a=54.938*g/mole);
  
  G4Material* Steel = new G4Material(name="StainlessSteel", density=7.74*g/cm3, 
                                     nel=4, kStateSolid, STP_Temperature);
  Steel->AddElement(Fe, fractionmass=0.7075);
  Steel->AddElement(Cr, fractionmass=0.19);
  Steel->AddElement(Ni, fractionmass=0.0925);
  Steel->AddElement(Mn, fractionmass=0.01);
  
  G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4Material* Vacuum = new G4Material("Vacuum", density=0.5e-13*g/cm3, 
  				      nel=1, kStateGas, temperature=STP_Temperature, 
  				      pressure=1.0e-10*bar);
  Vacuum->AddMaterial(Air, fractionmass=1.0);
  
  G4Material* Plastic = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");  

  // Print all the materials defined.
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

///////////////////////////////////////////////////////////////////////////
//////////////// Principal Geometry (sizes of solids) /////////////////////
///////////////////////////////////////////////////////////////////////////

  // Inch definition.
  static const G4double inch = 25.400*mm;
  
  // World Volume
  WorldMater = Vacuum; /* Vacuum */			// World Material = Vacuum

  // Target Volume
  fTargetLength = 1.0e-05*cm;                       	// Full length of Target
  TargetMater = Au;					// Target material = gold
  
  // Stainless Steel Scattering Chamber pieces
  G4double fMainChamberLength = 20.56*inch;		// Full length of main chamber 
  G4double MainChamberOR = 4.00*inch;			// Outer radius of chamber
  G4double MainChamberIR = 3.88*inch;			// Inner radius of chamber
  G4double fChamberFaceBigLength = 0.438*inch;		// Full length of end flange 
  G4double fChamberFaceLength = 0.94*inch;		// Full length of end flange
  G4double ChamberFaceIR = 0.685*inch; 			// Beam pipe outer radius
  G4double ChamberFaceBigRadius = 4.92*inch;		// Big outer radius
  G4double ChamberFaceRadius = 3.88*inch;		// Radius for part inside chamber
  G4double PortTubeOR = 1.5*inch;			// Port tube outer radius
  G4double PortTubeIR = 1.435*inch;			// Port tube inner radius
  G4double fPortTubeLength = 4.11*inch;			// Full length of port tube  
  
  // Stainless Steel Pump Chamber
  G4double PumpChamberIR = MainChamberIR;		// Inner radius of chamber
  G4double PumpChamberOR = MainChamberOR;		// Outer radius of chamber
  G4double fPumpChamberLength = 9.875*inch;		// Full length of chamber
    
  // Aluminum Extension Tube
  G4double ExtensionTubeIR = MainChamberIR;		// Inner radius of tube
  G4double ExtensionTubeOR = MainChamberOR;		// Outer radius of tube
  G4double fExtensionTubeLength = 58.125*inch;		// Full length of tube
    
  // Aluminum Scattering Chamber Liner. 
  G4double ChamberLinerOR = 3.75*inch;			// Outer radius of chamber liner 
  G4double ChamberLinerIR = 3.25*inch; 			// Inner radius of chamber liner
  G4double fChamberLinerLength = 11.56*inch;		// Full length of chamber liner
  
  // Aluminum Baffle/Collimator.
  G4double fAlBaffleLength = 1.0*inch - 0.04*inch;      // Thickness  
  G4double fAlBaffleWidth = 4.19*inch;		 	// "Square" face Dimension. 
  G4double AlBaffleBeamHoleRadius = 0.5*1.0*inch;      	// Center hole radius.
  G4double AlBaffleCollDSHoleLocation = 0.5*2.85*inch;	// Distance from center of baffle to center of DS collimator hole.
  G4double AlBaffleCollUSHoleLocation = 0.5*3.11*inch;	// Distance from center of baffle to center of US collimator hole.
    
  // Aluminum Nosepiece/collimator
  G4double AlNose1OutRadius = 0.75*inch;		// Front piece largest radius.
  G4double fAlNose1Length = 0.25*inch;			// Full front piece depth.
  G4double AlNose2OutRadius = 0.2495*inch;		// Same as front inner radius.
  G4double AlNoseInRadius = 0.19*inch;			// True inner radius of collimator. 
  G4double fAlNose2Length = 0.75*inch;			// Full length of collimator.
  G4Material* AlNoseMater = Al;                         // Nose is made of aluminium.
  
  // Lead Collimator
  G4double PbInRadius = 0.25*inch;			// Inner Radius of Lead Collimator
  G4double PbOutRadius = 1.425*inch;			// Normal Radius
  G4double PbBigRadius = 2.00*inch;			// Big part radius
  G4double PbSmallRadius = 1.365*inch; 			// Small Part Radius.	
  G4double fPbBigSize = 0.50*inch;			// Length of Big section.
  G4double PbBigPos = 0.25*inch;			// Big section is 1/4" back on lead.
  G4double fPbNormSize = 2.00*inch;			// Length of whole thing.
  G4double PbSmallSize = 0.31*inch;			// Part of the skinny that sticks out.
  
  // Air Column
  G4double fAirLength = 2.0*inch;			// Length of air column in front of the  
  G4double AirRadius = 0.19*inch; 
  
  // Detector Housing
  G4double HouseFrontFaceOR = 2.0625*inch;		// Front edge of detector housing
  G4double fHouseFrontFaceLength = 0.25*inch;		// Length of the front face with boltholes.
  G4double HouseFrontOR = 1.50*inch;			// Outer radius of section around lead.
  G4double HouseFrontIR = 1.435*inch;                   // Only +/- .01" clearance around lead.
  G4double fHouseFrontLength = 1.00*inch; 		// Length of skinny section.  
  G4double fHouseWeldLength = 0.125*inch; 		// Section joining to can to front
  G4double HouseCanOR = 1.75*inch;			// Long section outer radius.
  G4double HouseCanIR = 1.665*inch;			// Wall is 0.035" thick.
  G4double fHouseCanLength = 11.5625*inch; 		// Brings total housing length to 12+15/16"
  
  // Delta E detector specs
  G4double fDeltaElength = 1.0*mm;			// Full length of dE wafer
  G4double fDeltaEwidth = 1.0*inch;			// Full edge-length of transverse face 
  G4Material* DeltaEMater = Plastic;			// Delta E material = plastic

  // Energy Detector specs
  G4double fEDetLength = 2.50*inch;			// Full length of E detector
  G4double EDetRadius = 1.5*inch;			// Radius of E detector
  G4Material* EDetMater = Plastic;			// Energy Detector material = plastic  

  // "Photocathode"s
  G4Material* PhotoCathodeMaterial = Cs;
  
  // dE
  G4double fdEPhotoLength = 1.0*mm;
  G4double fdEPhotoWidth = 1.0*inch;
  
  // E
  G4double fEPhotoLength = 1.0*mm;
  G4double EPhotoRadius = 1.5*inch;
  
  // Scattering angle
  G4double ScatteringAngle = 172.6*deg;
  
///////////////////////////////////////////////////////////////////////////
////////////////////////////// World Def //////////////////////////////////
///////////////////////////////////////////////////////////////////////////
  
  //------------------------------ 
  // World
  //------------------------------ 

  G4double HalfWorldLength = 0.5*fWorldLength;
 
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(fWorldLength);
  //G4cout << "Computed tolerance = "
  //       << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
  //       << " mm" << G4endl;

  solidWorld= new G4Box("world",HalfWorldLength,HalfWorldLength,HalfWorldLength);
  logicWorld= new G4LogicalVolume( solidWorld, WorldMater, "World", 0, 0, 0);
  
  // Must place the World Physical volume unrotated at (0,0,0).
  physiWorld = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
                                 logicWorld,      // its logical volume
                                 "World",         // its name
                                 0,               // its mother  volume
                                 false,           // no boolean operations
                                 0);              // copy number  

///////////////////////////////////////////////////////////////////////////
////////////// Vacuum Components (Scattering chamber, etc.) ///////////////
///////////////////////////////////////////////////////////////////////////

  G4RotationMatrix* RotateForUp = new G4RotationMatrix();
                    RotateForUp->rotateX(ScatteringAngle);
  G4RotationMatrix* RotateForDown = new G4RotationMatrix();
                    RotateForDown->rotateX(-ScatteringAngle);
  G4RotationMatrix* RotateForLeft = new G4RotationMatrix();
                    RotateForLeft->rotateY(-ScatteringAngle);
  G4RotationMatrix* RotateForRight = new G4RotationMatrix();
                    RotateForRight->rotateY(ScatteringAngle);

  G4ThreeVector unionShift;

  //------------------------------ 
  // Target
  //------------------------------
  
  G4double targetSize = 0.5*fTargetLength;    	// Half length of the Target
  
  G4ThreeVector positionTarget = G4ThreeVector(TargetXpos,TargetYpos,TargetZpos);
   
  solidTarget = new G4Tubs("solidTarget",
                           0.0*cm,
                           2.54*cm,
                           targetSize,
                           0.0*deg,
                           360.0*deg);
  logicTarget = new G4LogicalVolume(solidTarget, TargetMater, "logicTarget");
  physiTarget = new G4PVPlacement(0,               // no rotation
                                  positionTarget,  // at (x,y,z)
                                  logicTarget,     // its logical volume
                                  "physiTarget",   // its name
                                  logicWorld,      // its mother volume
                                  false,           // no boolean operations
                                  0);              // copy number
  
  //G4cout << "Target is " << fTargetLength/cm << " cm of " 
  //       << TargetMater->GetName() << G4endl;

  //------------------------------ 
  // Aluminium Baffle Plate
  //------------------------------
  
  G4double AlBaffleLength = 0.5*fAlBaffleLength;
  G4double AlBaffleWidth = 0.5*fAlBaffleWidth;
  G4double BaffleZpos = -10.929*inch - AlBaffleLength;

  G4double CollimatorLocation = 0.5*( AlBaffleCollDSHoleLocation + AlBaffleCollUSHoleLocation );
  
  G4ThreeVector PositionAlBaffle = G4ThreeVector(0,0,BaffleZpos);
  
  G4Box* Baffle1 = new G4Box("Baffle1", AlBaffleWidth, AlBaffleWidth, AlBaffleLength);
  G4Tubs* Baffle2 = new G4Tubs("Baffle2", 0, AlBaffleBeamHoleRadius, AlBaffleLength+1.0*mm, 0*deg, 360*deg);
  
  unionShift = G4ThreeVector(0,0,0);
  G4SubtractionSolid* Baffle3 = new G4SubtractionSolid("Baffle3", Baffle1, Baffle2, 0, unionShift);
  
  G4Cons* CollimatorCutout = new G4Cons("Collimator Cutout",	// Name
                                        0,			// Min R of low-z end 
                                        0.064*inch,		// Max R of low-z end
                                        0,			// Min R of high-z end 
                                        0.190*inch,		// Min R of high-z end 
                                        1.00*inch,		// Half-Length in z
                                        0,			// Starting phi-angle
                                        360*deg);		// Ending phi-angle 
  
   
  unionShift = G4ThreeVector(CollimatorLocation,0,0);
  G4SubtractionSolid* Baffle4 = new G4SubtractionSolid("Baffle4", Baffle3, CollimatorCutout, RotateForLeft, unionShift);
  unionShift = G4ThreeVector(-CollimatorLocation,0,0);
  G4SubtractionSolid* Baffle5 = new G4SubtractionSolid("Baffle5", Baffle4, CollimatorCutout, RotateForRight, unionShift);
  unionShift = G4ThreeVector(0,CollimatorLocation,0);
  G4SubtractionSolid* Baffle6 = new G4SubtractionSolid("Baffle6", Baffle5, CollimatorCutout, RotateForUp, unionShift);
  unionShift = G4ThreeVector(0,-CollimatorLocation,0);
  G4SubtractionSolid* Baffle7 = new G4SubtractionSolid("Baffle7", Baffle6, CollimatorCutout, RotateForDown, unionShift);
  
  G4LogicalVolume* logicAlBaffle = new G4LogicalVolume(Baffle7, Al, "logicAlBaffle");
  new G4PVPlacement(0, PositionAlBaffle, logicAlBaffle, "PhysiAlBaffle", logicWorld, false, 0);

  //------------------------------ 
  // Stainless Steel Scattering Chamber
  //------------------------------

  // Half Lengths
  G4double MainChamberLength = 0.5*fMainChamberLength;
  
  G4double MainChamberOffset = -7.332*inch;  // calculated from drawing 39300-C-0052

  G4ThreeVector positionMainChamber = G4ThreeVector(0., 0., MainChamberOffset);
  
  // Main Chamber
  G4Tubs* solidMainChamber = new G4Tubs("solidMainChamber",
  					 MainChamberIR,
  					 MainChamberOR,
  					 MainChamberLength,
  					 0.0*deg,
  					 360.0*deg);
  G4LogicalVolume* logicMainChamber = new G4LogicalVolume(solidMainChamber, Al,
  							   "logicMainChamber");
  new G4PVPlacement(0, positionMainChamber, logicMainChamber, 
                    "physiMainChamber", logicWorld, false, 0);

  //------------------------------ 
  // Scattering Chamber End Flange & Port Tubes
  //------------------------------

  G4double ChamberFaceBigLength = 0.5*fChamberFaceBigLength;		// Half length of end flange. 
  G4double ChamberFaceLength = 0.5*fChamberFaceLength;    
  G4double PortTubeLength = 0.5*fPortTubeLength; 
 
  // Calculated from 
  G4double EndPlatePlacement = 17.612*inch + fChamberFaceBigLength - ChamberFaceLength;
  			
  G4ThreeVector positionEndFlange = G4ThreeVector(0. , 0., -EndPlatePlacement);
  
  // Chamber end flange with scatering ports. 					 
  G4Tubs* solidChamberEnd1 = new G4Tubs("solidChamberEnd1",
                                        ChamberFaceIR,
                                        ChamberFaceRadius,
                                        ChamberFaceLength,
                                        0.0*deg,
                                        360.0*deg);
  G4Tubs* solidChamberEnd2 = new G4Tubs("solidChamberEnd2",
                                        ChamberFaceRadius,
                                        ChamberFaceBigRadius,
                                        ChamberFaceBigLength,
                                        0.0*deg,
                                        360.0*deg);                                        
  G4Tubs* solidPortTubeCutout = new G4Tubs("solidPortTubeCutout",
  					  0.0*mm,
  					  PortTubeOR,
  					  PortTubeLength,
  					  0.0*deg,
  					  360.0*deg);
  G4Tubs* solidPortTube = new G4Tubs("solidPortTube",
  				    PortTubeIR,
  				    PortTubeOR,
  				    PortTubeLength,
  				    0.0*deg,
  				    360.0*deg);
  					  
  // Movement so they line up properly.  
  unionShift = G4ThreeVector(0, 0, -ChamberFaceLength+ChamberFaceBigLength);  
  G4UnionSolid* solidEndFlange1 = new G4UnionSolid("solidEnd1",
                                              	   solidChamberEnd1,
                                                   solidChamberEnd2,
                                                   0,
                                              	   unionShift);
  // Cut UP hole
  unionShift = G4ThreeVector(0, 2.354*inch, -ChamberFaceLength);
  G4SubtractionSolid* solidEndFlange2 = new G4SubtractionSolid("solidEnd2",
                                              	   		solidEndFlange1,
                                                   		solidPortTubeCutout,
                                                   		RotateForUp,
                                              	   		unionShift);
  // Cut DOWN hole
  unionShift = G4ThreeVector(0, -2.354*inch, -ChamberFaceLength);
  G4SubtractionSolid* solidEndFlange3 = new G4SubtractionSolid("solidEnd3",
                                              	   		solidEndFlange2,
                                                   		solidPortTubeCutout,
                                                   		RotateForDown,
                                              	   		unionShift);
  // Cut LEFT hole
  unionShift = G4ThreeVector(2.354*inch, 0, -ChamberFaceLength);
  G4SubtractionSolid* solidEndFlange4 = new G4SubtractionSolid("solidEnd4",
                                              	   		solidEndFlange3,
                                                   		solidPortTubeCutout,
                                                   		RotateForLeft,
                                              	   		unionShift);
  // Cut RIGHT hole
  unionShift = G4ThreeVector(-2.354*inch, 0, -ChamberFaceLength);
  G4SubtractionSolid* solidEndFlange5 = new G4SubtractionSolid("solidEnd5",
                                              	   		solidEndFlange4,
                                                   		solidPortTubeCutout,
                                                   		RotateForRight,
                                              	   		unionShift);
  // Add LEFT Port Tube
  unionShift = G4ThreeVector(2.55796*inch, 0.0, -2.06214*inch);  // Derived from drawing 39300-C-0051
  G4UnionSolid* solidEndFlange6 = new G4UnionSolid("solidEnd6",
  						  solidEndFlange5,
  						  solidPortTube,
  						  RotateForLeft,
  						  unionShift);
  // Add UP Port Tube
  unionShift = unionShift.rotateZ(90.0*deg);
  G4UnionSolid* solidEndFlange7 = new G4UnionSolid("solidEnd7",
  						  solidEndFlange6,
  						  solidPortTube,
  						  RotateForUp,
  						  unionShift);
  // Add RIGHT Port Tube
  unionShift = unionShift.rotateZ(90.0*deg);
  G4UnionSolid* solidEndFlange8 = new G4UnionSolid("solidEnd8",
  						  solidEndFlange7,
  						  solidPortTube,
  						  RotateForRight,
  						  unionShift);
  // Add DOWN Port Tube
  unionShift = unionShift.rotateZ(90.0*deg);
  G4UnionSolid* solidEndFlange9 = new G4UnionSolid("solidEnd9",
  						  solidEndFlange8,
  						  solidPortTube,
  						  RotateForDown,
  						  unionShift);
  
  G4LogicalVolume* logicEndFlange = new G4LogicalVolume(solidEndFlange9, Steel, "logicEndFlange");
  new G4PVPlacement(0, positionEndFlange, logicEndFlange, "physiEndFlange", logicWorld, false, 0);
                                             
  //------------------------------ 
  // Aluminum Scattering Chamber Liner
  //------------------------------

  G4double ChamberLinerLength = 0.5*fChamberLinerLength;
  G4double ChamberLinerOffset = -1.84*inch; 		// Position offset for the liner.  
  
  G4ThreeVector positionChamberLiner = G4ThreeVector(0., 0., ChamberLinerOffset);
  
  G4Tubs* solidChamberLiner = new G4Tubs("solidChamberLiner",
  					 ChamberLinerIR,
  					 ChamberLinerOR,
  					 ChamberLinerLength,
  					 0.0*deg,
  					 360*deg);
  G4LogicalVolume* logicChamberLiner = new G4LogicalVolume(solidChamberLiner, Al,
  							   "logicChamberLiner");
  new G4PVPlacement(0, positionChamberLiner, logicChamberLiner, 
                    "physiChamberLiner", logicWorld, false, 0);    

  //------------------------------
  // Steel Pump Chamber
  //------------------------------
  
  G4double PumpChamberLength = 0.5*fPumpChamberLength;
  G4double PumpChamberZpos = MainChamberOffset + (MainChamberLength+PumpChamberLength);
  G4ThreeVector positionPumpChamber = G4ThreeVector(0., 0., PumpChamberZpos);
  
  G4Tubs* solidPumpChamber = new G4Tubs("solidTarget",
  					PumpChamberIR,
  					PumpChamberOR,
  					PumpChamberLength,
  					0.0*deg,
  					360.0*deg);
  G4LogicalVolume* logicPumpChamber = new G4LogicalVolume(solidPumpChamber, Steel, "logicPumpChamber");
  new G4PVPlacement(0, positionPumpChamber, logicPumpChamber, "physiPumpChamber", logicWorld, false, 0);
  
  //------------------------------
  // Aluminum Extension Tube
  //------------------------------
  
  G4double ExtensionTubeLength = 0.5*fExtensionTubeLength;
  G4double ExtensionTubeZpos = PumpChamberZpos + (PumpChamberLength + ExtensionTubeLength);
  G4ThreeVector positionExtensionTube = G4ThreeVector(0., 0., ExtensionTubeZpos);
  
  G4Tubs* solidExtensionTube = new G4Tubs("solidExtensionTube",
  					  ExtensionTubeIR,
  					  ExtensionTubeOR,
  					  ExtensionTubeLength,
  					  0.0*deg,
  					  360.0*deg);
  G4LogicalVolume* logicExtensionTube = new G4LogicalVolume(solidExtensionTube, Al, "logicExtensionTube");
  new G4PVPlacement(0, positionExtensionTube, logicExtensionTube, "physiExtensionTube", logicWorld, false, 0);

  //------------------------------
  // Be/Cu Dump Plate
  //------------------------------

  const G4int nLayersBe = 6;
  const G4int nLayersCu = 18;
  const G4int nRings = 20;
  G4double LayerDepth = (127.0/120.0)*mm; 		// Layer depth
  G4double HalfLayerDepth = 0.5*LayerDepth;		// Half-Layer depth
  G4double radialStep = 5.08*mm;

  G4double BeDumpZpos = ExtensionTubeZpos + (ExtensionTubeLength + HalfLayerDepth*nLayersBe ) + 1.53*inch;	// Center Z position of first layer 
  G4double CuDumpZpos = BeDumpZpos + HalfLayerDepth*(nLayersBe+nLayersCu);
  
  G4ThreeVector BeDumpPosition(0,0,BeDumpZpos);
  G4ThreeVector CuDumpPosition(0,0,CuDumpZpos);
  
  // Be Mother Volume
  G4Tubs* solidBeDump = new G4Tubs("solidBeDump",0, nRings*radialStep, nLayersBe*HalfLayerDepth, 0.0*deg, 360.0*deg);
  G4LogicalVolume* logicBeDump = new G4LogicalVolume(solidBeDump, Vacuum, "logicBeDump");
  new G4PVPlacement(0, BeDumpPosition, logicBeDump, "physiBeDump", logicWorld, false, 0);
  
  // Be Daughter Volumes (Parameterised)
  G4Tubs* solidBeDumpSegment = new G4Tubs("solidBeDumpSegment", 0, radialStep, HalfLayerDepth, 0.0*deg, 360.0*deg);
  G4LogicalVolume* logicBeDumpSegment = new G4LogicalVolume(solidBeDumpSegment, Be, "logicBeDumpSegment");
  G4VPVParameterisation* BeDumpParam = new MottDumpParameterisation(nLayersBe,
  								    -(nLayersBe-1)*HalfLayerDepth,
  								    LayerDepth,
  								    LayerDepth,
  								    nRings,
  								    radialStep);
  new G4PVParameterised("physiBeDumpSegments", logicBeDumpSegment, logicBeDump, kZAxis, nLayersBe*nRings, BeDumpParam);
  
  // Cu Mother Volume
  G4Tubs* solidCuDump = new G4Tubs("solidCuDump",0, nRings*radialStep, nLayersCu*HalfLayerDepth, 0.0*deg, 360.0*deg);
  G4LogicalVolume* logicCuDump = new G4LogicalVolume(solidCuDump, Vacuum, "logicCuDump");
  new G4PVPlacement(0, CuDumpPosition, logicCuDump, "physiCuDump", logicWorld, false, 0);
  
  // Cu Daugter Volumes (Parameterised)
  G4Tubs* solidCuDumpSegment = new G4Tubs("solidCuDumpSegment", 0, radialStep, HalfLayerDepth, 0.0*deg, 360.0*deg);
  G4LogicalVolume* logicCuDumpSegment = new G4LogicalVolume(solidCuDumpSegment, Cu, "logicCuDumpSegment");
  G4VPVParameterisation* CuDumpParam = new MottDumpParameterisation(nLayersCu,
  								    -(nLayersCu-1)*HalfLayerDepth,
  								    LayerDepth,
  								    LayerDepth,
  								    nRings,
  								    radialStep);
  new G4PVParameterised("physiCuDumpPlate", logicCuDumpSegment, logicCuDump, kZAxis, nLayersCu*nRings, CuDumpParam);

///////////////////////////////////////////////////////////////////////////
//////////////////////// DETECTOR PACKAGE(S) //////////////////////////////
///////////////////////////////////////////////////////////////////////////
  
  // Front of the Aluminium nosepiece is this distance from the Origin.
  // This sets the distance to the Detector Package as a whole.
  G4double distToDetPackage = 22.18*inch;

  //------------------------------
  // Aluminium Detector Collimator/Nosepeice
  //------------------------------
  
  G4double AlNose1Length = 0.5*fAlNose1Length;
  G4double AlNose2Length = 0.5*fAlNose2Length;

  G4double AlNosePlacement = distToDetPackage + AlNose1Length; 
  
  // Make placements.
  G4ThreeVector positionDetectorUp = G4ThreeVector(0, 0, AlNosePlacement);
  		positionDetectorUp.rotateX(-ScatteringAngle);
  G4ThreeVector positionDetectorDown = G4ThreeVector(0, 0, AlNosePlacement);
  		positionDetectorDown.rotateX(ScatteringAngle);
  G4ThreeVector positionDetectorLeft = G4ThreeVector(0, 0, AlNosePlacement);
  		positionDetectorLeft.rotateY(ScatteringAngle);
  G4ThreeVector positionDetectorRight = G4ThreeVector(0, 0, AlNosePlacement);
  		positionDetectorRight.rotateY(-ScatteringAngle);

  // The first piece. Al plate out front.
  G4Tubs* solidAlNose1 = new G4Tubs("solidAlNose1", 
               			    AlNoseInRadius,
               			    AlNose1OutRadius,
               			    AlNose1Length,
               			    0.0*deg,
               			    360.0*deg);

  // The second piece. Al shim inside.
  G4Tubs* solidAlNose2 = new G4Tubs("solidAlNose2", 
               			    AlNoseInRadius,
               			    AlNose2OutRadius,
               			    AlNose2Length,
               			    0.0*deg,
               			    360.0*deg);
  
  unionShift = G4ThreeVector(0, 0, AlNose2Length-AlNose1Length);
  
  //Full Aluminium Nosepiece. 
  G4UnionSolid* solidAlNose = new G4UnionSolid("solidAlNose",	// Name of the unified solid 
                                               solidAlNose1,    // First solid sets coordinates.
                                               solidAlNose2,	// Second solid
                                               0,		// Rotation for second first's frame
                                               unionShift);	// Translation for second solid
  G4LogicalVolume* logicAlNose = new G4LogicalVolume(solidAlNose, AlNoseMater, "logicAlNose");
 
  //Place all 4 
  new G4PVPlacement(RotateForUp, positionDetectorUp, logicAlNose, 
                    "physiAlNoseUp", logicWorld, false, 0);
  new G4PVPlacement(RotateForDown, positionDetectorDown, logicAlNose, 
                    "physiAlNoseDown", logicWorld, false, 0);
  new G4PVPlacement(RotateForLeft, positionDetectorLeft, logicAlNose, 
                    "physiAlNoseLeft", logicWorld, false, 0);
  new G4PVPlacement(RotateForRight, positionDetectorRight, logicAlNose, 
                    "physiAlNoseRight", logicWorld, false, 0);
  
  //------------------------------
  // Lead Detector Cap
  //------------------------------

  // Make all the half-lengths necessary  
  G4double PbNormSize= 0.5*fPbNormSize;
  G4double PbBigSize = 0.5*fPbBigSize;
  G4double CutoutWidth = 0.69*inch; 			// Half-width of cutout solid
  G4double CutoutDepth = 0.51*inch; 			// Depth of cutout
  G4double CutoutLength = 1.60*inch;			// Half-length of cutout solid
  
  // Position it behind the Aluminum nose
  // This position refers to the center of the cap
  G4double PbCapZpos = AlNosePlacement + AlNose1Length + PbNormSize; 
  
  G4ThreeVector positionPbCapUp = G4ThreeVector(0,0, PbCapZpos);
                positionPbCapUp.rotateX(-ScatteringAngle);
  G4ThreeVector positionPbCapDown = G4ThreeVector(0,0, PbCapZpos);
                positionPbCapDown.rotateX(ScatteringAngle);
  G4ThreeVector positionPbCapLeft = G4ThreeVector(0,0, PbCapZpos);
                positionPbCapLeft.rotateY(ScatteringAngle);
  G4ThreeVector positionPbCapRight = G4ThreeVector(0,0, PbCapZpos);
                positionPbCapRight.rotateY(-ScatteringAngle); 
                                              
  // Basis: 2" long tube with most common OD
  G4Tubs* solidPbNorm = new G4Tubs("solidPbNorm", 
                                   PbInRadius,
                                   PbOutRadius,
                                   PbNormSize,
                                   0.0*deg,
                                   360.0*deg);
  
  // Hoop to go around the above volume. To be Unioned                                 
  G4Tubs* solidPbBig = new G4Tubs("solidPbBig", 
                                   PbOutRadius,
                                   PbBigRadius,
                                   PbBigSize,
                                   0.0*deg,
                                   360.0*deg);
  
  // Distance along the thing the hoop needs to go
  unionShift = G4ThreeVector(0, 0, PbBigSize+PbBigPos-PbNormSize);
  
  // Make the first union. 
  G4UnionSolid* solidPbOne = new G4UnionSolid("solidPbOne",
                                              solidPbNorm,
                                              solidPbBig,
                                              0,
                                              unionShift);
  
  // Cutout hoop from lead
  G4double PbShimOR = 1.430*inch;			// Necessary for Boolean Subtraction
  G4Tubs* solidPbShim = new G4Tubs("solidPbShim", 
                                   PbSmallRadius,
                                   PbShimOR,
                                   PbSmallSize,
                                   0.0*deg,
                                   360.0*deg);
  
  // Only cut off the trailing edge
  unionShift = G4ThreeVector(0, 0, PbNormSize);
  
  G4SubtractionSolid* solidPbTwo = new G4SubtractionSolid("solidPbTwo",
                                                          solidPbOne,
                                              		  solidPbShim,
                                              		  0,
		                                          unionShift);    
  
  // Rectangle to remove for dE detector 
  G4Box* solidCutout = new G4Box("solidCutout", CutoutWidth, CutoutLength, CutoutDepth);
  unionShift = G4ThreeVector(0, 1.0*inch, PbNormSize);
  G4SubtractionSolid* solidPbCap = new G4SubtractionSolid("solidPbCap",
                                                          solidPbTwo,
                                              		  solidCutout,
                                              		  0,
		                                          unionShift);
  
  // Make the thing all ready
  G4LogicalVolume* logicPbCap = new G4LogicalVolume(solidPbCap, Pb, "logicPbCap");
  new G4PVPlacement(RotateForUp, positionPbCapUp, logicPbCap, "physiPbCapUp", logicWorld, false, 0);
  new G4PVPlacement(RotateForDown, positionPbCapDown, logicPbCap, "physiPbCapDown", logicWorld, false, 0);
  new G4PVPlacement(RotateForLeft, positionPbCapLeft, logicPbCap, "physiPbCapDown", logicWorld, false, 0);
  new G4PVPlacement(RotateForRight, positionPbCapRight, logicPbCap, "physiPbCapRight", logicWorld, false, 0);  		      		    
  
  //------------------------------------------------
  // Detector Housing. 
  //------------------------------------------------
    
  G4double HouseFrontFaceLength = 0.5*fHouseFrontFaceLength;
  G4double HouseFrontLength = 0.5*fHouseFrontLength;
  G4double HouseWeldLength = 0.5*fHouseWeldLength;
  G4double HouseCanLength = 0.5*fHouseCanLength;
  
  // Position the front edge of the housing relative to the center of the lead.    
  G4double HouseZPos = PbCapZpos - HouseFrontFaceLength;  
  G4ThreeVector positionHouseUp = G4ThreeVector(0,0,HouseZPos);
                positionHouseUp.rotateX(-ScatteringAngle);
  G4ThreeVector positionHouseDown = G4ThreeVector(0,0,HouseZPos);
                positionHouseDown.rotateX(ScatteringAngle);
  G4ThreeVector positionHouseLeft = G4ThreeVector(0,0,HouseZPos);
                positionHouseLeft.rotateY(ScatteringAngle);
  G4ThreeVector positionHouseRight = G4ThreeVector(0,0,HouseZPos);
                positionHouseRight.rotateY(-ScatteringAngle);                
  
  // Make each solid peice individually
  G4Tubs* solidHouseFrontFace = new G4Tubs("solidHouseFrontFace",
                                           HouseFrontIR,
                                           HouseFrontFaceOR,
                                           HouseFrontFaceLength,
                                           0.0*deg,
                                           360.0*deg);
  G4Tubs* solidHouseFront = new G4Tubs("solidHouseFront",
                                       HouseFrontIR,
                                       HouseFrontOR,
                                       HouseFrontLength,
                                       0.0*deg,
                                       360.0*deg);                                         
  G4Tubs* solidHouseWeld = new G4Tubs("solidHouseWeld",
                                      HouseFrontIR,
                                      HouseCanOR,
                                      HouseWeldLength,
                                      0.0*deg,
                                      360.0*deg);
  G4Tubs* solidHouseCan = new G4Tubs("solidHouseCan",
                                     HouseCanIR,
                                     HouseCanOR,
                                     HouseCanLength,
                                     0.0*deg,
                                     360.0*deg);
  
  // Combine them in order, End on end. 
  unionShift = G4ThreeVector(0, 0, HouseFrontFaceLength + HouseFrontLength);
  G4UnionSolid* HousingOne = new G4UnionSolid("HousingOne",
  					      solidHouseFrontFace,
  					      solidHouseFront,
  					      0,
  					      unionShift);  
  unionShift = G4ThreeVector(0, 0, HouseFrontFaceLength + fHouseFrontLength + HouseWeldLength);
  G4UnionSolid* HousingTwo = new G4UnionSolid("HousingTwo",
  					      HousingOne,
  					      solidHouseWeld,
  					      0,
  					      unionShift);
  unionShift = G4ThreeVector(0, 0, HouseFrontFaceLength + fHouseFrontLength 
                                   + fHouseWeldLength + HouseCanLength);
  G4UnionSolid* solidDetectorHousing = new G4UnionSolid("solidDetectorHousing",
  					                 HousingTwo,
  					                 solidHouseCan,
  					      		 0,
  					      		 unionShift);

  // Make and place them. 
  G4LogicalVolume* logicDetectorHousing = new G4LogicalVolume(solidDetectorHousing, 
  							      Steel, "logicDetectorHousing");
  new G4PVPlacement(RotateForUp, positionHouseUp, logicDetectorHousing, 
  		    "physiDetectorHousingUp", logicWorld, false, 0);
  new G4PVPlacement(RotateForDown, positionHouseDown, logicDetectorHousing, 
  		    "physiDetectorHousingDown", logicWorld, false, 0);
  new G4PVPlacement(RotateForLeft, positionHouseLeft, logicDetectorHousing, 
  		    "physiDetectorHousingLeft", logicWorld, false, 0);
  new G4PVPlacement(RotateForRight, positionHouseRight, logicDetectorHousing, 
  		    "physiDetectorHousingRight", logicWorld, false, 0);
  
  //------------------------------ 
  // DeltaE Detector
  //------------------------------
  
  G4double DeltaElength = 0.5*fDeltaElength;
  G4double DeltaEwidth = 0.5*fDeltaEwidth;
  G4double DeltaEZpos = PbCapZpos + 0.625*inch;			// 5/8" behind center of Pb. 
  
  // Center of dE is 5/8" behind the center of the Lead.
  G4ThreeVector positionDeltaEUp = G4ThreeVector(0,0,DeltaEZpos);	
                positionDeltaEUp.rotateX(-ScatteringAngle);
  G4ThreeVector positionDeltaEDown = G4ThreeVector(0,0,DeltaEZpos);
  		positionDeltaEDown.rotateX(ScatteringAngle);
  G4ThreeVector positionDeltaELeft = G4ThreeVector(0,0,DeltaEZpos);	
                positionDeltaELeft.rotateY(ScatteringAngle);
  G4ThreeVector positionDeltaERight = G4ThreeVector(0,0,DeltaEZpos);
  		positionDeltaERight.rotateY(-ScatteringAngle);  		
  
  G4Box* solidDeltaE = new G4Box("solidDeltaE", DeltaEwidth, DeltaEwidth, DeltaElength);

  // Up
  G4LogicalVolume* logicDeltaEUp = new G4LogicalVolume(solidDeltaE, DeltaEMater, "logicDeltaEUp");
  new G4PVPlacement(RotateForUp, positionDeltaEUp, logicDeltaEUp, 
  		    "physiDeltaEUp", logicWorld, false, 0);
  
  // Down
  G4LogicalVolume* logicDeltaEDown = new G4LogicalVolume(solidDeltaE, DeltaEMater, "logicDeltaEDown");
  new G4PVPlacement(RotateForDown, positionDeltaEDown, logicDeltaEDown, 
  		    "physiDeltaEDown", logicWorld, false, 0);
  		    
  // Left
  G4LogicalVolume* logicDeltaELeft = new G4LogicalVolume(solidDeltaE, DeltaEMater, "logicDeltaELeft");
  new G4PVPlacement(RotateForLeft, positionDeltaELeft, logicDeltaELeft, 
  		    "physiDeltaELeft", logicWorld, false, 0);
  
  // Right
  G4LogicalVolume* logicDeltaERight = new G4LogicalVolume(solidDeltaE,DeltaEMater,"logicDeltaERight");
  new G4PVPlacement(RotateForRight, positionDeltaERight, logicDeltaERight, 
  		    "physiDeltaERight", logicWorld, false, 0);  		    
  
  //------------------------------ 
  // Energy Detector
  //------------------------------

  G4double EDetLength = 0.5*fEDetLength; 
  G4double EDetZpos = PbCapZpos + PbNormSize + EDetLength + 4.0*mm; // against the Housing	
  
  G4ThreeVector positionEDetUp = G4ThreeVector(0,0,EDetZpos);
  		positionEDetUp.rotateX(-ScatteringAngle);
  G4ThreeVector positionEDetDown = G4ThreeVector(0,0,EDetZpos);
  		positionEDetDown.rotateX(ScatteringAngle);
  G4ThreeVector positionEDetLeft = G4ThreeVector(0,0,EDetZpos);
  		positionEDetLeft.rotateY(ScatteringAngle);
  G4ThreeVector positionEDetRight = G4ThreeVector(0,0,EDetZpos);
  		positionEDetRight.rotateY(-ScatteringAngle);
  
  G4Tubs* solidEDet = new G4Tubs("solidEDet", 
  				 0,
  				 EDetRadius,
  				 EDetLength,
  				 0.0*deg,
  				 360.0*deg);
  // Up
  G4LogicalVolume* logicEDetUp = new G4LogicalVolume(solidEDet, EDetMater, "logicEDetUp");
  new G4PVPlacement(RotateForUp, positionEDetUp, logicEDetUp, "physiEDetUp", logicWorld, false, 0);
  
  // Down
  G4LogicalVolume* logicEDetDown = new G4LogicalVolume(solidEDet, EDetMater, "logicEDetDown");
  new G4PVPlacement(RotateForDown, positionEDetDown, logicEDetDown, "physiEDetDown", logicWorld, false, 0);
                    
  // Left
  G4LogicalVolume* logicEDetLeft = new G4LogicalVolume(solidEDet, EDetMater, "logicEDetLeft");
  new G4PVPlacement(RotateForLeft, positionEDetLeft, logicEDetLeft, "physiEDetLeft", logicWorld, false, 0);
  
  // Right
  G4LogicalVolume* logicEDetRight = new G4LogicalVolume(solidEDet, EDetMater, "logicEDetRight");
  new G4PVPlacement(RotateForRight, positionEDetRight, logicEDetRight, "physiEDetRight", logicWorld, false, 0);                    

  //------------------------------------------------		    
  //  E "Photocathode"
  //------------------------------------------------
  
  G4double EPhotoLength = 0.5*fEPhotoLength;
  
  G4double EPhotoZpos = EDetZpos + EDetLength + EPhotoLength;
  
  G4ThreeVector positionEPhotoUp = G4ThreeVector(0,0,EPhotoZpos);
                positionEPhotoUp.rotateX(-ScatteringAngle);
  G4ThreeVector positionEPhotoDown = G4ThreeVector(0,0,EPhotoZpos);
                positionEPhotoDown.rotateX(ScatteringAngle);
  G4ThreeVector positionEPhotoLeft = G4ThreeVector(0,0,EPhotoZpos);
                positionEPhotoLeft.rotateY(ScatteringAngle);
  G4ThreeVector positionEPhotoRight = G4ThreeVector(0,0,EPhotoZpos);
                positionEPhotoRight.rotateY(-ScatteringAngle);
                
  G4Tubs* solidEPhotoPlate = new G4Tubs("solidEPhotoPlate", 0.0*mm, EPhotoRadius, EPhotoLength, 0.0*deg, 360.0*deg);
  G4LogicalVolume* logicEPhotoPlate = new G4LogicalVolume(solidEPhotoPlate, PhotoCathodeMaterial, "logicEPhotoPlate");

  new G4PVPlacement(RotateForUp, positionEPhotoUp, logicEPhotoPlate, "physiEPhotoUp", logicWorld, false, 0);
  new G4PVPlacement(RotateForDown, positionEPhotoDown, logicEPhotoPlate, "physiEPhotoDown", logicWorld, false, 0);
  new G4PVPlacement(RotateForLeft, positionEPhotoLeft, logicEPhotoPlate, "physiEPhotoLeft", logicWorld, false, 0);
  new G4PVPlacement(RotateForRight, positionEPhotoRight, logicEPhotoPlate, "physiEPhotoRight", logicWorld, false, 0);
  
  //------------------------------------------------		    
  //  dE "Photocathode"
  //------------------------------------------------
  
  G4double dEPhotoLength = 0.5*fdEPhotoLength;
  G4double dEPhotoWidth = 0.5*fdEPhotoWidth;
  
  G4double dEPhotoZpos = DeltaEZpos;
  G4double dEPhotoYpos = dEPhotoLength + DeltaEwidth; 
  
  G4ThreeVector positiondEPhotoUp = G4ThreeVector(0,dEPhotoYpos,dEPhotoZpos);
                positiondEPhotoUp.rotateX(-ScatteringAngle);
  G4ThreeVector positiondEPhotoDown = G4ThreeVector(0,dEPhotoYpos,dEPhotoZpos);
                positiondEPhotoDown.rotateX(ScatteringAngle);
  G4ThreeVector positiondEPhotoLeft = G4ThreeVector(0,dEPhotoYpos,dEPhotoZpos);
                positiondEPhotoLeft.rotateY(ScatteringAngle);
  G4ThreeVector positiondEPhotoRight = G4ThreeVector(0,dEPhotoYpos,dEPhotoZpos);
                positiondEPhotoRight.rotateY(-ScatteringAngle);
                
  G4Box* solidDeltaEPhotoPlate = new G4Box("solidDeltaEPhotoPlate", dEPhotoWidth, dEPhotoLength, dEPhotoLength);
  G4LogicalVolume* logicdEPhotoPlate = new G4LogicalVolume(solidDeltaEPhotoPlate, PhotoCathodeMaterial, "logicDeltaEPhotoPlate");  
  
  new G4PVPlacement(RotateForUp, positiondEPhotoUp, logicdEPhotoPlate, "physidEPhotoUp", logicWorld, false, 0);
  new G4PVPlacement(RotateForDown, positiondEPhotoDown, logicdEPhotoPlate, "physidEPhotoDown", logicWorld, false, 0);
  new G4PVPlacement(RotateForLeft, positiondEPhotoLeft, logicdEPhotoPlate, "physidEPhotoLeft", logicWorld, false, 0);
  new G4PVPlacement(RotateForRight, positiondEPhotoRight, logicdEPhotoPlate, "physidEPhotoRight", logicWorld, false, 0);
  
  //------------------------------------------------		    
  // Air Column
  //------------------------------------------------
  
  G4double AirLength = 0.5*fAirLength;
  G4double AirZpos = DeltaEZpos - (AirLength+DeltaElength);
  
  G4ThreeVector positionAirUp = G4ThreeVector(0,0,AirZpos);
                positionAirUp.rotateX(-ScatteringAngle);
  G4ThreeVector positionAirDown = G4ThreeVector(0,0,AirZpos);
                positionAirDown.rotateX(ScatteringAngle);
  G4ThreeVector positionAirLeft = G4ThreeVector(0,0,AirZpos);
                positionAirLeft.rotateY(ScatteringAngle);
  G4ThreeVector positionAirRight = G4ThreeVector(0,0,AirZpos);
                positionAirRight.rotateY(-ScatteringAngle);
  
  G4Tubs* solidAir = new G4Tubs("solidAir", 0.0*mm, AirRadius, AirLength, 0.0*deg, 360.0*deg);  				
  G4LogicalVolume* logicAirColumn = new G4LogicalVolume(solidAir, Air, "logicAirColumn");
  
  new G4PVPlacement(RotateForUp, positionAirUp, logicAirColumn, "physiAirColumnUp", logicWorld, false, 0);
  new G4PVPlacement(RotateForDown, positionAirDown, logicAirColumn, "physiAirColumnDown", logicWorld, false, 0);
  new G4PVPlacement(RotateForLeft, positionAirLeft, logicAirColumn, "physiAirColumnLeft", logicWorld, false, 0);
  new G4PVPlacement(RotateForRight, positionAirRight, logicAirColumn, "physiAirColumnRight", logicWorld, false, 0);
  
  //------------------------------------------------		    
  // Aluminium Window
  //------------------------------------------------
  
  G4double AlWindowLength = 4.0e-03*inch;
  G4double AlWindowRadius = 1.0*inch;
  G4double AlWindowZPos = AirZpos - (AirLength+AlWindowLength);
  
  G4ThreeVector positionAlWindowUp = G4ThreeVector(0,0,AlWindowZPos);
                positionAlWindowUp.rotateX(-ScatteringAngle);
  G4ThreeVector positionAlWindowDown = G4ThreeVector(0,0,AlWindowZPos);
                positionAlWindowDown.rotateX(ScatteringAngle);
  G4ThreeVector positionAlWindowLeft = G4ThreeVector(0,0,AlWindowZPos);
                positionAlWindowLeft.rotateY(ScatteringAngle);
  G4ThreeVector positionAlWindowRight = G4ThreeVector(0,0,AlWindowZPos);
                positionAlWindowRight.rotateY(-ScatteringAngle);
  
  
  G4Tubs* solidAlWindow = new G4Tubs("solidAlWindow", 0.0*mm, AlWindowRadius, AlWindowLength, 0.0*deg, 360.0*deg);
  G4LogicalVolume* logicAlWindow = new G4LogicalVolume(solidAlWindow, Al, "logicAlWindow");

  new G4PVPlacement(RotateForUp, positionAlWindowUp, logicAlWindow, "physiAlWindowUp", logicWorld, false, 0);
  new G4PVPlacement(RotateForDown, positionAlWindowDown, logicAlWindow, "physiAlWindowDown", logicWorld, false, 0);
  new G4PVPlacement(RotateForLeft, positionAlWindowLeft, logicAlWindow, "physiAlWindowLeft", logicWorld, false, 0);
  new G4PVPlacement(RotateForRight, positionAlWindowRight, logicAlWindow, "physiAlWindowRight", logicWorld, false, 0);
  
//////////////////////////////////////////////////////////////////////////
///////////////////////// OPTICAL PROPERTIES /////////////////////////////
//////////////////////////////////////////////////////////////////////////
  
  // Spectrum for which I have data on the scintillators. 
  const G4int nentries = 14;
  G4double PhotonEnergy[nentries] = {2.48*eV, 	// 500 nm
  				     2.53*eV,	// 490 nm
  				     2.58*eV,   // 480 nm
  				     2.64*eV,
  				     2.70*eV,
  				     2.76*eV,
  				     2.82*eV,
  				     2.88*eV,   
  				     2.92*eV,  	// 425 nm (Max)
  				     2.95*eV,
  				     3.02*eV,
  				     3.10*eV,
  				     3.18*eV,   // 390 nm
  				     3.26*eV}; 	// 380 nm

  //----------------------------------
  // EJ-200/EJ-212 (It's not really either).
  //----------------------------------

  // Scintillation Optical properties
  G4double ScintillationSpectrum[nentries] = {0.06, 0.11, 0.18, 0.30, 0.42, 0.51, 
					      0.74, 0.92, 1.00, 0.86, 0.30, 0.02, 
					      0.00, 0.00};
  // Cerenkov properties
  G4double RefractiveIndex[nentries] = {1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 
					1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 
					1.58, 1.58};
  G4double AttenuationLength[nentries] = {250*cm, 250*cm, 250*cm, 250*cm, 250*cm, 
					  250*cm, 250*cm, 250*cm, 250*cm, 250*cm, 
					  250*cm, 250*cm, 250*cm, 250*cm};

  G4MaterialPropertiesTable* PlasticMPT = new G4MaterialPropertiesTable();  
  PlasticMPT->AddConstProperty("SCINTILLATIONYIELD", 10000.0/MeV);
  PlasticMPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  PlasticMPT->AddConstProperty("FASTTIMECONSTANT", 2.4*ns);
  PlasticMPT->AddConstProperty("YIELDRATIO", 1.0);
  PlasticMPT->AddProperty("FASTCOMPONENT", PhotonEnergy, ScintillationSpectrum, nentries);
  PlasticMPT->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex, nentries); 
  PlasticMPT->AddProperty("ABSLENGTH", PhotonEnergy, AttenuationLength, nentries);

  Plastic->SetMaterialPropertiesTable(PlasticMPT);
  
  //-------------------------------------
  // Air/Vacuum Optical properties (identical)
  //-------------------------------------

  G4double AirRefractiveIndex[nentries] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					   0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
					   0.0, 0.0};
  G4double AirAttenuationLength[nentries] = {1000*km,1000*km,1000*km,1000*km,1000*km,1000*km,
					     1000*km,1000*km,1000*km,1000*km,1000*km,1000*km,
					     1000*km,1000*km};
  
  G4MaterialPropertiesTable* AirMPT = new G4MaterialPropertiesTable();
  AirMPT->AddProperty("RINDEX", PhotonEnergy, AirRefractiveIndex, nentries);
  AirMPT->AddProperty("ABSLENGTH", PhotonEnergy, AirAttenuationLength, nentries);

  Air->SetMaterialPropertiesTable(AirMPT);
  WorldMater->SetMaterialPropertiesTable(AirMPT);

  //-------------------------------------
  // Photocathode MPT
  //-------------------------------------
  
  
  G4double PMTQuantumEfficiency[nentries] = {	0.103,  // 500 nm
 						0.105,  // 490 nm
						0.107,	// 480 nm
						0.108,
						0.110,	
						0.110,
						0.111,
						0.111,
						0.112,  // 425 nm
						0.112,
						0.113,
						0.113, 
						0.113,  // 390 nm
						0.112 };// 380 nm 
  /*
  G4double PMTQuantumEfficiency[nentries] = {	1.00,
  						1.00,
  						1.00,
  						1.00,
  						1.00,
  						1.00,
  						1.00,
  						1.00,
  						1.00,
  						1.00,
  						1.00,
  						1.00,
  						1.00,
  						1.00  };// 380 nm 
  */
  G4double PMTReflectivity[nentries] = { 0.25,
					 0.25,
					 0.25,
					 0.25,
					 0.25,
					 0.25,
					 0.25,
					 0.25,
					 0.25,
					 0.25,
					 0.25,
					 0.25,
					 0.25,
					 0.25 }; 

  G4MaterialPropertiesTable* PhotocathodeMPT = new G4MaterialPropertiesTable();
  PhotocathodeMPT->AddProperty("REFLECTIVITY", PhotonEnergy, PMTReflectivity,nentries);
  PhotocathodeMPT->AddProperty("EFFICIENCY", PhotonEnergy, PMTQuantumEfficiency,nentries);

  G4OpticalSurface* PhotocathodeOpticalSurface =  new G4OpticalSurface("PhotocathodeOS");
  PhotocathodeOpticalSurface->SetType(dielectric_metal); 
  PhotocathodeOpticalSurface->SetFinish(polished); 
  PhotocathodeOpticalSurface->SetModel(glisur);
  PhotocathodeOpticalSurface->SetMaterialPropertiesTable(PhotocathodeMPT);

  new G4LogicalSkinSurface("E_PMT_Surface",logicEPhotoPlate,PhotocathodeOpticalSurface);
  new G4LogicalSkinSurface("dE_PMT_Surface",logicdEPhotoPlate,PhotocathodeOpticalSurface);

//////////////////////////////////////////////////////////////////////////
/////////////////////// VISUALIZATION ATTRIBUTES /////////////////////////
//////////////////////////////////////////////////////////////////////////

  // World = White Wire Cube
  G4VisAttributes* WorldVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  logicWorld->SetVisAttributes(WorldVisAtt);

  // Target = Magenta Disk
  G4VisAttributes* targetVisAtt = new G4VisAttributes(G4Colour::Magenta());
  targetVisAtt->SetForceSolid(true);
  logicTarget->SetVisAttributes(targetVisAtt);
  // Also the "Photocathodes"
  logicdEPhotoPlate->SetVisAttributes(targetVisAtt);
  logicEPhotoPlate->SetVisAttributes(targetVisAtt);
  
  // Energy Detector = Solid Cyan Cylinders
  G4VisAttributes* EDetVisAtt = new G4VisAttributes(G4Colour::Cyan());
  EDetVisAtt->SetForceSolid(true);
  logicEDetUp->SetVisAttributes(EDetVisAtt);
  logicEDetDown->SetVisAttributes(EDetVisAtt);
  logicEDetLeft->SetVisAttributes(EDetVisAtt);
  logicEDetRight->SetVisAttributes(EDetVisAtt);
  
  // DeltaE detectors = Solid Blue Squares
  G4VisAttributes* DeltaEVisAtt = new G4VisAttributes(G4Colour::Blue());
  DeltaEVisAtt->SetForceSolid(true);
  logicDeltaEUp->SetVisAttributes(DeltaEVisAtt);
  logicDeltaEDown->SetVisAttributes(DeltaEVisAtt);
  logicDeltaELeft->SetVisAttributes(DeltaEVisAtt);
  logicDeltaERight->SetVisAttributes(DeltaEVisAtt);
  
  // Aluminum Components = Solid Yellow
  G4VisAttributes* AlVisAtt = new G4VisAttributes(G4Colour::Yellow());
  AlVisAtt->SetForceSolid(true);
  logicAlNose->SetVisAttributes(AlVisAtt);
  logicChamberLiner->SetVisAttributes(AlVisAtt);
  logicAlWindow->SetVisAttributes(AlVisAtt);
  logicAlBaffle->SetVisAttributes(AlVisAtt);
  logicExtensionTube->SetVisAttributes(AlVisAtt);
  logicBeDumpSegment->SetVisAttributes(AlVisAtt);
  
  // Stainless Steel Components = Solid White (appears light grey in OGL)
  G4VisAttributes* SteelVisAtt = new G4VisAttributes(G4Colour::White());
  SteelVisAtt->SetForceSolid(true);
  logicDetectorHousing->SetVisAttributes(SteelVisAtt);
  logicMainChamber->SetVisAttributes(SteelVisAtt);
  logicEndFlange->SetVisAttributes(SteelVisAtt);
  logicPumpChamber->SetVisAttributes(SteelVisAtt);
  
  // Lead Components = Solid Grey
  G4VisAttributes* PbVisAtt = new G4VisAttributes(G4Colour::Grey());
  PbVisAtt->SetForceSolid(true);
  logicPbCap->SetVisAttributes(PbVisAtt);
  
  // Copper Components = "Orange" 
  G4Colour orange (1.0, 0.5, 0.0);
  G4VisAttributes* CuVisAtt = new G4VisAttributes(orange);
  CuVisAtt->SetForceSolid(true);
  logicCuDumpSegment->SetVisAttributes(CuVisAtt);
  
  // Invisible Components
  G4VisAttributes* InVisAtt = new G4VisAttributes(false);
  logicCuDump->SetVisAttributes(InVisAtt);
  logicBeDump->SetVisAttributes(InVisAtt);
  
///////////////////////////////////////////////////////////////////////////
//////////////////////// SENSITIVE DETECTOR(S) ////////////////////////////
///////////////////////////////////////////////////////////////////////////
         
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  // Up
  G4String dE_1_name = "Mott/dE_1";
  MottTrackerSD* dE_1_SD = new MottTrackerSD( dE_1_name );
  SDman->AddNewDetector( dE_1_SD );
  logicDeltaEUp->SetSensitiveDetector( dE_1_SD );
  G4String E_1_name = "Mott/E_1";
  MottTrackerSD* E_1_SD = new MottTrackerSD( E_1_name );
  SDman->AddNewDetector( E_1_SD );
  logicEDetUp->SetSensitiveDetector( E_1_SD );
  
  // Down
  G4String dE_2_name = "Mott/dE_2";
  MottTrackerSD* dE_2_SD = new MottTrackerSD( dE_2_name );
  SDman->AddNewDetector( dE_2_SD );
  logicDeltaEDown->SetSensitiveDetector( dE_2_SD );
  G4String E_2_name = "Mott/E_2";
  MottTrackerSD* E_2_SD = new MottTrackerSD( E_2_name );
  SDman->AddNewDetector( E_2_SD );
  logicEDetDown->SetSensitiveDetector( E_2_SD );
  
  // Left
  G4String dE_3_name = "Mott/dE_3";
  MottTrackerSD* dE_3_SD = new MottTrackerSD( dE_3_name );
  SDman->AddNewDetector( dE_3_SD );
  logicDeltaELeft->SetSensitiveDetector( dE_3_SD );
  G4String E_3_name = "Mott/E_3";
  MottTrackerSD* E_3_SD = new MottTrackerSD( E_3_name );
  SDman->AddNewDetector( E_3_SD );
  logicEDetLeft->SetSensitiveDetector( E_3_SD );

  // Right
  G4String dE_4_name = "Mott/dE_4";
  MottTrackerSD* dE_4_SD = new MottTrackerSD( dE_4_name );
  SDman->AddNewDetector( dE_4_SD );
  logicDeltaERight->SetSensitiveDetector( dE_4_SD );
  G4String E_4_name = "Mott/E_4";
  MottTrackerSD* E_4_SD = new MottTrackerSD( E_4_name );
  SDman->AddNewDetector( E_4_SD );
  logicEDetRight->SetSensitiveDetector( E_4_SD );
  
  /* Dump Plate(s)
  // Beryllium
  G4String BeDumpSDName = "Mott/BeDump";
  MottTrackerSD* BeDumpSD = new MottTrackerSD(BeDumpSDName);
  SDman->AddNewDetector(BeDumpSD);
  logicBeDumpSegment->SetSensitiveDetector(BeDumpSD);
  
  // Copper
  G4String CuDumpSDName = "Mott/CuDump";
  MottTrackerSD* CuDumpSD = new MottTrackerSD(CuDumpSDName);
  SDman->AddNewDetector(CuDumpSD);
  logicCuDumpSegment->SetSensitiveDetector(CuDumpSD);
  */
  
///////////////////////////////////////////////////////////////////////////
////////////////////////// USER LIMITS(S) /////////////////////////////////
///////////////////////////////////////////////////////////////////////////

  // Sets a max Step length in the world region, with G4StepLimiter
  G4double maxStep = 0.5*HalfWorldLength;
  stepLimit = new G4UserLimits(maxStep);
  logicWorld->SetUserLimits(stepLimit);
  
  std::cout << "\tLeaving MottDetectorConstruction::Construct()" << std::endl;

  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MottDetectorConstruction::SetTargetMater(G4String newMater) {

  G4double a, z, density;

  G4Material* Au = new G4Material("Gold", z=79, a=196.97*g/mole, density=19.32*g/cm3);
  G4Material* Ag = new G4Material("Silver", z=47, a=107.8682*g/mole, density=10.49*g/cm3);

  if(newMater=="Au") {
    TargetMater = Au;
  } else if(newMater=="Ag") {
    TargetMater = Ag;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MottDetectorConstruction::SetTargetIn()
{
  G4cout << " ##### Calling MottDetectorConstruction::SetTargetIn() ##### " << G4endl;
 
  TargetXpos = 0.0*mm;		//
  TargetYpos = 0.0*mm;		// Target home position
  TargetZpos = -3.9878*mm;	//
  
  physiTarget->SetTranslation(G4ThreeVector(TargetXpos,TargetYpos,TargetZpos));
  
  G4cout << " ##### Leaving MottDetectorConstruction::SetTargetIn() ##### " << G4endl;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MottDetectorConstruction::SetTargetOut()
{
  G4cout << " ##### Calling MottDetectorConstruction::SetTargetOut() ##### " << G4endl;

  TargetXpos = 0.0*mm;		//
  TargetYpos = 1000.0*mm;	// Target away position
  TargetZpos = -3.9878*mm;	//
  
  physiTarget->SetTranslation(G4ThreeVector(TargetXpos,TargetYpos,TargetZpos));

  G4cout << " ##### Leaving MottDetectorConstruction::SetTargetOut() ##### " << G4endl;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void MottDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((stepLimit)&&(maxStep>0.)) stepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
