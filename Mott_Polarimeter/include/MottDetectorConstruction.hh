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
// $Id: MottDetectorConstruction.hh,v 3.6 2014/01/16 23:40:03 mjmchugh Exp mjmchugh $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef MottDetectorConstruction_h
#define MottDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UserLimits;
class MottDetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MottDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
     MottDetectorConstruction();
    ~MottDetectorConstruction();

  public:
  
     G4VPhysicalVolume* Construct();
     
     G4double GetTargetFullLength()  {return fTargetLength;};
     G4double GetWorldFullLength()   {return fWorldLength;}; 
     
     void SetTargetLength(G4double);
     void SetTargetIn();
     void SetTargetOut();
     
     void SetMaxStep (G4double);     
     
  private:

     // World Volume
     G4Box*             solidWorld;		// pointer to the solid envelope 
     G4LogicalVolume*   logicWorld;		// pointer to the logical envelope
     G4VPhysicalVolume* physiWorld;		// pointer to the physical envelope
    
     G4Material*        WorldMater;		// pointer to the chamber material 
     G4double fWorldLength;			// Full length of the world volume
     
     // Target Volume
     G4Tubs*            solidTarget;		// pointer to the solid Target
     G4LogicalVolume*   logicTarget;		// pointer to the logical Target
     G4VPhysicalVolume* physiTarget;		// pointer to the physical Target

     G4Material*        TargetMater;		// pointer to the target  material
     G4double fTargetLength;			// Full length of Target
     G4double TargetZpos;			//
     G4double TargetXpos;			// Position
     G4double TargetYpos;   			//
        
     G4UserLimits* stepLimit;			// pointer to user step limits
     
     MottDetectorMessenger* detectorMessenger;	// pointer to the Messenger
       

     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
