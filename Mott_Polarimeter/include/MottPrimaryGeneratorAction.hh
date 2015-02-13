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
// $Id: MottPrimaryGeneratorAction.hh,v 3.3 2014/01/16 23:40:03 mjmchugh Exp mjmchugh $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#ifndef MottPrimaryGeneratorAction_h
#define MottPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include <vector>

class MottPrimaryGeneratorMessenger;
class MottDetectorConstruction;
class MottEventAction;
class G4ParticleGun;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
class MottPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    MottPrimaryGeneratorAction(MottDetectorConstruction*);    
   ~MottPrimaryGeneratorAction();

    void GeneratePrimaries(G4Event*);

    void SetBeamEnergy(G4double);
    G4double GetBeamEnergy() { return beamEnergy; };
    
    void SetEnergySpread(G4double spread) { energySpread = spread; };
    G4double GetEnergySpread() { return energySpread; };
    
    void SetBeamDiameter(G4double diameter) { beamDiameter = diameter; };
    G4double GetBeamDiameter() { return beamDiameter; };

    void ReadDataFiles();

    //Linear interpolation functions for arrays.
    G4double InterpolateCrossSection(G4double,G4double);					
    G4double InterpolateSherman(G4double,G4double);
    G4double InterpolateT(G4double,G4double);
    G4double InterpolateU(G4double,G4double);
    
  private:
  
    G4ParticleGun* particleGun;
    MottDetectorConstruction* myDetector;
    MottPrimaryGeneratorMessenger* myMessenger;
    MottEventAction* pEventAction;
    
    // Target Properties
    G4double TargetZ;

    // Beam Properties
    G4double beamEnergy;
    G4double energySpread;	// Std. Dev.
    G4double beamDiameter;  	// FWHM 
    
    std::vector < std::vector <G4double> > ThetaSc;
    std::vector < std::vector <G4double> > CrossSection;
    std::vector < std::vector <G4double> > SpinT;
    std::vector < std::vector <G4double> > SpinU;
    std::vector < std::vector <G4double> > Sherman;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


