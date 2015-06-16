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
#include <map>

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

    void SetThetaMin(G4double theta) { ThetaMin = theta; };
    G4double GetThetaMin() { return ThetaMin; };

    void SetThetaMax(G4double theta) { ThetaMax = theta; };
    G4double GetThetaMax() { return ThetaMax; };

    void SetPhiMin(G4double phi) { PhiMin = phi; };
    G4double GetPhiMin() { return PhiMin; };

    void SetPhiMax(G4double phi) { PhiMax = phi; };
    G4double GetPhiMax() { return PhiMax; };

    void SetThrowFromUpstream() { ThrowFromUpstream = true; };
    void SetThrowAtCollimators() { ThrowFromUpstream = false; ThrowAtCollimators = true; };
    void SetThrowInUserRange() { ThrowFromUpstream = false; ThrowAtCollimators = false; };

    void ReadDataFiles();

    int CalculateNewPol();

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

    typedef std::map <G4double, G4double> TargetELoss; 
    TargetELoss GoldELoss;    

    // Target Properties
    G4double TargetZ;		// Atomic Number of Target	

    // Beam Properties
    G4double beamEnergy;
    G4double energySpread;	// Std. Dev.
    G4double beamDiameter;  	// FWHM 
    
    // Primary Event Type
    G4bool ThrowFromUpstream;
    G4bool ThrowAtCollimators;

    // Angular Ranges
    G4double ThetaMin;
    G4double ThetaMax;
    G4double PhiMin;
    G4double PhiMax;

    std::vector < std::vector <G4double> > ThetaSc;
    std::vector < std::vector <G4double> > CrossSection;
    std::vector < std::vector <G4double> > SpinT;
    std::vector < std::vector <G4double> > SpinU;
    std::vector < std::vector <G4double> > Sherman;

    //Scattering vertex quantities
    G4double X;				// Position
    G4double Y;	
    G4double Z;
    G4double Px1;			// Incoming polarization
    G4double Py1;
    G4double Pz1;
    G4double Px2;			// Outgoing polarization
    G4double Py2;
    G4double Pz2;
    G4double Theta;			// Scattering Angle
    G4double Phi;
    G4double Energy;			// Energy at the vertex
    G4double CS;			// Cross Section
    G4double S;				// Sherman Function
    G4double T;				// Spin T
    G4double U;				// Spin U

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


