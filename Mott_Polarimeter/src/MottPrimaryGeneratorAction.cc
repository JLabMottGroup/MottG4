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
// $Id: MottPrimaryGeneratorAction.cc,v 3.9 2014/01/16 23:40:03 mjmchugh Exp mjmchugh $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MottPrimaryGeneratorAction.hh"
#include "MottDetectorConstruction.hh"
#include "MottPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottPrimaryGeneratorAction::MottPrimaryGeneratorAction(
                                               MottDetectorConstruction* myDC)
:myDetector(myDC)
{
  G4int n_particle = 1;

  myMessenger = new MottPrimaryGeneratorMessenger(this);

  particleGun = new G4ParticleGun(n_particle);
  
  // Set default 
  beamEnergy = 5.0*MeV;
  energySpread = 5e-03*MeV;
  beamDiameter = 1.0*mm;

  // Set default particle to electron
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("e-");

  // Anything set in this member function 
  // should not be set elsewhere and will 
  // remain constant through a run.
  
  G4ThreeVector direction = G4ThreeVector(0,0,1.0);     // no dispersion
  
  // Set the constant gun properties
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(direction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottPrimaryGeneratorAction::~MottPrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // All Gun settings changed here will change every event. 

  // Source Position
  G4double sigma = beamDiameter/(2.354820045*mm);
  G4double X = G4RandGauss::shoot(0.0,sigma)*mm;
  G4double Y = G4RandGauss::shoot(0.0,sigma)*mm;
  G4double Z = -10.0*cm;
  G4ThreeVector position = G4ThreeVector(X, Y, Z);
    
  // Beam Energy 
  G4double energy = G4RandGauss::shoot(beamEnergy/(1.0*MeV), energySpread/(1.0*MeV))*MeV;
 
  // Set variable gun properties
  particleGun->SetParticleEnergy(energy);
  particleGun->SetParticlePosition(position);
  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  // Thrown Angle
  //G4double ScatteringAngle = 172.7*deg;			// Average acceptance angle.
  //G4double Theta = ScatteringAngle - 1.0*deg + 2.0*G4UniformRand()*deg;
  //G4double Phi = 10.0*G4UniformRand()*deg - 5.0*deg;  
  //G4ThreeVector direction;
  //              direction.setRThetaPhi(1.0,Theta,Phi);   
  //particleGun->SetParticleMomentumDirection(direction);  


