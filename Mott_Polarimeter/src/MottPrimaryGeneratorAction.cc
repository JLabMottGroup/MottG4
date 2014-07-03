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
#include <iostream>
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottPrimaryGeneratorAction::MottPrimaryGeneratorAction(
                                               MottDetectorConstruction* myDC)
:myDetector(myDC)
{
  G4cout << "\tEntering MottPrimaryGeneratorAction::MottPrimaryGeneratorAction()" <<G4endl; 

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
  particleGun->SetParticleDefinition(particle);
  
  //Read in Xavi's data.
  std::ifstream inFile;
  inFile.open("/home/mjmchugh/Mott/MottG4/CrossSections/Mott-DWBA-Au-5MeV.out");
  
  G4double theta, i, t, u, s;
  for(G4int nLines=1; nLines<=427; nLines++) {
    inFile >> theta >> i >> t >> u >> s;
    ThetaSc.push_back(theta);
    CrossSection.push_back(i);
    SpinT.push_back(t);
    SpinU.push_back(u);
    Sherman.push_back(s);
    //G4cout << "\t" << nLines << "\t" << theta << "\t" << i << "\t" << t << "\t" << u << "\t" << s << G4endl;
  }
  
  inFile.close();
  
  G4cout << "\tLeaving MottPrimaryGeneratorAction::MottPrimaryGeneratorAction()" <<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottPrimaryGeneratorAction::~MottPrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // std::cout << "\tEntering MottPrimaryGeneratorAction::GeneratePrimaries()" << std::endl;

  // All Gun settings changed here will change every event. 

  // Source Position
  G4int ThrowFromUpstream = 0;
  
  G4ThreeVector gunPosition = G4ThreeVector(0,0,0);
  
  // Incident electron's polarization set here
  G4double Sx = 0.0;
  G4double Sy = 0.0;
  G4double Sz = 1.0;
  G4ThreeVector preScatteredPolarization = G4ThreeVector(Sx,Sy,Sz);

  G4double Phi = 0; 
  G4double Theta = 0;
  if(ThrowFromUpstream) {
    G4double sigma = beamDiameter/(2.354820045*mm);
    G4double X = G4RandGauss::shoot(0.0,sigma)*mm;
    G4double Y = G4RandGauss::shoot(0.0,sigma)*mm;
    G4double Z = -10.0*cm;
    gunPosition = G4ThreeVector(X, Y, Z);
  } else {
    G4int goodThrow = 0;
    while (goodThrow==0) { 
      G4double ScatteringAngle = 172.6*deg;					// Average acceptance angle.
      Theta = ScatteringAngle - 0.6*deg + 1.2*G4UniformRand()*deg;	// Throw from 172.0 to 173.2 degrees
      G4double quadrantRoll = G4UniformRand();					// Pick a quadrant							
      if(0.0<=quadrantRoll&&quadrantRoll<0.25) {				// Throw a 5 degree window in phi around each aperature
        Phi = 5.0*G4UniformRand()*deg - 2.5*deg;        
      } else if(0.25<=quadrantRoll&&quadrantRoll<0.5) {
        Phi = 5.0*G4UniformRand()*deg + 87.5*deg;
      } else if(0.5<=quadrantRoll&&quadrantRoll<0.75) {
        Phi = 5.0*G4UniformRand()*deg + 177.5*deg;
      } else if(0.75<=quadrantRoll&&quadrantRoll<=1.0) {
        Phi = 5.0*G4UniformRand()*deg + 267.5*deg;
      }
      G4double cs = InterpolateCrossSection(Theta/deg)*(1 + InterpolateSherman(Theta/deg)*cos(Phi));
      G4double rejectionThrow = 8.4602e-26*G4UniformRand();
      if(rejectionThrow<=cs) goodThrow = 1; 
    }
      
    G4ThreeVector direction;
                  direction.setRThetaPhi(1.0,Theta,Phi);   
    particleGun->SetParticleMomentumDirection(direction);
    //G4cout << "\t" << Theta/deg << "\t" << Phi/deg << G4endl;
  }
  
  // Beam Energy 
  G4double energy = G4RandGauss::shoot(beamEnergy/(1.0*MeV), energySpread/(1.0*MeV))*MeV;
 
  // Set variable gun properties
  particleGun->SetParticleEnergy(energy);
  particleGun->SetParticlePosition(gunPosition);
  particleGun->GeneratePrimaryVertex(anEvent);
  
  // Thrown Angle
  
  // std::cout << "\tLeaving MottPrimaryGeneratorAction::GeneratePrimaries()" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MottPrimaryGeneratorAction::InterpolateCrossSection(G4double x) {

  G4double y = -1;

  if(x<3.7||179.5<x) {
    G4cout << "Angle outside of range " << x << G4endl;
    return y;
  }
  
  G4int i=0;
  while(ThetaSc[i]<x) i++;
  G4double x1 = ThetaSc[i];	G4double y1 = CrossSection[i];
  G4double x0 = ThetaSc[i-1];	G4double y0 = CrossSection[i-1];
  
  y = y0 + (y1-y0)*(x-x0)/(x1-x0);

  return y;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double MottPrimaryGeneratorAction::InterpolateSherman(G4double x) {

  G4double y = -1;

  if(x<3.7||179.5<x) {
    G4cout << "Angle outside of range " << x << G4endl;
    return y;
  }
  
  G4int i=0;
  while(ThetaSc[i]<x) i++;

  G4double x1 = ThetaSc[i];	G4double y1 = Sherman[i];
  G4double x0 = ThetaSc[i-1];	G4double y0 = Sherman[i-1];
  
  y = y0 + (y1-y0)*(x-x0)/(x1-x0);

  return y;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double MottPrimaryGeneratorAction::InterpolateT(G4double x) {

  G4double y = -1;

  if(x<3.7||179.5<x) {
    G4cout << "Angle outside of range " << x << G4endl;
    return y;
  }
  
  G4int i=0;
  while(ThetaSc[i]<x) i++;

  G4double x1 = ThetaSc[i];	G4double y1 = SpinT[i];
  G4double x0 = ThetaSc[i-1];	G4double y0 = SpinT[i-1];
  
  y = y0 + (y1-y0)*(x-x0)/(x1-x0);

  return y;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double MottPrimaryGeneratorAction::InterpolateU(G4double x) {

  G4double y = -1;

  if(x<3.7||179.5<x) {
    G4cout << "Angle outside of range " << x << G4endl;
    return y;
  }
  
  G4int i=0;
  while(ThetaSc[i]<x) i++;

  G4double x1 = ThetaSc[i];	G4double y1 = SpinU[i];
  G4double x0 = ThetaSc[i-1];	G4double y0 = SpinU[i-1];
  
  y = y0 + (y1-y0)*(x-x0)/(x1-x0);

  return y;
}  


