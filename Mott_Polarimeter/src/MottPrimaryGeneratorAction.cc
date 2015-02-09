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

  // Resize data vectors
  ThetaSc.resize(5);
  CrossSection.resize(5);
  SpinT.resize(5);
  SpinU.resize(5);
  Sherman.resize(5);

  ReadDataFiles();

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
  // Random aspects must be input before calling GeneratePrimaryVertex();
  
  // Source Position
  G4int ThrowFromUpstream = 0;

  // Gausian beam profile.  
  G4double TargetLength = myDetector->GetTargetFullLength();
  G4double sigma = beamDiameter/(2.354820045*mm);

  G4double X = G4RandGauss::shoot(0.0,sigma)*mm;
  G4double Y = G4RandGauss::shoot(0.0,sigma)*mm;
  G4double Z = (G4UniformRand() - 0.5)*TargetLength + myDetector->GetTargetZPosition();

  G4ThreeVector gunPosition = G4ThreeVector(X,Y,Z);
  G4ThreeVector gunDirection = G4ThreeVector(0.0,0.0,1.0);
  
  // Incident electron's polarization set here
  G4double Sx = 0.0;
  G4double Sy = 1.0;
  G4double Sz = 0.0;
  G4ThreeVector preScatteredPolarization = G4ThreeVector(Sx,Sy,Sz);

  // Beam Energy 
  G4double energy = G4RandGauss::shoot(beamEnergy/(1.0*MeV), energySpread/(1.0*MeV))*MeV;

  G4double Phi = 0; 
  G4double Theta = 0;
  G4double MaxThrow = 8.0e-26;
  if(ThrowFromUpstream) {
    Z = -10.0*cm;
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
      G4double cs = InterpolateCrossSection(Theta/deg,energy/MeV);
      G4double s = InterpolateSherman(Theta/deg,energy/MeV);
      //G4double t = InterpolateT(Theta/deg,energy/MeV);
      //G4double u = InterpolateU(Theta/deg,energy/MeV);
      cs = cs*(1 + s*cos(Phi));
      G4double rejectionThrow = MaxThrow*G4UniformRand();
      if(rejectionThrow<=cs) {
      	goodThrow = 1;  
      	//G4cout << cs << "\t" << s << "\t" << t << "\t" << u << "\t" << G4endl;
      }
    }
      
    gunDirection.setRThetaPhi(1.0,Theta,Phi);   
    //G4cout << energy/MeV << "\t" << Theta/deg << "\t" << Phi/deg << "\t" << X/mm << "\t" << Y/mm << "\t" << Z/mm << G4endl;
  }
 
  // Set variable gun properties
  particleGun->SetParticleEnergy(energy);
  particleGun->SetParticlePosition(gunPosition);
  particleGun->SetParticleMomentumDirection(gunDirection); 
  particleGun->GeneratePrimaryVertex(anEvent);
  
  // std::cout << "\tLeaving MottPrimaryGeneratorAction::GeneratePrimaries()" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottPrimaryGeneratorAction::ReadDataFiles() {
  
  G4cout << "\tEntering MottPrimaryGeneratorAction::ReadDataFiles()" << G4endl;
  
  //Clear old data
  for(G4int i=0; i<5; i++) {
    ThetaSc[i].clear();
    CrossSection[i].clear();
    SpinT[i].clear();
    SpinU[i].clear();
    Sherman[i].clear();
  }
  
  //Read in Xavi's data.
  std::ifstream dcsFile;				// Cross section files
  std::ifstream spfFile;				// Spin transfer functions
  
  for(G4int i=0; i<5; i++) {  
  
    G4String skipLine;
    G4double theta, mu, dcs, dcs_angstrom;
    G4double s, t, u, unity;
    
    G4double fileEnergy = beamEnergy + 0.05*i-0.10;
    G4int energyInt = round(beamEnergy/MeV);
    G4int energyDecimal = 50*i-100;
    if(i<2) { 
      energyInt = energyInt-1;
      energyDecimal += 1000;
    }
    
    char dcsFileName[250], spfFileName[250]; 
    sprintf(dcsFileName,"/home/mjmchugh/Mott/MottG4/NewMottPhysics/Z=79_E=%4.2fMeV/dcs_%1dp%03de06.dat", fileEnergy, energyInt, energyDecimal);
    sprintf(spfFileName,"/home/mjmchugh/Mott/MottG4/NewMottPhysics/Z=79_E=%4.2fMeV/spf.dat", fileEnergy);

    dcsFile.open(dcsFileName);
    for(G4int nLines=1; nLines<27; nLines++) {
       std::getline(dcsFile, skipLine);
    }
  
    for(G4int nLines=27; nLines<=632; nLines++) {
      dcsFile >> theta >> mu >> dcs >> dcs_angstrom;
      ThetaSc[i].push_back(theta);
      CrossSection[i].push_back(dcs);
    }
    dcsFile.close();
  				
    spfFile.open(spfFileName);
    for(G4int nLines=1; nLines<7; nLines++) {
      std::getline(spfFile, skipLine);
    }

    for(G4int nLines=7; nLines<=612; nLines++) {
      spfFile >> theta >> mu >> u >> s >> t >> unity;
      SpinU[i].push_back(u);
      Sherman[i].push_back(s);
      SpinT[i].push_back(t);
    }
    spfFile.close();
  
  }  
  
  G4cout << "\tLeaving MottPrimaryGeneratorAction::ReadDataFiles()" << G4endl;  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MottPrimaryGeneratorAction::InterpolateCrossSection(G4double theta, G4double energy) {

  G4double F_xy = -1;

  if(theta<0.0||180.0<theta) {
    G4cout << "Angle outside of range " << theta << G4endl;
    return F_xy;
  }
  
  G4double E_lo=-1, E_hi=-1;
  G4double theta_lo=-1, theta_hi=-1;
  G4int j_lo=-1, j_hi=-1;
  G4int i_lo=-1, i_hi=-1;

  G4int j;
  G4double testEnergy;
  if (energy < beamEnergy - 0.1) {
    j_lo = 0;			j_hi = 1;
    E_lo = beamEnergy - 0.10;	E_hi = beamEnergy -0.05;
  } else if (energy > beamEnergy + 0.1) {
    j_lo = 3;			j_hi = 4;
    E_lo = beamEnergy + 0.05;	E_hi = beamEnergy + 0.10;
  }  else {
    for(j=1; j<5; j++) {
      testEnergy = beamEnergy + (0.05*j - 0.10);
      if(energy < testEnergy) {
        E_lo = testEnergy - 0.05;	E_hi = testEnergy;
        j_lo = j-1;		j_hi = j;
        break;
      }
    }
  }
  
  G4int i=0;
  while(ThetaSc[j_lo][i]<theta) i++;
  i_lo = i-1;			i_hi = i;
  theta_lo = ThetaSc[j_lo][i-1];	theta_hi = ThetaSc[j_lo][i];
   
  G4double F_11 = CrossSection[j_lo][i_lo]*std::abs(E_hi-energy)*std::abs(theta_hi-theta);
  G4double F_12 = CrossSection[j_lo][i_hi]*std::abs(E_hi-energy)*std::abs(theta_lo-theta);
  G4double F_21 = CrossSection[j_hi][i_lo]*std::abs(E_lo-energy)*std::abs(theta_hi-theta);
  G4double F_22 = CrossSection[j_hi][i_hi]*std::abs(E_lo-energy)*std::abs(theta_lo-theta);
  F_xy = (F_11 + F_12 + F_21 + F_22)/( (E_hi - E_lo)*(theta_hi-theta_lo) );
  
  return F_xy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double MottPrimaryGeneratorAction::InterpolateSherman(G4double theta, G4double energy) {

  G4double F_xy = -1;

  if(theta<0.0||180.0<theta) {
    G4cout << "Angle outside of range " << theta << G4endl;
    return F_xy;
  }
  
  G4double E_lo=-1, E_hi=-1;
  G4double theta_lo=-1, theta_hi=-1;
  G4int j_lo=-1, j_hi=-1;
  G4int i_lo=-1, i_hi=-1;

  G4int j;
  G4double testEnergy;
  if (energy < beamEnergy - 0.1) {
    j_lo = 0;			j_hi = 1;
    E_lo = beamEnergy - 0.10;	E_hi = beamEnergy -0.05;
  } else if (energy > beamEnergy + 0.1) {
    j_lo = 3;			j_hi = 4;
    E_lo = beamEnergy + 0.05;	E_hi = beamEnergy + 0.10;
  }  else {
    for(j=1; j<5; j++) {
      testEnergy = beamEnergy + (0.05*j - 0.10);
      if(energy < testEnergy) {
        E_lo = testEnergy - 0.05;	E_hi = testEnergy;
        j_lo = j-1;		j_hi = j;
        break;
      }
    }
  }
  
  G4int i=0;
  while(ThetaSc[j_lo][i]<theta) i++;
  i_lo = i-1;			i_hi = i;
  theta_lo = ThetaSc[j_lo][i-1];	theta_hi = ThetaSc[j_lo][i];
   
  G4double F_11 = Sherman[j_lo][i_lo]*std::abs(E_hi-energy)*std::abs(theta_hi-theta);
  G4double F_12 = Sherman[j_lo][i_hi]*std::abs(E_hi-energy)*std::abs(theta_lo-theta);
  G4double F_21 = Sherman[j_hi][i_lo]*std::abs(E_lo-energy)*std::abs(theta_hi-theta);
  G4double F_22 = Sherman[j_hi][i_hi]*std::abs(E_lo-energy)*std::abs(theta_lo-theta);
  F_xy = (F_11 + F_12 + F_21 + F_22)/( (E_hi - E_lo)*(theta_hi-theta_lo) );

  //G4cout << "\t\t" << Sherman[j_lo][i_lo] << "\t" <<  Sherman[j_lo][i_hi] << "\t" <<  Sherman[j_hi][i_lo] << "\t" <<  Sherman[j_hi][i_hi] << G4endl;
  //G4cout << "i and j: " << i_lo << " " << i_hi << " " << j_lo << " " << j_hi << G4endl;
  
  return F_xy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double MottPrimaryGeneratorAction::InterpolateT(G4double theta, G4double energy) {

  G4double F_xy = -1;

  if(theta<0.0||180.0<theta) {
    G4cout << "Angle outside of range " << theta << G4endl;
    return F_xy;
  }
  
  G4double E_lo=-1, E_hi=-1;
  G4double theta_lo=-1, theta_hi=-1;
  G4int j_lo=-1, j_hi=-1;
  G4int i_lo=-1, i_hi=-1;

  G4int j;
  G4double testEnergy;
  if (energy < beamEnergy - 0.1) {
    j_lo = 0;			j_hi = 1;
    E_lo = beamEnergy - 0.10;	E_hi = beamEnergy -0.05;
  } else if (energy > beamEnergy + 0.1) {
    j_lo = 3;			j_hi = 4;
    E_lo = beamEnergy + 0.05;	E_hi = beamEnergy + 0.10;
  }  else {
    for(j=1; j<5; j++) {
      testEnergy = beamEnergy + (0.05*j - 0.10);
      if(energy < testEnergy) {
        E_lo = testEnergy - 0.05;	E_hi = testEnergy;
        j_lo = j-1;			j_hi = j;
        break;
      }
    }
  }
  
  G4int i=0;
  while(ThetaSc[j_lo][i]<theta) i++;
  i_lo = i-1;				i_hi = i;
  theta_lo = ThetaSc[j_lo][i-1];	theta_hi = ThetaSc[j_lo][i];
   
  G4double F_11 = SpinT[j_lo][i_lo]*std::abs(E_hi-energy)*std::abs(theta_hi-theta);
  G4double F_12 = SpinT[j_lo][i_hi]*std::abs(E_hi-energy)*std::abs(theta_lo-theta);
  G4double F_21 = SpinT[j_hi][i_lo]*std::abs(E_lo-energy)*std::abs(theta_hi-theta);
  G4double F_22 = SpinT[j_hi][i_hi]*std::abs(E_lo-energy)*std::abs(theta_lo-theta);
  F_xy = (F_11 + F_12 + F_21 + F_22)/( (E_hi - E_lo)*(theta_hi-theta_lo) );
  
  return F_xy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double MottPrimaryGeneratorAction::InterpolateU(G4double theta, G4double energy) {

  G4double F_xy = -1;

  if(theta<0.0||180.0<theta) {
    G4cout << "Angle outside of range " << theta << G4endl;
    return F_xy;
  }
  
  G4double E_lo=-1, E_hi=-1;
  G4double theta_lo=-1, theta_hi=-1;
  G4int j_lo=-1, j_hi=-1;
  G4int i_lo=-1, i_hi=-1;

  G4int j;
  G4double testEnergy;
  if (energy < beamEnergy - 0.1) {
    j_lo = 0;			j_hi = 1;
    E_lo = beamEnergy - 0.10;	E_hi = beamEnergy -0.05;
  } else if (energy > beamEnergy + 0.1) {
    j_lo = 3;			j_hi = 4;
    E_lo = beamEnergy + 0.05;	E_hi = beamEnergy + 0.10;
  }  else {
    for(j=1; j<5; j++) {
      testEnergy = beamEnergy + (0.05*j - 0.10);
      if(energy < testEnergy) {
        E_lo = testEnergy - 0.05;	E_hi = testEnergy;
        j_lo = j-1;			j_hi = j;
        break;
      }
    }
  }
  
  G4int i=0;
  while(ThetaSc[j_lo][i]<theta) i++;
  i_lo = i-1;				i_hi = i;
  theta_lo = ThetaSc[j_lo][i-1];	theta_hi = ThetaSc[j_lo][i];
    
  G4double F_11 = SpinU[j_lo][i_lo]*std::abs(E_hi-energy)*std::abs(theta_hi-theta);
  G4double F_12 = SpinU[j_lo][i_hi]*std::abs(E_hi-energy)*std::abs(theta_lo-theta);
  G4double F_21 = SpinU[j_hi][i_lo]*std::abs(E_lo-energy)*std::abs(theta_hi-theta);
  G4double F_22 = SpinU[j_hi][i_hi]*std::abs(E_lo-energy)*std::abs(theta_lo-theta);
  F_xy = (F_11 + F_12 + F_21 + F_22)/( (E_hi - E_lo)*(theta_hi-theta_lo) );
  
  return F_xy;
}    

