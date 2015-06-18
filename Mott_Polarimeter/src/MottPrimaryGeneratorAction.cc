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
#include "MottEventAction.hh"

#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <iostream>
#include <cmath>
#include <map>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottPrimaryGeneratorAction::MottPrimaryGeneratorAction(
                                               MottDetectorConstruction* myDC)
:myDetector(myDC)
{
  //G4cout << "\tEntering MottPrimaryGeneratorAction::MottPrimaryGeneratorAction()" <<G4endl; 

  G4int n_particle = 1;

  myMessenger = new MottPrimaryGeneratorMessenger(this);
  particleGun = new G4ParticleGun(n_particle);
  pEventAction = (MottEventAction*) G4RunManager::GetRunManager()->GetUserEventAction();
 
  ThrowFromUpstream = false;
  ThrowAtCollimators = true;

  TargetZ = 79;

  GoldELoss[3.0*MeV] = 27.02*MeV/cm;
  GoldELoss[5.0*MeV] = 32.54*MeV/cm;
  GoldELoss[6.0*MeV] = 35.26*MeV/cm;
  GoldELoss[8.0*MeV] = 40.65*MeV/cm;

  // Set default 
  beamEnergy = 5.0*MeV; 	// kinetic energy
  energySpread = 25.0*keV;
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

  //G4cout << "\tLeaving MottPrimaryGeneratorAction::MottPrimaryGeneratorAction()" <<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottPrimaryGeneratorAction::~MottPrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// All Gun settings changed here will change every event. 
// Random aspects must be input before calling GeneratePrimaryVertex();
void MottPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //std::cout << "\tEntering MottPrimaryGeneratorAction::GeneratePrimaries()" << std::endl;

  pEventAction = (MottEventAction*) G4RunManager::GetRunManager()->GetUserEventAction();
  
  // Initialize all variables to 0
  X = 0;	Y = 0;		Z = 0;
  Px1 = 0; 	Py1 = 0;	Pz1 = 0;
  Px2 = 0; 	Py2 = 0;	Pz2 = 0;
  Theta = 0;	Phi = 0;	Energy = 0;
  CS = 0; 	S = 0; 		T = 0; 		U = 0;

  // Gausian beam profile.  
  G4double TargetLength = myDetector->GetTargetFullLength();
  G4double sigma = beamDiameter/(2.354820045*mm);
  X = G4RandGauss::shoot(0.0,sigma)*mm;
  Y = G4RandGauss::shoot(0.0,sigma)*mm;
  G4double depth = G4UniformRand()*TargetLength;
  Z = depth - 0.5*TargetLength + myDetector->GetTargetZPosition();
  
  G4ThreeVector gunPosition = G4ThreeVector(X,Y,Z);
  G4ThreeVector gunDirection = G4ThreeVector(0.0,0.0,1.0);
  
  // Incident electron's polarization set here
  G4double Py1 = 1.0;
  G4ThreeVector preScatteredPolarization = G4ThreeVector(Px1,Py1,Pz1);

  // Electron Energy at the scattering vertex (Gaussian minus ~linear ELoss up to the scattering vertex)
  G4double Energy = G4RandGauss::shoot(beamEnergy/(1.0*MeV), energySpread/(1.0*MeV))*MeV 
                    - GoldELoss[beamEnergy]*depth;

  if(ThrowFromUpstream) {
    Z = -10.0*cm;
    gunPosition = G4ThreeVector(X, Y, Z);
  } else if(ThrowAtCollimators) {
    G4int goodThrow = 0;
    G4double MaxThrow = 1.6*CrossSection[0][589];				// 1.6 times the lowest energy dcs at 172 deg.
    while (goodThrow==0) { 
      G4double ScatteringAngle = 172.6*deg;					// Average acceptance angle.
      Theta = ScatteringAngle - 0.6*deg + 1.2*G4UniformRand()*deg;	        // Throw from 172.0 to 173.2 degrees
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
      CS = InterpolateCrossSection(Theta/deg,Energy/MeV);
      S = InterpolateSherman(Theta/deg,Energy/MeV);
      T = InterpolateT(Theta/deg,Energy/MeV);
      U = InterpolateU(Theta/deg,Energy/MeV);
      CS = CS*(1 + S*cos(Phi));
      G4double rejectionThrow = MaxThrow*G4UniformRand();
      if(rejectionThrow<=CS) {
      	goodThrow = 1;  
      	//G4cout << cs << "\t" << s << "\t" << t << "\t" << u << "\t" << G4endl;
      }
    }      
    gunDirection.setRThetaPhi(1.0,Theta,Phi);
  } else {
    Theta = ThetaMin + (ThetaMax-ThetaMin)*G4UniformRand();
    Phi = PhiMin + (PhiMax-PhiMin)*G4UniformRand();
    CS = InterpolateCrossSection(Theta/deg,Energy/MeV);
    S = InterpolateSherman(Theta/deg,Energy/MeV);
    T = InterpolateT(Theta/deg,Energy/MeV);
    U = InterpolateU(Theta/deg,Energy/MeV);
    CS = CS*(1 + S*cos(Phi));
  }

  // Calculate the outgoing electron's polarization (post scattering)
  CalculateNewPol();

  //G4cout << ThrowFromUpstream << " " << ThrowAtCollimators << G4endl;
  G4cout << Energy << " " << X/mm << " " << Y/mm << " " << Z/mm << " " << Theta/deg << " " << Phi/deg << " "
         << Px2 << " " << Py2 << " " << Pz2 << " " << CS << " " << S << " " << T << " " << U << G4endl;

  // Primary verted quantitites to store in rootfile
  pEventAction->SetKEPrime(Energy);			// Energy
  pEventAction->SetXPos(X);				// location
  pEventAction->SetYPos(Y);
  pEventAction->SetZPos(Z);
  pEventAction->SetTheta(Theta);			// Scattering angle
  pEventAction->SetPhi(Phi);  				// Azimuthal angle
  pEventAction->SetXPol(Px2);				// Outgoing Polarization
  pEventAction->SetYPol(Py2);
  pEventAction->SetZPol(Pz2);
  pEventAction->SetCS(CS);				// Differential Cross Section
  pEventAction->SetS(S);				// Sherman Function 
  pEventAction->SetT(T);				// SpinT
  pEventAction->SetU(U);				// SpinU

  //G4cout << energy << "\t" << X << "\t" << Y << "\t" << Z << "\t" << Theta << "\t" << Phi << "\t" << G4endl;

  // Set variable gun properties
  particleGun->SetParticleEnergy(Energy);			// scattered electron KE
  particleGun->SetParticlePosition(gunPosition);
  particleGun->SetParticleMomentumDirection(gunDirection); 
  particleGun->GeneratePrimaryVertex(anEvent);
  
  //std::cout << "\tLeaving MottPrimaryGeneratorAction::GeneratePrimaries()" << std::endl;
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
  
  char* MOTTG4DIR = getenv("MOTTG4DIR");

  //Read in Xavi's data.
  std::ifstream dcsFile;				// Cross section files
  std::ifstream spfFile;				// Spin transfer functions

  if(myDetector->GetTargetMater()!=NULL) {
    G4cout << "Target info "  << myDetector->GetTargetMater()->GetName() << " " << myDetector->GetTargetMater()->GetZ() << G4endl;
    TargetZ = myDetector->GetTargetMater()->GetZ();
  }  

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
    sprintf(dcsFileName,"%s/NewMottPhysics/Z=%2.0f_E=%4.2fMeV/dcs_%1dp%03de06.dat", MOTTG4DIR, TargetZ, fileEnergy, energyInt, energyDecimal);
    sprintf(spfFileName,"%s/NewMottPhysics/Z=%2.0f_E=%4.2fMeV/spf.dat", MOTTG4DIR, TargetZ, fileEnergy);

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
    E_lo = beamEnergy - 0.10;	E_hi = beamEnergy - 0.05;
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
  i_lo = i-1;			i_hi = i;
  theta_lo = ThetaSc[j_lo][i-1];	theta_hi = ThetaSc[j_lo][i];

  if( beamEnergy-0.10<=energy && energy<=beamEnergy+0.10) {   
    G4double F_11 = CrossSection[j_lo][i_lo]*std::abs(E_hi-energy)*std::abs(theta_hi-theta);
    G4double F_12 = CrossSection[j_lo][i_hi]*std::abs(E_hi-energy)*std::abs(theta_lo-theta);
    G4double F_21 = CrossSection[j_hi][i_lo]*std::abs(E_lo-energy)*std::abs(theta_hi-theta);
    G4double F_22 = CrossSection[j_hi][i_hi]*std::abs(E_lo-energy)*std::abs(theta_lo-theta);
    F_xy = (F_11 + F_12 + F_21 + F_22)/( (E_hi - E_lo)*(theta_hi-theta_lo) );  
  } else if ( energy < beamEnergy-0.10 ) {
    G4double G_1 = (CrossSection[j_hi][i_lo]-CrossSection[j_lo][i_lo])/(E_hi-E_lo)*(energy - E_lo) +  CrossSection[j_lo][i_lo];
    G4double G_2 = (CrossSection[j_hi][i_hi]-CrossSection[j_lo][i_hi])/(E_hi-E_lo)*(energy - E_lo) +  CrossSection[j_lo][i_hi];
    F_xy = (theta_hi-theta)/(theta_hi-theta_lo)*G_1 + (theta-theta_lo)/(theta_hi-theta_lo)*G_2;    
  } else if ( energy > beamEnergy+0.10 ) {
    G4double G_1 = (CrossSection[j_hi][i_lo]-CrossSection[j_lo][i_lo])/(E_hi-E_lo)*(energy - E_hi) +  CrossSection[j_hi][i_lo];
    G4double G_2 = (CrossSection[j_hi][i_hi]-CrossSection[j_lo][i_hi])/(E_hi-E_lo)*(energy - E_hi) +  CrossSection[j_hi][i_hi];
    F_xy = (theta_hi-theta)/(theta_hi-theta_lo)*G_1 + (theta-theta_lo)/(theta_hi-theta_lo)*G_2;
  }

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
   
  if( beamEnergy-0.10<=energy && energy<=beamEnergy+0.10) {   
    G4double F_11 = Sherman[j_lo][i_lo]*std::abs(E_hi-energy)*std::abs(theta_hi-theta);
    G4double F_12 = Sherman[j_lo][i_hi]*std::abs(E_hi-energy)*std::abs(theta_lo-theta);
    G4double F_21 = Sherman[j_hi][i_lo]*std::abs(E_lo-energy)*std::abs(theta_hi-theta);
    G4double F_22 = Sherman[j_hi][i_hi]*std::abs(E_lo-energy)*std::abs(theta_lo-theta);
    F_xy = (F_11 + F_12 + F_21 + F_22)/( (E_hi - E_lo)*(theta_hi-theta_lo) );  
  } else if ( energy < beamEnergy-0.10 ) {
    G4double G_1 = (Sherman[j_hi][i_lo]-Sherman[j_lo][i_lo])/(E_hi-E_lo)*(energy - E_lo) +  Sherman[j_lo][i_lo];
    G4double G_2 = (Sherman[j_hi][i_hi]-Sherman[j_lo][i_hi])/(E_hi-E_lo)*(energy - E_lo) +  Sherman[j_lo][i_hi];
    F_xy = (theta_hi-theta)/(theta_hi-theta_lo)*G_1 + (theta-theta_lo)/(theta_hi-theta_lo)*G_2;    
  } else if ( energy > beamEnergy+0.10 ) {
    G4double G_1 = (Sherman[j_hi][i_lo]-Sherman[j_lo][i_lo])/(E_hi-E_lo)*(energy - E_hi) +  Sherman[j_hi][i_lo];
    G4double G_2 = (Sherman[j_hi][i_hi]-Sherman[j_lo][i_hi])/(E_hi-E_lo)*(energy - E_hi) +  Sherman[j_hi][i_hi];
    F_xy = (theta_hi-theta)/(theta_hi-theta_lo)*G_1 + (theta-theta_lo)/(theta_hi-theta_lo)*G_2;
  }

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
   
  if( beamEnergy-0.10<=energy && energy<=beamEnergy+0.10) {   
    G4double F_11 = SpinT[j_lo][i_lo]*std::abs(E_hi-energy)*std::abs(theta_hi-theta);
    G4double F_12 = SpinT[j_lo][i_hi]*std::abs(E_hi-energy)*std::abs(theta_lo-theta);
    G4double F_21 = SpinT[j_hi][i_lo]*std::abs(E_lo-energy)*std::abs(theta_hi-theta);
    G4double F_22 = SpinT[j_hi][i_hi]*std::abs(E_lo-energy)*std::abs(theta_lo-theta);
    F_xy = (F_11 + F_12 + F_21 + F_22)/( (E_hi - E_lo)*(theta_hi-theta_lo) );  
  } else if ( energy < beamEnergy-0.10 ) {
    G4double G_1 = (SpinT[j_hi][i_lo]-SpinT[j_lo][i_lo])/(E_hi-E_lo)*(energy - E_lo) +  SpinT[j_lo][i_lo];
    G4double G_2 = (SpinT[j_hi][i_hi]-SpinT[j_lo][i_hi])/(E_hi-E_lo)*(energy - E_lo) +  SpinT[j_lo][i_hi];
    F_xy = (theta_hi-theta)/(theta_hi-theta_lo)*G_1 + (theta-theta_lo)/(theta_hi-theta_lo)*G_2;    
  } else if ( energy > beamEnergy+0.10 ) {
    G4double G_1 = (SpinT[j_hi][i_lo]-SpinT[j_lo][i_lo])/(E_hi-E_lo)*(energy - E_hi) +  SpinT[j_hi][i_lo];
    G4double G_2 = (SpinT[j_hi][i_hi]-SpinT[j_lo][i_hi])/(E_hi-E_lo)*(energy - E_hi) +  SpinT[j_hi][i_hi];
    F_xy = (theta_hi-theta)/(theta_hi-theta_lo)*G_1 + (theta-theta_lo)/(theta_hi-theta_lo)*G_2;
  }
  
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
    
  if( beamEnergy-0.10<=energy && energy<=beamEnergy+0.10) {   
    G4double F_11 = SpinU[j_lo][i_lo]*std::abs(E_hi-energy)*std::abs(theta_hi-theta);
    G4double F_12 = SpinU[j_lo][i_hi]*std::abs(E_hi-energy)*std::abs(theta_lo-theta);
    G4double F_21 = SpinU[j_hi][i_lo]*std::abs(E_lo-energy)*std::abs(theta_hi-theta);
    G4double F_22 = SpinU[j_hi][i_hi]*std::abs(E_lo-energy)*std::abs(theta_lo-theta);
    F_xy = (F_11 + F_12 + F_21 + F_22)/( (E_hi - E_lo)*(theta_hi-theta_lo) );  
  } else if ( energy < beamEnergy-0.10 ) {
    G4double G_1 = (SpinU[j_hi][i_lo]-SpinU[j_lo][i_lo])/(E_hi-E_lo)*(energy - E_lo) +  SpinU[j_lo][i_lo];
    G4double G_2 = (SpinU[j_hi][i_hi]-SpinU[j_lo][i_hi])/(E_hi-E_lo)*(energy - E_lo) +  SpinU[j_lo][i_hi];
    F_xy = (theta_hi-theta)/(theta_hi-theta_lo)*G_1 + (theta-theta_lo)/(theta_hi-theta_lo)*G_2;    
  } else if ( energy > beamEnergy+0.10 ) {
    G4double G_1 = (SpinU[j_hi][i_lo]-SpinU[j_lo][i_lo])/(E_hi-E_lo)*(energy - E_hi) +  SpinU[j_hi][i_lo];
    G4double G_2 = (SpinU[j_hi][i_hi]-SpinU[j_lo][i_hi])/(E_hi-E_lo)*(energy - E_hi) +  SpinU[j_hi][i_hi];
    F_xy = (theta_hi-theta)/(theta_hi-theta_lo)*G_1 + (theta-theta_lo)/(theta_hi-theta_lo)*G_2;
  }
  
  return F_xy;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MottPrimaryGeneratorAction::SetBeamEnergy(G4double energy) {

  G4double E = energy/MeV;

  if ( !(E==3.0||E==5.0||E==6.0||E==8.0) ) {
    G4cout << "Please enter one of the four acceptable kinetic energies." << G4endl;
    return;
  }

  beamEnergy = energy;

  ReadDataFiles();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int MottPrimaryGeneratorAction::CalculateNewPol() {
  
  // Define the basis vectors used for the new polarization
  // normal perpendicular to the scattering plane n
  G4double e1x = -sin(Phi);
  G4double e1y = cos(Phi);
  G4double e1z = 0;
  
  // n x P 
  G4double e2x = 0;
  G4double e2y = 0;
  G4double e2z = -sin(Phi);

  // n x ( P x n )
  G4double e3x = sin(Phi)*cos(Phi);
  G4double e3y = sin(Phi)*sin(Phi);
  G4double e3z = 0;

  // Normalization of the new polarization vector
  G4double norm = 1 + cos(Phi)*S;
  
  if(norm > 0) {
    // New Polarization
    Px2 = ((cos(Phi) + S)*e1x + U*e2x + T*e3x)/norm;
    Py2 = ((cos(Phi) + S)*e1y + U*e2y + T*e3y)/norm;
    Pz2 = ((cos(Phi) + S)*e1z + U*e2z + T*e3z)/norm;
  } else {
    Px2 = 0;
    Py2 = 0;
    Pz2 = 0;
  }

  return 0;
}
