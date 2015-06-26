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
 
  EventType = 1;

  TargetZ = 79;

  ThetaMin = 0;
  ThetaMax = pi;
  PhiMin = 0;
  PhiMax = 2*pi;  

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
  Px1 = 0; 	Py1 = 0;	Pz1 = 0;
  X1 = 0;	Y1 = 0;		Z1 = 0;
  Theta1 = 0;	Phi1 = 0;	Energy1 = 0;
  CS1 = 0; 	S1 = 0; 	T1 = 0; 	U1 = 0;
  Px2 = 0; 	Py2 = 0;	Pz2 = 0;
  X2 = 0;	Y2 = 0;		Z2 = 0;
  Theta2 = 0;	Phi2 = 0;	Energy2 = 0;
  CS2 = 0; 	S2 = 0; 	T2 = 0; 	U2 = 0;
  
  // Gausian beam profile.  
  G4double TargetLength = myDetector->GetTargetFullLength();
  G4double sigma = beamDiameter/(2.354820045*mm);
  X1 = G4RandGauss::shoot(0.0,sigma)*mm;
  Y1 = G4RandGauss::shoot(0.0,sigma)*mm;
  G4double depth = G4UniformRand()*TargetLength;
  Z1 = depth - 0.5*TargetLength + myDetector->GetTargetZPosition();
  G4ThreeVector gunPosition = G4ThreeVector(X1,Y1,Z1);
  G4ThreeVector gunDirection = G4ThreeVector(0.0,0.0,1.0);
  
  // Incident electron's polarization set here
  G4double Py1 = 1.0;
  G4ThreeVector P_1 = G4ThreeVector(Px1,Py1,Pz1);

  // Electron Energy at the scattering vertex (Gaussian minus ~linear ELoss up to the scattering vertex)
  Energy1 = G4RandGauss::shoot(beamEnergy/(1.0*MeV), energySpread/(1.0*MeV))*MeV;
  Energy1 -= CalculateTotalELoss(depth, Energy1, TargetZ);
  G4double Energy = 0;

  if(EventType == 0) {  							// Throw from upstream

    Z1 = -10.0*cm;
    gunPosition = G4ThreeVector(X1, Y1, Z1);

  } else if(EventType == 1) {							// Throw single scattered electrons

    G4int goodThrow = 0;
    G4double MaxThrow = 1.0e-25;
    while (goodThrow==0) {
      G4double x_0 = 0.0;
      G4double y_0 = 0.0*mm;
      G4double z_0 = -277.597*mm; 						// Front face of collimator
      G4double CoinToss = G4UniformRand();
      if(0.0<=CoinToss&&CoinToss<0.5) {						// Pick L or R detector
        x_0 = - 2.850*25.4*mm/2.0;
      } else if(0.5<=CoinToss&&CoinToss<=1.0) {
        x_0 = 2.850*25.4*mm/2.0;
      } 
      G4double R = 0.192*25.4*mm/2.0;						// radius of collimator opening
      G4double ph = 2*pi*G4UniformRand();
      G4double u = R*(G4UniformRand() + G4UniformRand());
      if (u > R) u = 2*R-u;
      G4ThreeVector throwTo = G4ThreeVector(u*cos(ph)+x_0, u*sin(ph)+y_0, z_0); // location in collimator opening
      G4ThreeVector d_1 = throwTo - gunPosition;
      d_1 = d_1.unit();
      G4ThreeVector n_1 = gunDirection.cross(d_1);
      n_1 = n_1.unit();
      Theta1 = d_1.theta();
      Phi1 = d_1.phi();
      CS1 = InterpolateCrossSection(Theta1/deg,Energy1/MeV);
      S1 = InterpolateSherman(Theta1/deg,Energy1/MeV);
      T1 = InterpolateT(Theta1/deg,Energy1/MeV);
      U1 = InterpolateU(Theta1/deg,Energy1/MeV);
      CS1 = CS1*(1 + S1*n_1.dot(P_1));
      G4double rejectionThrow = MaxThrow*G4UniformRand();
      if(rejectionThrow<=CS1) {
      	goodThrow = 1;
        G4ThreeVector P_2 = CalculateNewPol(n_1,P_1, S1, T1, U1);
        Px2 = P_2.x();
        Py2 = P_2.y();
        Pz2 = P_2.z();
        Energy = Energy1;
      }
    }      
    gunDirection.setRThetaPhi(1.0, Theta1, Phi1);

  } else if (EventType == 2) {						// Throw double scattered electrons		
									// at the collimator
    G4int goodThrow = 0;
    G4double MaxThrow = 2.0e-46;
    while(goodThrow == 0) {
      // Kinematic Considerations
      // Pick point for second scattering (within a disk within a given radius of the first scattering)
      G4double R_2 = G4UniformRand()*0.157*mm;
      G4double ph_2 = 2*pi*G4UniformRand();
      X2 = R_2*cos(ph_2)+X1;
      Y2 = R_2*sin(ph_2)+Y1;
      Z2 = (G4UniformRand() - 0.5)*TargetLength + myDetector->GetTargetZPosition();
      G4ThreeVector x_2 = G4ThreeVector(X2, Y2, Z2);
      // pick point in collimator acceptance to throw to
      G4double x_0 = 0.0;
      G4double CoinToss = G4UniformRand();
      if(0.0<=CoinToss && CoinToss<0.5) {					// Pick L or R detector
        x_0 = - 2.850*25.4*mm/2.0;
      } else if(0.5<=CoinToss && CoinToss<=1.0) {
        x_0 = 2.850*25.4*mm/2.0;
      } 
      G4double y_0 = 0.0*mm;
      G4double z_0 = -277.597*mm; 						// Front face of collimator
      G4double R = 0.192*25.4*mm/2.0;						// radius of collimator opening
      G4double ph = 2*pi*G4UniformRand();
      G4double u = G4UniformRand() + G4UniformRand();
      if (u > 1.0) {
        u = R*(2-u);
      } else {
        u = R*u;
      }
      G4ThreeVector x_3 = G4ThreeVector(u*cos(ph)+x_0, u*sin(ph)+y_0, z_0); // location in collimator opening    
      // calculate things for first scattering.
      G4ThreeVector d_1 = x_2 - gunPosition;
      G4double d_1_length = d_1.mag()/mm;
      G4ThreeVector n_1 = gunDirection.cross(d_1);
      n_1 = n_1.unit();
      Theta1 = d_1.theta();
      Phi1 = d_1.phi();
      G4double CS1 = InterpolateCrossSection(Theta1/deg,Energy1/MeV);
      G4double S1 = InterpolateSherman(Theta1/deg,Energy1/MeV);
      G4double T1 = InterpolateT(Theta1/deg,Energy1/MeV);
      G4double U1 = InterpolateU(Theta1/deg,Energy1/MeV);
      CS1 = CS1*(1 + S1*n_1.dot(P_1));    
      G4ThreeVector P_2 = CalculateNewPol(n_1,P_1,S1,T1,U1);
      Px2 = P_2.x();
      Py2 = P_2.y();
      Pz2 = P_2.z();
      // Calculate things for second scattering
      G4ThreeVector d_2 = x_3-x_2;
      G4double Theta2 = d_1.angle(d_2);
      G4ThreeVector n_2 = d_1.cross(d_2);
      n_2 = n_2.unit();
      G4double Energy2 = Energy1 - CalculateTotalELoss(d_1_length, Energy1, TargetZ);
      G4double CS2 = InterpolateCrossSection(Theta2/deg,Energy2/MeV);
      G4double S2 = InterpolateSherman(Theta2/deg,Energy2/MeV);
      G4double T2 = InterpolateT(Theta2/deg,Energy2/MeV);
      G4double U2 = InterpolateU(Theta2/deg,Energy2/MeV);
      CS2 = CS2*(1 + S2*n_2.dot(P_2));
      G4double CS = CS2*CS1;
      G4double rejectionThrow = MaxThrow*G4UniformRand();
      if(rejectionThrow<=CS) {
        goodThrow = 1;
        gunDirection = d_2.unit();
        Energy = Energy2;
      }
    }

  } else {					// Throw uniformly across the user specified angular range 

    Theta1 = acos( (cos(ThetaMin)-cos(ThetaMax))*G4UniformRand() - cos(ThetaMax) );
    Phi1 = PhiMin + (PhiMax-PhiMin)*G4UniformRand();
    CS1 = InterpolateCrossSection(Theta1/deg,Energy1/MeV);
    S1 = InterpolateSherman(Theta1/deg,Energy1/MeV);
    T1 = InterpolateT(Theta1/deg,Energy1/MeV);
    U1 = InterpolateU(Theta1/deg,Energy1/MeV);
    CS1 = CS1*(1 + S1*cos(Phi1));

  }

  // Primary vertex quantitites to store in rootfile
  pEventAction->SetKEPrime(Energy1, 0);			// Energy
  pEventAction->SetXPos(X1, 0);				// location
  pEventAction->SetYPos(Y1, 0);
  pEventAction->SetZPos(Z1, 0);
  pEventAction->SetTheta(Theta1, 0);			// Scattering angle
  pEventAction->SetPhi(Phi1, 0);  			// Azimuthal angle
  pEventAction->SetXPol(Px1, 0);			// Incoming Polarization
  pEventAction->SetYPol(Py1, 0);
  pEventAction->SetZPol(Pz1, 0);
  pEventAction->SetCS(CS1, 0);				// Differential Cross Section
  pEventAction->SetS(S1, 0);				// Sherman Function 
  pEventAction->SetT(T1, 0);				// SpinT
  pEventAction->SetU(U1, 0);				// SpinU
  // Secondary scattering info
  pEventAction->SetKEPrime(Energy2, 1);			// Energy
  pEventAction->SetXPos(X2, 1);				// location
  pEventAction->SetYPos(Y2, 1);
  pEventAction->SetZPos(Z2, 1);
  pEventAction->SetTheta(Theta2, 1);			// Scattering angle
  pEventAction->SetPhi(Phi2, 1);  			// Azimuthal angle
  pEventAction->SetXPol(Px2, 1);			// Incoming Polarization
  pEventAction->SetYPol(Py2, 1);
  pEventAction->SetZPol(Pz2, 1);
  pEventAction->SetCS(CS2, 1);				// Differential Cross Section
  pEventAction->SetS(S2, 1);				// Sherman Function 
  pEventAction->SetT(T2, 1);				// SpinT
  pEventAction->SetU(U2, 1);				// SpinU

  // Set variable gun properties
  particleGun->SetParticleEnergy(Energy);			// scattered electron KE
  particleGun->SetParticlePosition(gunPosition);
  particleGun->SetParticleMomentumDirection(gunDirection); 
  particleGun->GeneratePrimaryVertex(anEvent);

  //G4cout << "\t " << anEvent->GetEventID() << G4endl;
  
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
G4ThreeVector MottPrimaryGeneratorAction::CalculateNewPol(G4ThreeVector n, G4ThreeVector P, G4double s, G4double t, G4double u) {
  
  // Define the basis vectors used for the new polarization
  // normal perpendicular to the scattering plane n
  G4ThreeVector e1 = n;
  G4ThreeVector e2 = n.cross(P);	// n x P
  G4ThreeVector e3 = n.cross(-1.0*e2);	// n x ( P x n )

  // Normalization of the new polarization vector
  G4double norm = 1 + s*n.dot(P);
  
  G4ThreeVector P_new = G4ThreeVector(0.0, 0.0, 0.0);
  P_new = ((n.dot(P) + s)*e1 + u*e2 + t*e3)/norm;
  
  return P_new;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double MottPrimaryGeneratorAction::CalculateInstantaneousELoss(G4double E, G4int Z = 79) {
  
  G4double ELoss = 0.0;
  if (Z==79) {
    ELoss = 1.888*MeV/mm + 0.2723*E/mm;
  } else if (Z==47) {
    ELoss = 1.165*MeV/mm + 0.1142*E/mm; 
  }
  return ELoss;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double MottPrimaryGeneratorAction::CalculateTotalELoss(G4double x, G4double E_0, G4int Z) {

  G4double E = E_0;
  G4double stepLength = x/10.0;
  //G4cout << stepLength/mm << " mm";
  for (G4int i=0; i<10; i++) {
    E -= stepLength*CalculateInstantaneousELoss(E,Z);
    //G4cout << "\t" << stepLength*CalculateInstantaneousELoss(E,Z);
  }
  G4double DeltaE = E_0 - E; 
  //G4cout << "\t" << DeltaE << G4endl;
  return DeltaE;
}
