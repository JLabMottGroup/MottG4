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
// $Id: MottEventAction.hh,v 3.6 2013/12/02 21:53:35 mjmchugh Exp $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#ifndef MottEventAction_h
#define MottEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
class G4Event;
class MottEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MottEventAction : public G4UserEventAction
{
  public:
    MottEventAction();
   ~MottEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    
    inline void SetEdep(G4double edep, G4int i) { Edep[i] = edep; }; 
    inline void SetdEdep(G4double edep, G4int i) { dEdep[i] = edep; };
    inline void SetEtrackL(G4double dl, G4int i) { EtrackL[i] = dl; }; 
    inline void SetdEtrackL(G4double dl, G4int i) { dEtrackL[i] = dl; };
    inline void SetNumEPhotons(G4int nhits, G4int i) { NumEPhotons[i] = nhits; };
    inline void SetNumdEPhotons(G4int nhits, G4int i) { NumdEPhotons[i] = nhits; };
    
    /* Dump Plate
    inline void SetBeEnergyDeposited(G4double edep, G4int i) { BeEnergyDeposited[i] = edep; };
    inline void AddBeEnergyDeposited(G4double edep, G4int i) { BeEnergyDeposited[i] += edep; };
    G4double GetBeEnergyDeposited(G4int i) { return BeEnergyDeposited[i]; };
    inline void SetCuEnergyDeposited(G4double edep, G4int i) { CuEnergyDeposited[i] = edep; };
    inline void AddCuEnergyDeposited(G4double edep, G4int i) { CuEnergyDeposited[i] += edep; };
    G4double GetCuEnergyDeposited(G4int i) { return CuEnergyDeposited[i]; };
    */
    
    G4double GetEdep(G4int i) { return Edep[i]; };
    G4double GetdEdep(G4int i) { return dEdep[i]; };
    G4double GetEtrackL(G4int i) { return EtrackL[i]; };
    G4double GetdEtrackL(G4int i) { return dEtrackL[i]; };
    G4int GetNumEPhotons(G4int i) { return NumEPhotons[i]; };    
    G4int GetNumdEPhotons(G4int i) { return NumdEPhotons[i]; };  
    G4int GetNEPE(G4int i) { return NEPE[i]; };  
    G4int GetndEPE(G4int i) { return NdEPE[i]; };

    // Vertex Set()/Get()
    inline void SetKEPrime(G4double KE, G4int i) { KEPrime[i] = KE; };
    G4double GetKEPrime(G4int i) { return KEPrime[i]; };
    inline void SetXPos(G4double X, G4int i) { XPos[i] = X; };
    G4double GetXPos(G4int i) { return XPos[i]; };
    inline void SetYPos(G4double Y, G4int i) { YPos[i] = Y; };
    G4double GetYPos(G4int i) { return YPos[i]; };
    inline void SetZPos(G4double Z, G4int i) { ZPos[i] = Z; };
    G4double GetZPos(G4int i) { return ZPos[i]; };
    inline void SetTheta(G4double theta, G4int i) { Theta[i] = theta; };
    G4double GetTheta(G4int i) { return Theta[i]; };
    inline void SetPhi(G4double phi, G4int i) { Phi[i] = phi; };
    G4double GetPhi(G4int i) { return Phi[i]; };
    inline void SetXPol(G4double X, G4int i) { XPol[i] = X; };
    G4double GetXPol(G4int i) { return XPol[i]; };
    inline void SetYPol(G4double Y, G4int i) { YPol[i] = Y; };
    G4double GetYPol(G4int i) { return YPol[i]; };
    inline void SetZPol(G4double Z, G4int i) { ZPol[i] = Z; };
    G4double GetZPol(G4int i) { return ZPol[i]; };
    inline void SetCS(G4double cs, G4int i) { CS[i] = cs; };
    G4double GetCS(G4int i) { return CS[i]; };
    inline void SetS(G4double s, G4int i) { S[i] = s; };
    G4double GetS(G4int i) { return S[i]; };
    inline void SetT(G4double t, G4int i) { T[i] = t; };
    G4double GetT(G4int i) { return T[i]; };
    inline void SetU(G4double u, G4int i) { U[i] = u; };
    G4double GetU(G4int i) { return U[i]; };
    inline void SetScatTheta(G4double theta) { ScatTheta = theta; };
    G4double GetScatTheta() { return ScatTheta; };
    inline void SetScatPhi(G4double phi) { ScatPhi = phi; };
    G4double GetPhi() { return ScatPhi; };

    void SetStoreAll(G4int value = 1) { StoreAll = value; };

  private:

    MottEventActionMessenger* myMessenger;

    // Store all events flag
    G4int StoreAll;
  
    // Main Detectors
    G4double Edep[4];
    G4double dEdep[4];
    G4double EtrackL[4];
    G4double dEtrackL[4];
    G4int NumEPhotons[4]; 
    G4int NumdEPhotons[4];
    G4int NEPE[4];
    G4int NdEPE[4];     

    // Vertex Info
    G4double KEPrime[2];
    G4double XPos[2];	// Scattering Vertex location
    G4double YPos[2];
    G4double ZPos[2];
    G4double Theta[2];
    G4double Phi[2]; 
    G4double XPol[2];	// New Polarization vector
    G4double YPol[2];
    G4double ZPol[2];

    // Dynamic Variables
    G4double CS[2];
    G4double S[2];
    G4double T[2];
    G4double U[2];

    G4double ScatTheta;
    G4double ScatPhi;

    // Dump "Detectors"
    //G4double BeEnergyDeposited[6*20];
    //G4double CuEnergyDeposited[18*20];

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
