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
// $Id: MottTrackerHit.hh,v 3.5 2013/12/02 21:53:35 mjmchugh Exp $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef MottTrackerHit_h
#define MottTrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MottTrackerHit : public G4VHit
{
  public:

      MottTrackerHit();
     ~MottTrackerHit();
      MottTrackerHit(const MottTrackerHit&);
      const MottTrackerHit& operator=(const MottTrackerHit&);
      G4int operator==(const MottTrackerHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:
  
      void SetCopyNo(G4int copy) { CopyNo = copy; };
      void SetTrackID(G4int track) { trackID = track; };
      void SetEdep(G4double de) { edep = de; };
      void SetTrackLength(G4double dl) { tracklength = dl; };
      void SetPos(G4ThreeVector xyz) { pos = xyz; };
      
      G4int GetCopyNo() { return CopyNo; };		
      G4int GetTrackID() { return trackID; };
      G4double GetEdep() { return edep; };
      G4double GetTrackLength() { return tracklength; };      
      G4ThreeVector GetPos() { return pos; };
      
  private:
  
      G4int CopyNo;
      G4int trackID;
      G4double edep;
      G4double tracklength;
      G4ThreeVector pos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<MottTrackerHit> MottTrackerHitsCollection;

extern G4Allocator<MottTrackerHit> MottTrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* MottTrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) MottTrackerHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void MottTrackerHit::operator delete(void *aHit)
{
  MottTrackerHitAllocator.FreeSingle((MottTrackerHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
