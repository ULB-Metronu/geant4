//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ProductionCutsTable.hh,v 1.2 2003/11/03 02:18:44 kurasige Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//
// ------------------------------------------------------------
//   First Implementation          05 Oct. 2002  M.Asai    
// ------------------------------------------------------------

#ifndef G4ProductionCutsTable_h 
#define G4ProductionCutsTable_h 1

class G4RegionStore;
class G4VRangeToEnergyConverter;
class G4LogicalVolume;
class G4ProductionCuts;

#include "globals.hh"
#include "G4ios.hh"
#include <vector>
#include "G4MaterialCutsCouple.hh"
#include "G4Region.hh"


class G4ProductionCutsTable  
{
  // Class Description
  //  G4ProductionCutsTable is a static singleton class of a table of
  //  G4ProductionCuts objects. This class also manages tables of
  //  production cut and energy cut for each particle type.

  public: // with description
    static G4ProductionCutsTable* GetProductionCutsTable();
    // This static method returns the singleton pointer of this class object.
    // At the first invokation of this method, the singleton object is instantiated.

  protected:
    G4ProductionCutsTable();
  private:
    G4ProductionCutsTable(const G4ProductionCutsTable& right);

  public:
    virtual ~G4ProductionCutsTable();

  public: // with description
    void UpdateCoupleTable();
    // This method triggers an update of the table of G4ProductionCuts objects.

    void SetEnergyRange(G4double lowedge, G4double highedge);
    // This method sets the limits of energy cuts for all particles.

    G4double GetLowEdgeEnergy() const;
    G4double GetHighEdgeEnergy() const;
    // These methods get the limits of energy cuts for all particles.

    void DumpCouples() const;
    // Display a list of registored couples

  private:

   static G4ProductionCutsTable* fG4ProductionCutsTable;

   typedef std::vector<G4MaterialCutsCouple*> G4CoupleTable;
   typedef std::vector<G4MaterialCutsCouple*>::const_iterator CoupleTableIterator;
   typedef std::vector<G4double> G4CutVectorForAParticle;
   typedef std::vector<G4CutVectorForAParticle*> G4CutTable;
   G4CoupleTable coupleTable;
   G4CutTable rangeCutTable;
   G4CutTable energyCutTable;

   G4RegionStore* fG4RegionStore;
   G4VRangeToEnergyConverter* converters[NumberOfG4CutIndex]; 

   G4ProductionCuts* defaultProductionCuts;

// These two vectors are for the backward comparibility
   G4double* rangeDoubleVector[NumberOfG4CutIndex];
   G4double* energyDoubleVector[NumberOfG4CutIndex];

  public: 
   const std::vector<G4double>* GetRangeCutsVector(size_t pcIdx) const;
   const std::vector<G4double>* GetEnergyCutsVector(size_t pcIdx) const;

// These two vectors are for the backward comparibility
   G4double* GetRangeCutsDoubleVector(size_t pcIdx) const;
   G4double* GetEnergyCutsDoubleVector(size_t pcIdx) const;
  
  public: // with description
   size_t GetTableSize() const;
     // This method returns the size of the couple table.

   const G4MaterialCutsCouple* GetMaterialCutsCouple(G4int i) const;
    // This method returns the pointer to the couple.

   const G4MaterialCutsCouple*
     GetMaterialCutsCouple(const G4Material* aMat,
                           const G4ProductionCuts* aCut) const;
    // This method returns the pointer to the couple.

   G4int GetCoupleIndex(const G4MaterialCutsCouple* aCouple) const;
   G4int GetCoupleIndex(const G4Material* aMat,
                           const G4ProductionCuts* aCut) const;
    // These methods return the index of the couple.
    // -1 is returned if index is not found.

   G4bool IsModified() const;
    // This method returns TRUE if at least one production cut value is modified.
  
   void PhysicsTableUpdated();
    // This method resets the status of IsModified(). This method must
    // be exclusively used by RunManager when physics tables are built.

   G4ProductionCuts* GetDefaultProductionCuts() const;
    // This method returns the default production cuts.
   
  private:
    void ScanAndSetCouple(G4LogicalVolume* aLV,
			  G4MaterialCutsCouple* aCouple,
			  G4Region* aRegion);

    bool IsCoupleUsedInTheRegion(const G4MaterialCutsCouple* aCouple,
				 const G4Region* aRegion) const;

  public: // with description
   // Store cuts and material information in files under the specified directory.
  G4bool  StoreCutsTable(const G4String& directory, 
			 G4bool          ascii = false);
  
  // Retrieve cuts values information in files under the specified directory.
  G4bool  RetrieveCutsTable(const G4String& directory,
			    G4bool          ascii = false);
  
  // check stored material and cut values are consistent with the current detector setup. 
  G4bool CheckForRetrieveCutsTable(const G4String& directory, 
				   G4bool          ascii = false);

  protected:

  // Store material information in files under the specified directory.
  virtual G4bool  StoreMaterialInfo(const G4String& directory, 
				    G4bool          ascii = false);

  // check stored material is consistent with the current detector setup. 
  virtual G4bool  CheckMaterialInfo(const G4String& directory, 
				    G4bool          ascii = false);
  
  // Store materialCutsCouple information in files under the specified directory.
  virtual G4bool  StoreMaterialCutsCoupleInfo(const G4String& directory, 
				    G4bool          ascii = false);
  // check stored materialCutsCouple is consistent with the current detector setup. 
  virtual G4bool  CheckMaterialCutsCoupleInfo(const G4String& directory, 
				    G4bool          ascii = false);

  // Store cut values information in files under the specified directory.
  virtual G4bool  StoreCutsInfo(const G4String& directory, 
                                G4bool          ascii = false);
  
  // Retrieve cut values information in files under the specified directory.
  virtual G4bool  RetrieveCutsInfo(const G4String& directory,
                                   G4bool          ascii = false);

  private:
   G4bool firstUse;
   enum { FixedStringLengthForStore = 32 }; 

  public: // with description  
      void  SetVerboseLevel(G4int value);
      G4int GetVerboseLevel() const;
      // controle flag for output message
      //  0: Silent
      //  1: Warning message
      //  2: More

  private:
   G4int verboseLevel;

};

inline 
 const std::vector<G4double>* G4ProductionCutsTable::GetRangeCutsVector(size_t pcIdx) const
{ 
  return rangeCutTable[pcIdx]; 
}

inline 
 const std::vector<G4double>* G4ProductionCutsTable::GetEnergyCutsVector(size_t pcIdx) const
{
 return energyCutTable[pcIdx]; 
}

inline 
 size_t G4ProductionCutsTable::GetTableSize() const
{
  return coupleTable.size(); 
}

inline 
 const G4MaterialCutsCouple* G4ProductionCutsTable::GetMaterialCutsCouple(G4int i) const
{ 
  return coupleTable[size_t(i)]; 
}

inline 
 G4bool G4ProductionCutsTable::IsModified() const
{ 
  if(firstUse) return true;
  for(G4ProductionCutsTable::CoupleTableIterator itr=coupleTable.begin();
    itr!=coupleTable.end();itr++){ 
    if((*itr)->IsRecalcNeeded())
    {
      return true; 
    }
  }
  return false;
}
  
inline 
 void G4ProductionCutsTable::PhysicsTableUpdated()
{ 
  for(G4ProductionCutsTable::CoupleTableIterator itr=coupleTable.begin();itr!=coupleTable.end();itr++){ 
    (*itr)->PhysicsTableUpdated(); 
  }
}

inline
 G4double* G4ProductionCutsTable::GetRangeCutsDoubleVector(size_t pcIdx) const
{ return rangeDoubleVector[pcIdx]; }

inline
 G4double* G4ProductionCutsTable::GetEnergyCutsDoubleVector(size_t pcIdx) const
{ return energyDoubleVector[pcIdx]; }

inline
 G4ProductionCuts* G4ProductionCutsTable::GetDefaultProductionCuts() const
{ return defaultProductionCuts; }

inline
bool G4ProductionCutsTable::IsCoupleUsedInTheRegion(
                                 const G4MaterialCutsCouple* aCouple,
                                 const G4Region* aRegion) const
{
  G4ProductionCuts* fProductionCut = aRegion->GetProductionCuts();
  std::vector<G4Material*>::const_iterator mItr = aRegion->GetMaterialIterator();
  size_t nMaterial = aRegion->GetNumberOfMaterials();
  for(size_t iMate=0;iMate<nMaterial;iMate++, mItr++){
    if(aCouple->GetMaterial()==(*mItr) &&
       aCouple->GetProductionCuts()==fProductionCut){
      return true;
    }
  }
  return false;
}

inline
const G4MaterialCutsCouple* 
     G4ProductionCutsTable::GetMaterialCutsCouple(const G4Material* aMat,
                           const G4ProductionCuts* aCut) const
{ 
  for(CoupleTableIterator cItr=coupleTable.begin();cItr!=coupleTable.end();cItr++)
  { 
    if((*cItr)->GetMaterial()!=aMat) continue;
    if((*cItr)->GetProductionCuts()==aCut) return (*cItr);
  }
  return 0;
}

inline
G4int G4ProductionCutsTable::GetCoupleIndex(const G4MaterialCutsCouple* aCouple) const
{ 
  G4int idx = 0;
  for(CoupleTableIterator cItr=coupleTable.begin();cItr!=coupleTable.end();cItr++)
  {
    if((*cItr)==aCouple) return idx;
    idx++;
  }
  return -1;
}

inline
G4int G4ProductionCutsTable:: GetCoupleIndex(const G4Material* aMat,
                           const G4ProductionCuts* aCut) const
{
  const G4MaterialCutsCouple* aCouple = GetMaterialCutsCouple(aMat,aCut);
  return GetCoupleIndex(aCouple);
}

inline 
 void G4ProductionCutsTable::SetVerboseLevel(G4int value)
{
   verboseLevel = value;
}

inline 
 G4int G4ProductionCutsTable::GetVerboseLevel() const
{
   return verboseLevel;
}


#endif





