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
//
// 
// John Allison  5th April 2001
// A template for a simplest possible graphics driver.
//?? Lines or sections marked like this require specialisation for your driver.

#ifndef G4VTKSCENEHANDLER_HH
#define G4VTKSCENEHANDLER_HH


#include "G4VSceneHandler.hh"

#include <map>
#include <vector>

#include "G4VtkViewer.hh"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra-semi"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkLine.h"
#include "vtkNamedColors.h"
#include "vtkProperty.h"
#include "vtkTextProperty.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkFollower.h"
#include "vtkTextActor.h"
#include "vtkBillboardTextActor3D.h"
#include "vtkRegularPolygonSource.h"
#include "vtkSphereSource.h"
#include "vtkFeatureEdges.h"
#include "vtkGlyph2D.h"
#include "vtkGlyph3D.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkMatrix4x4.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkTriangleFilter.h"
#include "vtkTensorGlyph.h"

#pragma GCC diagnostic pop

class G4VtkSceneHandler: public G4VSceneHandler {
  friend class G4VtkViewer;

public:
  G4VtkSceneHandler(G4VGraphicsSystem& system, const G4String& name);
  virtual ~G4VtkSceneHandler();

  ////////////////////////////////////////////////////////////////
  // Required implementation of pure virtual functions...

  void AddPrimitive(const G4Polyline&);
  void AddPrimitive(const G4Text&);
  void AddPrimitive(const G4Circle&);
  void AddPrimitive(const G4Square&);
  void AddPrimitive(const G4Polyhedron&);
  // Further optional AddPrimitive methods.  Explicitly invoke base
  // class methods if not otherwise defined to avoid warnings about
  // hiding of base class methods.
  void AddPrimitive(const G4Polymarker& polymarker)
  {
#ifdef G4VTKDEBUG
    G4cout << "=================================" << G4endl;
    G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polymarker& polymarker) called." << G4endl;
#endif
    G4VSceneHandler::AddPrimitive (polymarker);
  }
  void AddPrimitive(const G4Scale& scale)
  {
#ifdef G4VTKDEBUG
    G4cout << "=================================" << G4endl;
    G4cout << "G4VtkSceneHandler::AddScale(const G4Scale& scale) called." << G4endl;
#endif
    G4VSceneHandler::AddPrimitive(scale);
  }

  void Modified() {
    for (auto it = polylineDataMap.begin(); it != polylineDataMap.end(); it++)
      it->second->Modified();
    for (auto it = polylineLineMap.begin(); it != polylineLineMap.end(); it++)
      it->second->Modified();

    for (auto it = circleDataMap.begin(); it != circleDataMap.end(); it++) {
      it->second->Modified();
    }

    for (auto it = squareDataMap.begin(); it != squareDataMap.end(); it++) {
      it->second->Modified();
    }

    for(auto it = instancePositionMap.begin(); it != instancePositionMap.end(); it++) {
      it->second->Modified();
    }

    for(auto it = instanceRotationMap.begin(); it != instanceRotationMap.end(); it++) {
      it->second->Modified();
    }

    for(auto it = instanceActorMap.begin(); it != instanceActorMap.end(); it++) {
      it->second->Modified();
    }

    for(auto it = tgMap.begin(); it != tgMap.end(); it++) {
      it->second->Update();
    }


    G4cout << "G4VtkSceneHandler::Modified()    polyline styles: " << polylineVisAttributesMap.size() << G4endl;
    for (auto it = polylineLineMap.begin(); it != polylineLineMap.end(); it++)
      G4cout << "G4VtkSceneHandler::Modified()   polyline segments: " << it->second->GetNumberOfCells() << G4endl;

    G4cout << "G4VtkSceneHandler::Modified()       circle styles: " << circleVisAttributesMap.size() << G4endl;
    for (auto it = circleDataMap.begin(); it != circleDataMap.end(); it++)
      G4cout << "G4VtkSceneHandler::Modified()             circles: " << it->second->GetNumberOfPoints() << G4endl;

    G4cout << "G4VtkSceneHandler::Modified()      square styles: " << squareVisAttributesMap.size() << G4endl;
    for (auto it = squareDataMap.begin(); it != squareDataMap.end(); it++)
      G4cout << "G4VtkSceneHanler::Modified()           squares: " << it->second->GetNumberOfPoints() << G4endl;

    G4cout << "G4VtkSceneHandler::Modified()   unique polyhedra: " << polyhedronDataMap.size() << G4endl;

    for (auto it = polyhedronPolyDataMap.begin(); it != polyhedronPolyDataMap.end(); it++) {
      G4cout << "G4VtkSceneHandler::Modified()  polyhedronPolyData: " << it->second->GetPoints()->GetNumberOfPoints() << " " << it->second->GetPolys()->GetNumberOfCells() << " " << polyhedronPolyDataCountMap[it->first] <<G4endl;
    }

   // G4cout << "G4VtkSceneHandler::Modified()    placed polyhedra: " << polyhedronActorVector.size() << G4endl;


  }

  void Clear() {
    polylineVisAttributesMap.clear();
    polylineDataMap.clear();
    polylineLineMap.clear();
    polylinePolyDataMap.clear();
    polylinePolyDataMapperMap.clear();
    polylinePolyDataActorMap.clear();

    for (auto it = polylineDataMap.begin(); it != polylineDataMap.end(); it++)
      it->second->Reset();
    for (auto it = polylineLineMap.begin(); it != polylineLineMap.end(); it++)
      it->second->Reset();

    circleVisAttributesMap.clear();
    circleDataMap.clear();
    circlePolyDataMap.clear();
    circleFilterMap.clear();
    circlePolyDataMapperMap.clear();
    circlePolyDataActorMap.clear();

    for (auto it = squareDataMap.begin(); it != squareDataMap.end(); it++)
      it->second->Reset();

    squareVisAttributesMap.clear();
    squareDataMap.clear();
    squarePolyDataMap.clear();
    squareFilterMap.clear();
    squarePolyDataMapperMap.clear();
    squarePolyDataActorMap.clear();

    for (auto it = squareDataMap.begin(); it != squareDataMap.end(); it++)
      it->second->Reset();

    polyhedronVisAttributesMap.clear();
    polyhedronPolyDataMap.clear();

    // polyhedronMapperMap.clear();
    // polyhedronActorVector.clear();

  }

protected:
  static G4int         fSceneIdCount;  // Counter for Vtk scene handlers.

  // maps for polylines
  std::map<std::size_t, const G4VisAttributes*>             polylineVisAttributesMap;
  std::map<std::size_t, vtkSmartPointer<vtkPoints>>         polylineDataMap;
  std::map<std::size_t, vtkSmartPointer<vtkCellArray>>      polylineLineMap;
  std::map<std::size_t, vtkSmartPointer<vtkPolyData>>       polylinePolyDataMap;
  std::map<std::size_t, vtkSmartPointer<vtkPolyDataMapper>> polylinePolyDataMapperMap;
  std::map<std::size_t, vtkSmartPointer<vtkActor>>          polylinePolyDataActorMap;

  // maps for circles
  std::map<std::size_t, const G4VisAttributes*>                circleVisAttributesMap;
  std::map<std::size_t, vtkSmartPointer<vtkPoints>>            circleDataMap;
  std::map<std::size_t, vtkSmartPointer<vtkPolyData>>          circlePolyDataMap;
  std::map<std::size_t, vtkSmartPointer<vtkVertexGlyphFilter>> circleFilterMap;
  std::map<std::size_t, vtkSmartPointer<vtkPolyDataMapper>>    circlePolyDataMapperMap;
  std::map<std::size_t, vtkSmartPointer<vtkActor>>             circlePolyDataActorMap;

  // maps for squares
  std::map<std::size_t, const G4VisAttributes*>                squareVisAttributesMap;
  std::map<std::size_t, vtkSmartPointer<vtkPoints>>            squareDataMap;
  std::map<std::size_t, vtkSmartPointer<vtkPolyData>>          squarePolyDataMap;
  std::map<std::size_t, vtkSmartPointer<vtkVertexGlyphFilter>> squareFilterMap;
  std::map<std::size_t, vtkSmartPointer<vtkPolyDataMapper>>    squarePolyDataMapperMap;
  std::map<std::size_t, vtkSmartPointer<vtkActor>>             squarePolyDataActorMap;

  // map for polyhedra
  std::map<std::size_t, const G4VisAttributes*>                polyhedronVisAttributesMap;
  std::map<std::size_t, vtkSmartPointer<vtkPoints>>            polyhedronDataMap;
  std::map<std::size_t, vtkSmartPointer<vtkCellArray>>         polyhedronPolyMap;
  std::map<std::size_t, vtkSmartPointer<vtkPolyData>>          polyhedronPolyDataMap;
  std::map<std::size_t, std::size_t>                           polyhedronPolyDataCountMap;

  //std::map<std::size_t, vtkSmartPointer<vtkPolyDataMapper>>    polyhedronMapperMap;
  //std::map<std::size_t, vtkSmartPointer<vtkActor>>             polyhedronActorMap;

  std::map<std::size_t, vtkSmartPointer<vtkPoints>>            instancePositionMap;
  std::map<std::size_t, vtkSmartPointer<vtkDoubleArray>>       instanceRotationMap;
  std::map<std::size_t, vtkSmartPointer<vtkActor>>             instanceActorMap;
  std::map<std::size_t, vtkSmartPointer<vtkTensorGlyph>>       tgMap;

private:
#ifdef G4VTKDEBUG
  void PrintThings();
#endif

};

#endif
