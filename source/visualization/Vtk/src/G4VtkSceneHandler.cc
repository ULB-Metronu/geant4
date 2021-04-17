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

#include "G4VtkSceneHandler.hh"

#include "G4PhysicalVolumeModel.hh"
#include "G4LogicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polyhedron.hh"
#include "G4UnitsTable.hh"

G4int G4VtkSceneHandler::fSceneIdCount = 0;
// Counter for XXX scene handlers.

G4VtkSceneHandler::G4VtkSceneHandler(G4VGraphicsSystem& system,
                                     const G4String& name) :
                                     G4VSceneHandler(system, fSceneIdCount++, name)
{}

G4VtkSceneHandler::~G4VtkSceneHandler() {}

#ifdef G4VTKDEBUG
void G4VtkSceneHandler::PrintThings() {
  G4cout << "  with transformation " << fObjectTransformation.xx() << G4endl;
  if (fpModel) {
    G4cout << " from " << fpModel->GetCurrentDescription()
           << " (tag " << fpModel->GetCurrentTag()
           << ')';
  }
  else {
    G4cout << "(not from a model)";
  }
  G4PhysicalVolumeModel* pPVModel = dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  if (pPVModel) {
    G4cout << "\n  current physical volume:        " << pPVModel->GetCurrentPV()->GetName()
           << "\n  current logical volume :        " << pPVModel->GetCurrentLV()->GetName() // There might be a problem with the LV pointer if this is a G4LogicalVolumeModel
           << "\n  current depth of geometry tree: " << pPVModel->GetCurrentDepth();
  }
  G4cout << G4endl;
}
#endif

void G4VtkSceneHandler::AddPrimitive(const G4Polyline& polyline) {
#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyline& polyline) called." << G4endl;
  // PrintThings();
#endif

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  const size_t nLines = polyline.size();

  for(size_t i=0; i < nLines; ++i) {
    points->InsertNextPoint(polyline[i].x(),
                            polyline[i].y(),
                            polyline[i].z());
  }

  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

  for (unsigned int i = 0; i < nLines-1; i++)
  {
    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, i);
    line->GetPointIds()->SetId(1, i+1);
    lines->InsertNextCell(line);
  }

  // Create a polydata to store everything in
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();

  // Add the points to the dataset
  polyData->SetPoints(points);

  // Add the lines to the dataset
  polyData->SetLines(lines);

  // Setup actor and mapper
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(polyData);

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  vtkSmartPointer<vtkMatrix4x4> transform = vtkSmartPointer<vtkMatrix4x4>::New();
  transform->SetElement(0,0,fObjectTransformation.xx());
  transform->SetElement(0,1,fObjectTransformation.xy());
  transform->SetElement(0,2,fObjectTransformation.xz());

  transform->SetElement(1,0,fObjectTransformation.yx());
  transform->SetElement(1,1,fObjectTransformation.yy());
  transform->SetElement(1,2,fObjectTransformation.yz());

  transform->SetElement(2,0,fObjectTransformation.zx());
  transform->SetElement(2,1,fObjectTransformation.zy());
  transform->SetElement(2,2,fObjectTransformation.zz());

  transform->SetElement(0,3,fObjectTransformation.dx());
  transform->SetElement(1,3,fObjectTransformation.dy());
  transform->SetElement(2,3,fObjectTransformation.dz());
  transform->SetElement(3,3,1.);

  // actor->SetUserMatrix(transform);

  // Setup actor and mapper
  vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();
  actor->GetProperty()->SetLineWidth(4);
  actor->GetProperty()->SetColor(colors->GetColor3d("Peacock").GetData());
  actor->SetVisibility(true);
  actor->GetProperty()->BackfaceCullingOff();
  actor->GetProperty()->FrontfaceCullingOff();

  G4VtkViewer* pVtkViewer = dynamic_cast<G4VtkViewer*>(fpViewer);
  pVtkViewer->renderer->AddActor(actor);

}

void G4VtkSceneHandler::AddPrimitive(const G4Text& text) {

  MarkerSizeType sizeType;
  GetMarkerSize(text,sizeType);
  const G4VisAttributes *fpVisAttribs = text.GetVisAttributes();
  const G4Colour &c = GetTextColour(text);
  G4double opacity = c.GetAlpha();

#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Text& text) called " << text.GetText() << " " << sizeType << G4endl;
  // G4cout << text << G4endl;
  // PrintThings();
  switch (sizeType) {
    default:
      G4cout << "default" << G4endl;
      break;
    case screen:
      // Draw in screen coordinates.
      G4cout << "screen" << G4endl;
      break;
    case world:
      G4cout << "world" << G4endl;
      break;
  }

#endif


  auto position = fObjectTransformation*G4Translate3D(text.GetPosition());

  double x = text.GetPosition().x();
  double y = text.GetPosition().y();
  double z = text.GetPosition().z();

  vtkSmartPointer<vtkBillboardTextActor3D> actor = vtkSmartPointer<vtkBillboardTextActor3D>::New();
  actor->SetInput(text.GetText().c_str());
  actor->SetPosition(x,y,z);
  actor->GetTextProperty()->SetFontSize (text.GetScreenSize());
  actor->GetTextProperty()->SetColor(c.GetRed(), c.GetBlue(), c.GetGreen());
  actor->GetTextProperty()->SetOpacity(opacity);

  G4VtkViewer* pVtkViewer = dynamic_cast<G4VtkViewer*>(fpViewer);
  pVtkViewer->renderer->AddActor(actor);
}

void G4VtkSceneHandler::AddPrimitive(const G4Circle& circle) {
#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Circle& circle) called." << G4endl;
  // PrintThings();
#endif

  MarkerSizeType sizeType;
  G4double size = GetMarkerSize (circle, sizeType);


  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->InsertNextPoint(circle.GetPosition().x(),
                          circle.GetPosition().y(),
                          circle.GetPosition().z());

  vtkNew<vtkPolyData> polydata;
  polydata->SetPoints(points);


  // Draw in world coordinates.
  vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
  polygonSource->SetNumberOfSides(50);
  polygonSource->SetRadius(0.5);

  auto position = fObjectTransformation*G4Translate3D(circle.GetPosition());

  vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
  glyph3D->SetSourceConnection(polygonSource->GetOutputPort());
  glyph3D->SetInputData(polydata);
  glyph3D->Update();

  // Visualize
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(glyph3D->GetOutputPort());;

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  G4VtkViewer* pVtkViewer = dynamic_cast<G4VtkViewer*>(fpViewer);
  pVtkViewer->renderer->AddActor(actor);

  switch (sizeType) {
    default:
    case screen:
      // Draw in screen coordinates.
      break;
    case world:
      // Draw in world coordinates.
      vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
      polygonSource->SetNumberOfSides(50);
      polygonSource->SetRadius(0.5);

      auto position = fObjectTransformation*G4Translate3D(circle.GetPosition());

      G4cout << circle.GetPosition().x() << " " << circle.GetPosition().y() << " " << circle.GetPosition().z() << G4endl;
      polygonSource->SetCenter(circle.GetPosition().x(),
                               circle.GetPosition().y(),
                               circle.GetPosition().z());

      // Visualize
      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper->SetInputConnection(polygonSource->GetOutputPort());;

      vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
      actor->SetMapper(mapper);

      G4VtkViewer* pVtkViewer = dynamic_cast<G4VtkViewer*>(fpViewer);
      pVtkViewer->renderer->AddActor(actor);

  }
}

void G4VtkSceneHandler::AddPrimitive(const G4Square& square) {
#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Square& square) called" << G4endl;
  //PrintThings();
#endif

  MarkerSizeType sizeType;
  G4double size = GetMarkerSize (square, sizeType);

  // Draw in world coordinates.
  vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
  polygonSource->SetNumberOfSides(4);
  polygonSource->SetRadius(5);

  auto position = fObjectTransformation*G4Translate3D(square.GetPosition());

  G4cout << square.GetPosition().x() << " " << square.GetPosition().y() << " " << square.GetPosition().z() << G4endl;
  polygonSource->SetCenter(square.GetPosition().x(),
                           square.GetPosition().y(),
                           square.GetPosition().z());

  // Visualize
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(polygonSource->GetOutputPort());;

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  G4VtkViewer* pVtkViewer = dynamic_cast<G4VtkViewer*>(fpViewer);
  pVtkViewer->renderer->AddActor(actor);

  switch (sizeType) {
    default:
    case screen:
      // Draw in screen coordinates.
      G4cout << "screen";
      break;
    case world:
      // Draw in world coordinates.
      G4cout << "world";
      break;
  }
}

void G4VtkSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) {
#ifdef G4VTKDEBUG
  G4cout << "G4VtkSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) called." << G4endl;
  // PrintThings();
#endif

  //Get colour, etc..
  if (polyhedron.GetNoFacets() == 0) return;

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA = fpViewer -> GetApplicableVisAttributes(polyhedron.GetVisAttributes ());
  G4Color colour    = pVA->GetColour();
  G4bool  isVisible = pVA->IsVisible();

  // Get view parameters that the user can force through the vis attributes, thereby over-riding the current view parameter.
  G4ViewParameters::DrawingStyle drawing_style = GetDrawingStyle(pVA);
  //G4bool isAuxEdgeVisible = GetAuxEdgeVisible (pVA);

  vtkSmartPointer<vtkPolyData>  polydata  = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints>    points    = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> polys     = vtkSmartPointer<vtkCellArray>::New();

  G4bool notLastFace;
  int    iVert  = 0;
  do {
    G4Point3D  vertex[4];
    G4int      edgeFlag[4];
    G4Normal3D normals[4];
    G4int      nEdges;
    notLastFace = polyhedron.GetNextFacet(nEdges, vertex, edgeFlag, normals);

    vtkSmartPointer<vtkIdList> poly = vtkSmartPointer<vtkIdList>::New();
    // loop over vertices
    for(int i=0; i < nEdges; i++) {
      points->InsertNextPoint(vertex[i].x(), vertex[i].y(), vertex[i].z());
      poly->InsertNextId(iVert);
      iVert++;
    }
    polys->InsertNextCell(poly);

  } while (notLastFace);

  polydata->SetPoints(points);
  polydata->SetPolys(polys);

  /*
   */
  vtkSmartPointer<vtkFeatureEdges> featureEdges = vtkSmartPointer<vtkFeatureEdges>::New();
  featureEdges->SetInputData(polydata);
  featureEdges->SetColoring(true);
  featureEdges->BoundaryEdgesOff();
  featureEdges->NonManifoldEdgesOff();
  featureEdges->ManifoldEdgesOff();
  featureEdges->FeatureEdgesOn();
  featureEdges->SetFeatureAngle(1);
  featureEdges->Update();

  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(polydata);
  // mapper->SetInputConnection(featureEdges->GetOutputPort());
  // mapper->SetScalarModeToUseCellData();

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  vtkSmartPointer<vtkMatrix4x4> transform = vtkSmartPointer<vtkMatrix4x4>::New();
  transform->SetElement(0,0,fObjectTransformation.xx());
  transform->SetElement(0,1,fObjectTransformation.xy());
  transform->SetElement(0,2,fObjectTransformation.xz());

  transform->SetElement(1,0,fObjectTransformation.yx());
  transform->SetElement(1,1,fObjectTransformation.yy());
  transform->SetElement(1,2,fObjectTransformation.yz());

  transform->SetElement(2,0,fObjectTransformation.zx());
  transform->SetElement(2,1,fObjectTransformation.zy());
  transform->SetElement(2,2,fObjectTransformation.zz());

  transform->SetElement(0,3,fObjectTransformation.dx());
  transform->SetElement(1,3,fObjectTransformation.dy());
  transform->SetElement(2,3,fObjectTransformation.dz());
  transform->SetElement(3,3,1.);

  actor->SetUserMatrix(transform);

  actor->GetProperty()->SetColor(colour.GetRed(), colour.GetGreen(), colour.GetBlue());
  actor->GetProperty()->SetOpacity(colour.GetAlpha());
  actor->SetVisibility(isVisible);

  G4VtkViewer* pVtkViewer = dynamic_cast<G4VtkViewer*>(fpViewer);
  pVtkViewer->renderer->AddActor(actor);


  // Initial action depending on drawing style.
  switch (drawing_style) {
  case (G4ViewParameters::hsr):
    {
      break;
    }
  case (G4ViewParameters::hlr):
    {
      break;
    }
  case (G4ViewParameters::wireframe):
    {
      break;
    }
  default:
    {
      break;
    }     
  }
}
