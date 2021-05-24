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

#include "G4VtkViewer.hh"

#include "G4VSceneHandler.hh"
#include "G4VtkSceneHandler.hh"

G4VtkViewer::G4VtkViewer (G4VSceneHandler& sceneHandler, const G4String& name) :
G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name) {

  vtkObject::GlobalWarningDisplayOff();

  G4cout << "G4VtkViewer::G4VtkViewer" << G4endl;

  G4cout << "G4VtkViewer::G4VtkViewer> " << fVP.GetWindowSizeHintX() << " " << fVP.GetWindowSizeHintY() << G4endl;
  G4cout << "G4VtkViewer::G4VtkViewer> " << fVP.GetWindowLocationHintX() << " " << fVP.GetWindowLocationHintY() << G4endl;

  renderWindow->SetSize(fVP.GetWindowSizeHintX(), fVP.GetWindowSizeHintY());
  renderWindow->SetPosition(fVP.GetWindowLocationHintX(), fVP.GetWindowLocationHintY());
  renderWindow->SetWindowName("Vtk viewer");

  renderWindow->AddRenderer(renderer);
  renderWindowInteractor->SetRenderWindow(renderWindow);

  // TODO proper camera parameter settings
  camera->SetPosition(0, 0, 1000);
  camera->SetFocalPoint(0, 0, 0);
  renderer->SetActiveCamera(camera);

  // Set callback to match VTK parameters to Geant4
  geant4Callback->SetGeant4ViewParameters(&fVP);
  renderer->AddObserver(vtkCommand::EndEvent, geant4Callback);

  // vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  // renderWindowInteractor->SetInteractorStyle(style);
}

G4VtkViewer::~G4VtkViewer() {}

void G4VtkViewer::SetView() {
#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::SetView() called>" << G4endl;
  G4cout << fVP << G4endl;
#endif

  // culling
  // G4bool isCulling = fVP.IsCulling();
  // G4bool isCullingInvisible = fVP.IsCullingInvisible();
  // G4bool isDensityCulling = fVP.IsDensityCulling();
  // G4bool isCullingCovered = fVP.IsCullingCovered();

  //G4cout << "G4VtkViewer::SetView() called> isCulling:         " << isCulling << G4endl;
  //G4cout << "G4VtkViewer::SetView() called> isCullingInvisible:" << isCullingInvisible << G4endl;
  //G4cout << "G4VtkViewer::SetView() called> isDensityCulling:  " << isDensityCulling << G4endl;
  //G4cout << "G4VtkViewer::SetView() called> isCullingCovered:  " << isCullingCovered << G4endl;

  // background colour
  const G4Colour backgroundColour = fVP.GetBackgroundColour();
  G4cout << "G4VtkViewer::SetView() called>     setting backgroundColour: " << backgroundColour << G4endl;
  renderer->SetBackground(backgroundColour.GetRed(), backgroundColour.GetGreen(), backgroundColour.GetBlue());

  // Camera
  G4double radius = fSceneHandler.GetExtent().GetExtentRadius();
  if (radius <= 0.) radius = 1.;
  G4double      cameraDistance = fVP.GetCameraDistance(2 * radius);
  G4double      fieldHalfAngle = fVP.GetFieldHalfAngle();
  G4double          zoomFactor = fVP.GetZoomFactor();
  // G4double               dolly = fVP.GetDolly();
  G4Point3D viewpointDirection = fVP.GetViewpointDirection();
  G4Point3D        targetPoint = fVP.GetCurrentTargetPoint();

  G4Point3D cameraPosition = targetPoint + viewpointDirection.unit() * cameraDistance;

  //G4cout << "G4VtkViewer::SetView() called>     CameraDistance: " << cameraDistance << G4endl;
  //G4cout << "G4VtkViewer::SetView() called>         ZoomFactor: " << zoomFactor << G4endl;
  //G4cout << "G4VtkViewer::SetView() called>              Dolly: " << dolly << G4endl;
  //G4cout << "G4VtkViewer::SetView() called> ViewpointDirection: " << viewpointDirection.x() << " "
  //       << viewpointDirection.y() << " " << viewpointDirection.z() << G4endl;
  //G4cout << "G4VtkViewer::SetView() called> CurrentTargetPoint: " << targetPoint.x() << " " << targetPoint.y() << " "
  //       << targetPoint.z() << G4endl;

  // projection type
  if(fieldHalfAngle == 0) {
    G4cout << "G4VtkViewer::SetView() called>           setting orthogonal: " << G4endl;
    camera->SetParallelProjection(true);
  }
  else {
    G4cout << "G4VtkViewer::SetView() called>          setting perspective: " << G4endl;
    camera->SetParallelProjection(false);
  }

  // view angle
  if(fieldHalfAngle != 0) {
    G4cout << "G4VtkViewer::SetView() called>       setting fieldHalfAngle: " << fieldHalfAngle << G4endl;
    camera->SetViewAngle(2*fieldHalfAngle/M_PI*180);
  }

  // zoom factor
  camera->Zoom(zoomFactor);

  // target and camera positions
  G4cout << "G4VtkViewer::SetView() called>          setting targetPoint: " << targetPoint.x() << " "
                                                                            << targetPoint.y() << " "
                                                                            << targetPoint.z() << G4endl;
  camera->SetFocalPoint(targetPoint.x(), targetPoint.y(), targetPoint.z());
  G4cout << "G4VtkViewer::SetView() called>       setting cameraPosition: " << cameraPosition.x() << " "
                                                                            << cameraPosition.y() << " "
                                                                            << cameraPosition.z() << G4endl;
  camera->SetPosition(cameraPosition.x(), cameraPosition.y(), cameraPosition.z());

  // Light
  const G4Vector3D lightDirection = fVP.GetLightpointDirection();
  G4bool     lightsMoveWithCamera = fVP.GetLightsMoveWithCamera();
  G4Vector3D        lightPosition = targetPoint + lightDirection.unit() * cameraDistance;
  G4cout << "G4VtkViewer::SetView() called>       setting LightDirection: " << lightDirection.x() << " "
                                                                            << lightDirection.y() << " "
                                                                            << lightDirection.z() << G4endl;

  G4cout << "G4VtkViewer::SetView() called> setting lightsMoveWithCamera: " << lightsMoveWithCamera << G4endl;


  if (lightsMoveWithCamera) {
    G4cout << "G4VtkViewer::SetView() called> setting lightsMoveWithCamera: " << lightsMoveWithCamera << G4endl;
    light = vtkSmartPointer<vtkLight>::New();
    light->SetLightTypeToCameraLight();
    light->SetPosition(camera->GetPosition());
    renderer->RemoveAllLights();
    renderer->AddLight(light);
    renderer->SetLightFollowCamera(true);
    renderWindowInteractor->SetLightFollowCamera(true);
  }
  else {
    G4cout << "G4VtkViewer::SetView() called> setting lightsMoveWithCamera: " << lightsMoveWithCamera << G4endl;
    light = vtkSmartPointer<vtkLight>::New();
    light->SetLightTypeToSceneLight();
    light->SetPosition(lightDirection.unit()*cameraDistance*100);
    renderer->RemoveAllLights();
    renderer->AddLight(light);
    renderer->SetLightFollowCamera(false);
    renderWindowInteractor->SetLightFollowCamera(false);
  }

  // Rotation style
  const G4Vector3D upVector       = fVP.GetUpVector();
  G4cout << "G4VtkViewer::SetView() called>             setting UpVector: " << upVector.x() << " "
                                                                            << upVector.y() << " "
                                                                            << upVector.z() << G4endl;

  // Interaction (rotation style)
  G4ViewParameters::RotationStyle rotationStyle  = fVP.GetRotationStyle();
  G4cout << "G4VtkViewer::SetView() called>        setting rotationStyle: " << rotationStyle << G4endl;
  if (rotationStyle == G4ViewParameters::RotationStyle::freeRotation) {
    G4cout << "G4VtkViewer::SetView() called>         setting freeRotation: " << G4endl;
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    renderWindowInteractor->SetInteractorStyle(style);
  }
  else if(rotationStyle == G4ViewParameters::RotationStyle::constrainUpDirection) {
    G4cout << "G4VtkViewer::SetView() called> setting constrainedUpDirection: " << G4endl;
    // camera->SetViewUp(upVector.x(), upVector.y(), upVector.z());
    vtkSmartPointer<vtkInteractorStyleTerrain> style = vtkSmartPointer<vtkInteractorStyleTerrain>::New();
    renderWindowInteractor->SetInteractorStyle(style);
  }
}

void G4VtkViewer::ClearView() {
#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::ClearView() called." << G4endl;
#endif

  vtkActorCollection *actors = renderer->GetActors();
  vtkActor *actor = actors->GetLastActor();

  while(actor) {
    G4cout << "G4VtkViewer::ClearView() remove actor " << actor << G4endl;
    renderer->RemoveActor(actor);
    actor = actors->GetLastActor();
  }

  vtkPropCollection *props = renderer->GetViewProps();
  vtkProp *prop  = props->GetLastProp();

  while(prop) {
    G4cout << "G4VtkViewer::ClearView() remove prop " << prop << G4endl;
    renderer->RemoveViewProp(prop);
    prop = props->GetLastProp();
  }

  G4VtkSceneHandler& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  fVtkSceneHandler.Clear();

}

void G4VtkViewer::DrawView() {
#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::DrawView() called." << G4endl;
#endif

  // First, a view should decide when to re-visit the G4 kernel.
  // Sometimes it might not be necessary, e.g., if the scene is stored
  // in a graphical database (e.g., OpenGL's display lists) and only
  // the viewing angle has changed.  But graphics systems without a
  // graphical database will always need to visit the G4 kernel.

  NeedKernelVisit();  // Default is - always visit G4 kernel.
  // Note: this routine sets the fNeedKernelVisit flag of *all* the
  // views of the scene.

  ProcessView();      // The basic logic is here.

  // Then a view may have more to do, e.g., display the graphical
  // database.  That code should come here...

  // ...before finally...
  FinishView();       // Flush streams and/or swap buffers.
}

void G4VtkViewer::ShowView() {
#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::ShowView() called." << G4endl;
  // static_cast<G4VtkSceneHandler&>(fSceneHandler).PrintStores();
#endif

  G4VtkSceneHandler& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  fVtkSceneHandler.Modified();

  infoTextActor->GetTextProperty()->SetFontSize(28);
  G4Colour colour = fVP.GetBackgroundColour();

  // make sure text is always visible
  infoTextActor->GetTextProperty()->SetColor(fmod(colour.GetRed()   + 0.5, 1.0),
                                             fmod(colour.GetGreen() + 0.5, 1.0),
                                             fmod(colour.GetBlue()  + 0.5, 1.0));
  infoCallback->SetTextActor(infoTextActor);
  renderer->AddObserver(vtkCommand::EndEvent,infoCallback);
  renderer->AddActor(infoTextActor);

  renderWindow->Render();
  renderWindowInteractor->Initialize();
  renderWindowInteractor->Start();
}
