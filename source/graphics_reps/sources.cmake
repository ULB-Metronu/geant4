#------------------------------------------------------------------------------
# Module : G4graphics_reps
# Package: Geant4.src.G4graphics_reps
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#

geant4_define_module(NAME G4graphics_reps
  HEADERS
    G4AttCheck.hh
    G4AttDef.hh
    G4AttDefStore.hh
    G4AttDefT.hh
    G4AttHolder.hh
    G4AttUtils.hh
    G4AttValue.hh
    G4Circle.hh
    G4Color.hh
    G4Colour.hh
    G4ConversionFatalError.hh
    G4ConversionUtils.hh
    G4CreatorFactoryT.hh
    G4DimensionedDouble.hh
    G4DimensionedThreeVector.hh
    G4DimensionedType.hh
    G4PlacedPolyhedron.hh
    G4Point3DList.hh
    G4Polyhedron.hh
    G4PolyhedronArbitrary.hh
    G4Polyline.hh
    G4Polymarker.hh
    G4Polymarker.icc
    G4Scale.hh
    G4Scale.icc
    G4SmartFilter.hh
    G4Square.hh
    G4Square.icc
    G4Text.hh
    G4Text.icc
    G4TypeKey.hh
    G4TypeKeyT.hh
    G4VFilter.hh
    G4VGraphicsScene.hh
    G4VMarker.hh
    G4VMarker.icc
    G4VVisManager.hh
    G4VisAttributes.hh
    G4VisAttributes.icc
    G4VisExtent.hh
    G4Visible.hh
    G4Visible.icc
    HepPolyhedron.h
    HepPolyhedronProcessor.h
    graphics_reps_defs.hh
  SOURCES
    BooleanProcessor.src
    G4AttCheck.cc
    G4AttDef.cc
    G4AttDefStore.cc
    G4AttHolder.cc
    G4AttUtils.cc
    G4Circle.cc
    G4Colour.cc
    G4DimensionedTypeUtils.cc
    G4PlacedPolyhedron.cc
    G4Point3DList.cc
    G4Polyhedron.cc
    G4PolyhedronArbitrary.cc
    G4Polyline.cc
    G4Polymarker.cc
    G4Scale.cc
    G4Square.cc
    G4Text.cc
    G4VGraphicsScene.cc
    G4VMarker.cc
    G4VVisManager.cc
    G4VisAttributes.cc
    G4VisExtent.cc
    G4Visible.cc
    HepPolyhedron.cc
    HepPolyhedronProcessor.src
  GRANULAR_DEPENDENCIES
    G4globman
    G4intercoms
  GLOBAL_DEPENDENCIES
    G4global
    G4intercoms
  LINK_LIBRARIES
    ${GMP_LIBRARY}
    ${MPFR_LIBRARY}
)

# List any source specific properties here
geant4_module_include_directories(G4graphics_reps PUBLIC $<BUILD_INTERFACE:${CGAL_INCLUDE_DIRS}>)

if(GEANT4_USE_CGAL)
    geant4_module_compile_definitions(G4graphics_reps PUBLIC G4VIS_USE_CGAL)
endif()