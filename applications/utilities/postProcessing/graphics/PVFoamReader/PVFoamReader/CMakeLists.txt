# Set the project/library name here.  Classes should be declared as
# "class VTK_FoamReader_EXPORT vtkFoo", where PVFoamReader is the name of the
# project set here.
PROJECT(PVFoamReader)

# Set your list of sources here.  Do not change the name of the
# PVFoamReader_SRCS variable.
SET(PVFoamReader_SRCS
  vtkFoamReader.cxx
  vtkFoamData.cxx
)

# List the kits from VTK that are needed by this project.  Do not
# change the name of the PVFoamReader_LIBS variable.
SET(PVFoamReader_LIBS
  vtkGraphics
)

INCLUDE_DIRECTORIES(
    $ENV{WM_PROJECT_DIR}/src/OpenFOAM/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{ParaView_INST_DIR}/include
    ../vtkFoam/lnInclude
)
ADD_DEFINITIONS(-DWM_$ENV{WM_PRECISION_OPTION})

#-----------------------------------------------------------------------------
# Most users should not need to change anything below this line.

# Need to include class headers and the configuration header.
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

# Import ParaView build settings.
FIND_PACKAGE(ParaView)
IF(NOT ParaView_FOUND)
  MESSAGE(FATAL_ERROR "ParaView is required.  Set ParaView_DIR.")
ENDIF(NOT ParaView_FOUND)
IF(NOT PARAVIEW_BUILD_SHARED_LIBS)
  MESSAGE(FATAL_ERROR "ParaView must be built with BUILD_SHARED_LIBS ON.")
ENDIF(NOT PARAVIEW_BUILD_SHARED_LIBS)
INCLUDE(${PARAVIEW_USE_FILE})

# Configure output directories.
SET (LIBRARY_OUTPUT_PATH $ENV{FOAM_LIBBIN} CACHE INTERNAL
     "Single output directory for building all libraries.")
SET (EXECUTABLE_OUTPUT_PATH $ENV{FOAM_APPBIN} CACHE INTERNAL
     "Single output directory for building all executables.")

# Create the configuration header.
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/PVFoamReaderConfigure.h.in
               ${PROJECT_BINARY_DIR}/vtk${PROJECT_NAME}Configure.h
               @ONLY IMMEDIATE)

# Must be defined before ADD_LIBRARY  for cmake version 2.3.4-20060317
LINK_DIRECTORIES($ENV{FOAM_LIBBIN})
LINK_DIRECTORIES($ENV{FOAM_USER_LIBBIN})

# Create vtk client/server wrappers for the classes.
VTK_WRAP_ClientServer(${PROJECT_NAME} PVFoamReaderCS_SRCS "${PVFoamReader_SRCS}")

# Build the package as a plugin for ParaView.
ADD_LIBRARY(${PROJECT_NAME} MODULE ${PVFoamReader_SRCS} ${PVFoamReaderCS_SRCS})
LINK_DIRECTORIES($ENV{ParaView_INST_DIR}/lib)
TARGET_LINK_LIBRARIES(${PROJECT_NAME} vtkClientServer)
FOREACH(c ${PVFoamReader_LIBS})
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${c}CS)
ENDFOREACH(c)

# Place the package configuration file into the build tree.
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/PVFoamReader.xml.in
               ${PROJECT_BINARY_DIR}/../${PROJECT_NAME}.xml @ONLY IMMEDIATE)

# Place the package configuration file into the build tree.
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/PVFoamReader.pvsm.in
               ${PROJECT_BINARY_DIR}/../${PROJECT_NAME}.pvsm @ONLY IMMEDIATE)

TARGET_LINK_LIBRARIES(${PROJECT_NAME} OpenFOAM finiteVolume vtkFoam)
#-----------------------------------------------------------------------------
