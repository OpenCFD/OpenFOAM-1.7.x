/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::vtkFoam

Description

SourceFiles
    vtkFoam.C
    vtkFoamInsertNextPoint.H
    vtkFoamAddFields.H
    vtkFoamAddInternalMesh.H
    vtkFoamConvertFields.H
    vtkFoamConvertVolField.H
    vtkFoamConvertPatchFaceField.H
    vtkFoamConvertPointField.H
    vtkFoamConvertPatchPointField.H

\*---------------------------------------------------------------------------*/

#ifndef vtkFoam_H
#define vtkFoam_H

#include "className.H"
#include "fileName.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

// VTK class forward declarations
class vtkFoamReader;
class vtkUnstructuredGrid;
class vtkPoints;
class vtkDataArraySelection;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Foam class forward declarations
class argList;
class Time;
class fvMesh;
class IOobjectList;
class polyPatch;


/*---------------------------------------------------------------------------*\
                           Class vtkFoam Declaration
\*---------------------------------------------------------------------------*/

class vtkFoam
{
    // Private data

        //- Access to the controlling vtkFoamReader
        vtkFoamReader *reader_;

        autoPtr<argList> argsPtr_;
        autoPtr<Time> dbPtr_;
        fvMesh* meshPtr_;

        //- Selected regions, [0] = internal mesh, [1-nPatches] = patches
        boolList selectedRegions_;

        //- Lables of cell-centres used as additional points when decomposing
        // polyhedra
        labelList addPointCellLabels_;

        //- Label of original cell the decomposed cells are split from
        labelList superCells_;


    // Private Member Functions

        //- Pad-out the time name to avoid bug in the GUI redraw
        static string padTimeString(const string&);

        //- Pad-out the patch name
        static string padPatchString(const string&);

        //- Find and set the selected time from all the methods of selection
        static void setSelectedTime
        (
            Time& runTime,
            vtkFoamReader* reader
        );

        //- Update the selected regions
        void updateSelectedRegions();

        //- Convert the mesh according to the list of selected regions
        void convertMesh();

        //- Add the internal mesh to the set of Outputs if selected
        void addInternalMesh(const fvMesh&, vtkUnstructuredGrid*);

        //- Add the internal patch to the set of Outputs if selected
        void addPatch(const polyPatch&, vtkUnstructuredGrid*);

        //- Add the fields in th selested time directory to the selection lists
        template<class GeoField>
        void addFields
        (
            vtkDataArraySelection *fieldSelection,
            const IOobjectList& objects
        );

        //- Convert the selected volFields
        template<class Type>
        void convertVolFields
        (
            const fvMesh& mesh,
            const volPointInterpolation& pInterp,
            const IOobjectList& objects,
            vtkDataArraySelection *fieldSelection
        );

        template<class Type>
        void convertVolField
        (
            const GeometricField<Type, fvPatchField, volMesh>& tf
        );

        template<class Type>
        void convertPatchFaceField
        (
            const word& name,
            const Field<Type>& tf,
            const label regioni
        );

        //- Convert the selected pointFields
        template<class Type>
        void convertPointFields
        (
            const fvMesh& mesh,
            const IOobjectList& objects,
            vtkDataArraySelection *fieldSelection
        );

        template<class Type>
        void convertPointField
        (
            const GeometricField<Type, pointPatchField, pointMesh>& ptf,
            const GeometricField<Type, fvPatchField, volMesh>& tf
        );

        template<class Type>
        void convertPatchPointField
        (
            const word& name,
            const Field<Type>& tf,
            const label regioni
        );

        //- Set the name of the Output vtkUnstructuredGrid
        void SetName(vtkUnstructuredGrid *vtkMesh, const char* name);


        //- Disallow default bitwise copy construct
        vtkFoam(const vtkFoam&);

        //- Disallow default bitwise assignment
        void operator=(const vtkFoam&);


public:

    // Static data members

        ClassName("vtkFoam");


    // Constructors

        //- Construct from components
        vtkFoam(const char* const FileName, vtkFoamReader* reader);


    // Destructor

        ~vtkFoam();


    // Member Functions

        void UpdateInformation();
        void Update();
};


// * * * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * //

template<>
void vtkFoam::convertVolField
(
    const GeometricField<scalar, fvPatchField, volMesh>& sf
);

template<>
void vtkFoam::convertPatchFaceField
(
    const word& name,
    const Field<scalar>& sf,
    const label regioni
);

template<>
void vtkFoam::convertPointField
(
    const GeometricField<scalar, pointPatchField, pointMesh>& psf,
    const GeometricField<scalar, fvPatchField, volMesh>& sf
);

template<>
void vtkFoam::convertPatchPointField
(
    const word& name,
    const Field<scalar>& sf,
    const label regioni
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
