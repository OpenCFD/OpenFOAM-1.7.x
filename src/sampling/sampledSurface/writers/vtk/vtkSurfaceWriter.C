/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "vtkSurfaceWriter.H"

#include "OFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::vtkSurfaceWriter<Type>::writeGeometry
(
    Ostream& os,
    const pointField& points,
    const faceList& faces
)
{
    // header
    os
        << "# vtk DataFile Version 2.0" << nl
        << "sampleSurface" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl;

    // Write vertex coords
    os  << "POINTS " << points.size() << " float" << nl;
    forAll(points, pointI)
    {
        const point& pt = points[pointI];
        os  << float(pt.x()) << ' '
            << float(pt.y()) << ' '
            << float(pt.z()) << nl;
    }
    os  << nl;


    // Write faces
    label nNodes = 0;
    forAll(faces, faceI)
    {
        nNodes += faces[faceI].size();
    }

    os  << "POLYGONS " << faces.size() << ' '
        << faces.size() + nNodes << nl;

    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        os << f.size();
        forAll(f, fp)
        {
            os << ' ' << f[fp];
        }
        os << nl;
    }
}


namespace Foam
{

    // Write scalarField in vtk format
    template<>
    void Foam::vtkSurfaceWriter<Foam::scalar>::writeData
    (
        Ostream& os,
        const Field<Foam::scalar>& values
    )
    {
        os << "1 " << values.size() << " float" << nl;

        forAll(values, elemI)
        {
            if (elemI)
            {
                if (elemI % 10)
                {
                    os << ' ';
                }
                else
                {
                    os << nl;
                }
            }

            const scalar& v = values[elemI];
            os << float(v);
        }
        os << nl;
    }

    // Write vectorField in vtk format
    template<>
    void Foam::vtkSurfaceWriter<Foam::vector>::writeData
    (
        Ostream& os,
        const Field<Foam::vector>& values
    )
    {
        os << "3 " << values.size() << " float" << nl;

        forAll(values, elemI)
        {
            const vector& v = values[elemI];
            os  << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
                << nl;
        }
    }


    // Write sphericalTensorField in vtk format
    template<>
    void Foam::vtkSurfaceWriter<Foam::sphericalTensor>::writeData
    (
        Ostream& os,
        const Field<sphericalTensor>& values
    )
    {
        os << "1 " << values.size() << " float" << nl;

        forAll(values, elemI)
        {
            const sphericalTensor& v = values[elemI];
            os << float(v[0]) << nl;
        }
    }


    // Write symmTensorField in vtk format
    template<>
    void Foam::vtkSurfaceWriter<Foam::symmTensor>::writeData
    (
        Ostream& os,
        const Field<symmTensor>& values
    )
    {
        os << "6 " << values.size() << " float" << nl;

        forAll(values, elemI)
        {
            const symmTensor& v = values[elemI];
            os  << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
                << float(v[3]) << ' ' << float(v[4]) << ' ' << float(v[5])
                << nl;

        }
    }


    // Write tensorField in vtk format
    template<>
    void Foam::vtkSurfaceWriter<Foam::tensor>::writeData
    (
        Ostream& os,
        const Field<tensor>& values
    )
    {
        os << "9 " << values.size() << " float" << nl;

        forAll(values, elemI)
        {
            const tensor& v = values[elemI];
            os  << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
                << float(v[3]) << ' ' << float(v[4]) << ' ' << float(v[5])
                << float(v[6]) << ' ' << float(v[7]) << ' ' << float(v[8])
                << nl;
        }
    }

}


// Write generic field in vtk format
template<class Type>
void Foam::vtkSurfaceWriter<Type>::writeData
(
    Ostream& os,
    const Field<Type>& values
)
{
    os << "1 " << values.size() << " float" << nl;

    forAll(values, elemI)
    {
        os << float(0) << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Type>
Foam::vtkSurfaceWriter<Type>::vtkSurfaceWriter()
:
    surfaceWriter<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::vtkSurfaceWriter<Type>::~vtkSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkSurfaceWriter<Type>::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const bool verbose
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    fileName fName(outputDir/surfaceName + ".vtk");

    if (verbose)
    {
        Info<< "Writing geometry to " << fName << endl;
    }

    OFstream os(fName);
    writeGeometry(os, points, faces);
}


template<class Type>
void Foam::vtkSurfaceWriter<Type>::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const fileName& fieldName,
    const Field<Type>& values,
    const bool verbose
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os
    (
        outputDir/fieldName + '_' + surfaceName + ".vtk"
    );

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << os.name() << endl;
    }

    writeGeometry(os, points, faces);

    // start writing data
    if (values.size() == points.size())
    {
        os  << "POINT_DATA ";
    }
    else
    {
        os  << "CELL_DATA ";
    }

    os  << values.size() << nl
        << "FIELD attributes 1" << nl
        << fieldName << " ";

    // Write data
    writeData(os, values);

}


// ************************************************************************* //
