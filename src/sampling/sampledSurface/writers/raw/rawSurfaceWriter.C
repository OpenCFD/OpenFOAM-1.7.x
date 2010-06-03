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

\*---------------------------------------------------------------------------*/

#include "rawSurfaceWriter.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::rawSurfaceWriter<Type>::writeGeometry
(
    const pointField& points,
    const label pointI,
    Ostream& os
)
{
    const point& pt = points[pointI];

    os << pt.x() << ' ' << pt.y() << ' ' << pt.z() << ' ';
}


template<class Type>
void Foam::rawSurfaceWriter<Type>::writeGeometry
(
    const pointField& points,
    const faceList& faces,
    const label faceI,
    Ostream& os
)
{
    const point& ct = faces[faceI].centre(points);

    os << ct.x() << ' ' << ct.y() << ' ' << ct.z() << ' ';
}


// Write scalarField in raw format
template<class Type>
void Foam::rawSurfaceWriter<Type>::writeData
(
    const fileName& fieldName,
    const pointField& points,
    const faceList& faces,
    const scalarField& values,
    Ostream& os
)
{
    // header
    os  << "#  x  y  z  " << fieldName << endl;

    // Write data
    if (values.size() == points.size())
    {
        forAll(values, elemI)
        {
            writeGeometry(points, elemI, os);
            os << values[elemI] << nl;
        }
    }
    else
    {
        forAll(values, elemI)
        {
            writeGeometry(points, faces, elemI, os);
            os << values[elemI] << nl;
        }
    }

    os << nl;
}


// Write vectorField in raw format
template<class Type>
void Foam::rawSurfaceWriter<Type>::writeData
(
    const fileName& fieldName,
    const pointField& points,
    const faceList& faces,
    const vectorField& values,
    Ostream& os
)
{
    // header
    os  << "#  x  y  z  "
        << fieldName << "_x  "
        << fieldName << "_y  "
        << fieldName << "_z  "
        << endl;

    // Write data
    if (values.size() == points.size())
    {
        forAll(values, elemI)
        {
            writeGeometry(points, elemI, os);

            const vector& v = values[elemI];
            os << v[0] << ' ' << v[1] << ' ' << v[2] << nl;
        }
    }
    else
    {
        forAll(values, elemI)
        {
            writeGeometry(points, faces, elemI, os);

            const vector& v = values[elemI];
            os << v[0] << ' ' << v[1] << ' ' << v[2] << nl;
        }
    }

}


// Write sphericalTensorField in raw format
template<class Type>
void Foam::rawSurfaceWriter<Type>::writeData
(
    const fileName& fieldName,
    const pointField& points,
    const faceList& faces,
    const sphericalTensorField& values,
    Ostream& os
)
{
    // header
    os  << "#  ii  ";
    os << fieldName << "_ii" << endl;

    // Write data
    if (values.size() == points.size())
    {
        forAll(values, elemI)
        {
            writeGeometry(points, elemI, os);

            const sphericalTensor& v = values[elemI];
            os  << v[0] << nl;
        }
    }
    else
    {
        forAll(values, elemI)
        {
            writeGeometry(points, faces, elemI, os);

            const sphericalTensor& v = values[elemI];
            os  << v[0] << nl;
        }
    }
}


// Write symmTensorField in raw format
template<class Type>
void Foam::rawSurfaceWriter<Type>::writeData
(
    const fileName& fieldName,
    const pointField& points,
    const faceList& faces,
    const symmTensorField& values,
    Ostream& os
)
{
    // header
    os  << "#  xx  xy  xz  yy  yz ";
    for(int i=0; i<6; i++)
    {
        os << fieldName << "_" << i << "  ";
    }
    os << endl;

    // Write data
    if (values.size() == points.size())
    {
        forAll(values, elemI)
        {
            writeGeometry(points, elemI, os);

            const symmTensor& v = values[elemI];

            os  << v[0] << ' ' << v[1] << ' ' << v[2] << ' '
                << v[3] << ' ' << v[4] << ' ' << v[5] << ' '
                << nl;
        }
    }
    else
    {
        forAll(values, elemI)
        {
            writeGeometry(points, faces, elemI, os);

            const symmTensor& v = values[elemI];

            os  << v[0] << ' ' << v[1] << ' ' << v[2] << ' '
                << v[3] << ' ' << v[4] << ' ' << v[5] << ' '
                << nl;
        }
    }
}


// Write tensorField in raw format
template<class Type>
void Foam::rawSurfaceWriter<Type>::writeData
(
    const fileName& fieldName,
    const pointField& points,
    const faceList& faces,
    const tensorField& values,
    Ostream& os
)
{
    // header
    os  << "#  xx  xy  xz  yx  yy  yz  zx  zy  zz";
    for (int i=0; i<9; ++i)
    {
        os << fieldName << "_" << i << "  ";
    }
    os << endl;

    // Write data
    if (values.size() == points.size())
    {
        forAll(values, elemI)
        {
            writeGeometry(points, elemI, os);

            const tensor& v = values[elemI];
            os  << v[0] << ' ' << v[1] << ' ' << v[2] << ' '
                << v[3] << ' ' << v[4] << ' ' << v[5] << ' '
                << v[6] << ' ' << v[7] << ' ' << v[8] << nl;
        }
    }
    else
    {
        forAll(values, elemI)
        {
            writeGeometry(points, faces, elemI, os);

            const tensor& v = values[elemI];
            os  << v[0] << ' ' << v[1] << ' ' << v[2] << ' '
                << v[3] << ' ' << v[4] << ' ' << v[5] << ' '
                << v[6] << ' ' << v[7] << ' ' << v[8] << nl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::rawSurfaceWriter<Type>::rawSurfaceWriter()
:
    surfaceWriter<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::rawSurfaceWriter<Type>::~rawSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::rawSurfaceWriter<Type>::write
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

    OFstream os
    (
        outputDir/surfaceName + ".raw"
    );

    if (verbose)
    {
        Info<< "Writing geometry to " << os.name() << endl;
    }


    // header
    os  << "# geometry NO_DATA " << faces.size() << nl
        << "#  x  y  z" << endl;

    // Write faces
    forAll(faces, elemI)
    {
        writeGeometry(points, faces, elemI, os);
        os << nl;
    }

    os << nl;
}


namespace Foam
{
    // bool fields aren't supported
    template<>
    void Foam::rawSurfaceWriter<bool>::write
    (
        const fileName& outputDir,
        const fileName& surfaceName,
        const pointField& points,
        const faceList& faces,
        const fileName& fieldName,
        const Field<bool>& values,
        const bool verbose
    ) const
    {}
}


template<class Type>
void Foam::rawSurfaceWriter<Type>::write
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
        outputDir/fieldName + '_' + surfaceName + ".raw"
    );

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << os.name() << endl;
    }


    // header
    os  << "# " << fieldName;
    if (values.size() == points.size())
    {
        os  << "  POINT_DATA ";
    }
    else
    {
        os  << "  FACE_DATA ";
    }

    os  << values.size() << nl;

    writeData(fieldName, points, faces, values, os);
}


// ************************************************************************* //
