/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Application
    transformPoints

Description
    Transforms the mesh points in the polyMesh directory according to the
    translate, rotate and scale options.

Usage
    Options are:

    -translate vector
        Translates the points by the given vector,

    -rotate (vector vector)
        Rotates the points from the first vector to the second,

     or -yawPitchRoll (yawdegrees pitchdegrees rolldegrees)
     or -rollPitchYaw (rolldegrees pitchdegrees yawdegrees)

    -scale vector
        Scales the points by the given vector.

    The any or all of the three options may be specified and are processed
    in the above order.

    With -rotateFields (in combination with -rotate/yawPitchRoll/rollPitchYaw)
    it will also read & transform vector & tensor fields.

    Note:
    yaw (rotation about z)
    pitch (rotation about y)
    roll (rotation about x)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ReadFields.H"
#include "pointFields.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "IStringStream.H"

using namespace Foam;
using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void readAndRotateFields
(
    PtrList<GeoField>& flds,
    const fvMesh& mesh,
    const tensor& T,
    const IOobjectList& objects
)
{
    ReadFields(mesh, objects, flds);
    forAll(flds, i)
    {
        Info<< "Transforming " << flds[i].name() << endl;
        dimensionedTensor dimT("t", flds[i].dimensions(), T);
        transform(flds[i], dimT, flds[i]);
    }
}


void rotateFields(const argList& args, const Time& runTime, const tensor& T)
{
#   include "createNamedMesh.H"

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.

    PtrList<volScalarField> vsFlds;
    readAndRotateFields(vsFlds, mesh, T, objects);

    PtrList<volVectorField> vvFlds;
    readAndRotateFields(vvFlds, mesh, T, objects);

    PtrList<volSphericalTensorField> vstFlds;
    readAndRotateFields(vstFlds, mesh, T, objects);

    PtrList<volSymmTensorField> vsymtFlds;
    readAndRotateFields(vsymtFlds, mesh, T, objects);

    PtrList<volTensorField> vtFlds;
    readAndRotateFields(vtFlds, mesh, T, objects);

    // Read surface fields.

    PtrList<surfaceScalarField> ssFlds;
    readAndRotateFields(ssFlds, mesh, T, objects);

    PtrList<surfaceVectorField> svFlds;
    readAndRotateFields(svFlds, mesh, T, objects);

    PtrList<surfaceSphericalTensorField> sstFlds;
    readAndRotateFields(sstFlds, mesh, T, objects);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    readAndRotateFields(ssymtFlds, mesh, T, objects);

    PtrList<surfaceTensorField> stFlds;
    readAndRotateFields(stFlds, mesh, T, objects);

    mesh.write();
}


//  Main program:

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"
    argList::validOptions.insert("translate", "vector");
    argList::validOptions.insert("rotate", "(vector vector)");
    argList::validOptions.insert("rollPitchYaw", "(roll pitch yaw)");
    argList::validOptions.insert("yawPitchRoll", "(yaw pitch roll)");
    argList::validOptions.insert("rotateFields", "");
    argList::validOptions.insert("scale", "vector");

#   include "setRootCase.H"
#   include "createTime.H"

    word regionName = polyMesh::defaultRegion;
    fileName meshDir;

    if (args.optionReadIfPresent("region", regionName))
    {
        meshDir = regionName/polyMesh::meshSubDir;
    }
    else
    {
        meshDir = polyMesh::meshSubDir;
    }

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(meshDir, "points"),
            meshDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );


    if (args.options().empty())
    {
        FatalErrorIn(args.executable())
            << "No options supplied, please use one or more of "
               "-translate, -rotate or -scale options."
            << exit(FatalError);
    }

    if (args.optionFound("translate"))
    {
        vector transVector(args.optionLookup("translate")());

        Info<< "Translating points by " << transVector << endl;

        points += transVector;
    }

    if (args.optionFound("rotate"))
    {
        Pair<vector> n1n2(args.optionLookup("rotate")());
        n1n2[0] /= mag(n1n2[0]);
        n1n2[1] /= mag(n1n2[1]);
        tensor T = rotationTensor(n1n2[0], n1n2[1]);

        Info<< "Rotating points by " << T << endl;

        points = transform(T, points);

        if (args.optionFound("rotateFields"))
        {
            rotateFields(args, runTime, T);
        }
    }
    else if (args.optionFound("rollPitchYaw"))
    {
        vector v(args.optionLookup("rollPitchYaw")());

        Info<< "Rotating points by" << nl
            << "    roll  " << v.x() << nl
            << "    pitch " << v.y() << nl
            << "    yaw   " << v.z() << endl;


        // Convert to radians
        v *= pi/180.0;

        quaternion R(v.x(), v.y(), v.z());

        Info<< "Rotating points by quaternion " << R << endl;
        points = transform(R, points);

        if (args.optionFound("rotateFields"))
        {
            rotateFields(args, runTime, R.R());
        }
    }
    else if (args.optionFound("yawPitchRoll"))
    {
        vector v(args.optionLookup("yawPitchRoll")());

        Info<< "Rotating points by" << nl
            << "    yaw   " << v.x() << nl
            << "    pitch " << v.y() << nl
            << "    roll  " << v.z() << endl;


        // Convert to radians
        v *= pi/180.0;

        scalar yaw = v.x();
        scalar pitch = v.y();
        scalar roll = v.z();

        quaternion R = quaternion(vector(0, 0, 1), yaw);
        R *= quaternion(vector(0, 1, 0), pitch);
        R *= quaternion(vector(1, 0, 0), roll);

        Info<< "Rotating points by quaternion " << R << endl;
        points = transform(R, points);

        if (args.optionFound("rotateFields"))
        {
            rotateFields(args, runTime, R.R());
        }
    }

    if (args.optionFound("scale"))
    {
        vector scaleVector(args.optionLookup("scale")());

        Info<< "Scaling points by " << scaleVector << endl;

        points.replace(vector::X, scaleVector.x()*points.component(vector::X));
        points.replace(vector::Y, scaleVector.y()*points.component(vector::Y));
        points.replace(vector::Z, scaleVector.z()*points.component(vector::Z));
    }

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    Info << "Writing points into directory " << points.path() << nl << endl;
    points.write();

    return 0;
}


// ************************************************************************* //
