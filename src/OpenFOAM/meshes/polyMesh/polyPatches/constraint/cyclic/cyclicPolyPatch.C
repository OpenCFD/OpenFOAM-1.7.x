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

#include "cyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "demandDrivenData.H"
#include "OFstream.H"
#include "patchZones.H"
#include "matchPoints.H"
#include "EdgeMap.H"
#include "Time.H"
#include "transformList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicPolyPatch, dictionary);


template<>
const char* NamedEnum<cyclicPolyPatch::transformType, 3>::names[] =
{
    "unknown",
    "rotational",
    "translational"
};

const NamedEnum<cyclicPolyPatch::transformType, 3>
    cyclicPolyPatch::transformTypeNames;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::cyclicPolyPatch::findMaxArea
(
    const pointField& points,
    const faceList& faces
)
{
    label maxI = -1;
    scalar maxAreaSqr = -GREAT;

    forAll(faces, faceI)
    {
        scalar areaSqr = magSqr(faces[faceI].normal(points));

        if (areaSqr > maxAreaSqr)
        {
            maxAreaSqr = areaSqr;
            maxI = faceI;
        }
    }
    return maxI;
}


void Foam::cyclicPolyPatch::calcTransforms()
{
    if (size())
    {
        const pointField& points = this->points();

        // Determine geometric quantities on the two halves
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        primitivePatch half0
        (
            SubList<face>
            (
                *this,
                size()/2
            ),
            points
        );

        pointField half0Ctrs(calcFaceCentres(half0, half0.points()));

        scalarField half0Tols(calcFaceTol(half0, half0.points(), half0Ctrs));

        primitivePatch half1
        (
            SubList<face>
            (
                *this,
                size()/2,
                size()/2
            ),
            points
        );
        pointField half1Ctrs(calcFaceCentres(half1, half1.points()));

        // Dump halves
        if (debug)
        {
            fileName casePath(boundaryMesh().mesh().time().path());

            fileName nm0(casePath/name()+"_half0_faces.obj");
            Pout<< "cyclicPolyPatch::calcTransforms : Writing half0"
                << " faces to OBJ file " << nm0 << endl;
            writeOBJ(nm0, half0, half0.points());

            fileName nm1(casePath/name()+"_half1_faces.obj");
            Pout<< "cyclicPolyPatch::calcTransforms : Writing half1"
                << " faces to OBJ file " << nm1 << endl;
            writeOBJ(nm1, half1, half1.points());

            OFstream str(casePath/name()+"_half0_to_half1.obj");
            label vertI = 0;
            Pout<< "cyclicPolyPatch::calcTransforms :"
                << " Writing coupled face centres as lines to " << str.name()
                << endl;
            forAll(half0Ctrs, i)
            {
                const point& p0 = half0Ctrs[i];
                str << "v " << p0.x() << ' ' << p0.y() << ' ' << p0.z() << nl;
                vertI++;
                const point& p1 = half1Ctrs[i];
                str << "v " << p1.x() << ' ' << p1.y() << ' ' << p1.z() << nl;
                vertI++;
                str << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }

        vectorField half0Normals(half0.size());
        vectorField half1Normals(half1.size());

        for (label facei = 0; facei < size()/2; facei++)
        {
            half0Normals[facei] = operator[](facei).normal(points);
            label nbrFacei = facei+size()/2;
            half1Normals[facei] = operator[](nbrFacei).normal(points);

            scalar magSf = mag(half0Normals[facei]);
            scalar nbrMagSf = mag(half1Normals[facei]);
            scalar avSf = (magSf + nbrMagSf)/2.0;

            if (magSf < ROOTVSMALL && nbrMagSf < ROOTVSMALL)
            {
                // Undetermined normal. Use dummy normal to force separation
                // check. (note use of sqrt(VSMALL) since that is how mag
                // scales)
                half0Normals[facei] = point(1, 0, 0);
                half1Normals[facei] = half0Normals[facei];
            }
            else if (mag(magSf - nbrMagSf)/avSf > coupledPolyPatch::matchTol)
            {
                FatalErrorIn
                (
                    "cyclicPolyPatch::calcTransforms()"
                )   << "face " << facei << " area does not match neighbour "
                    << nbrFacei << " by "
                    << 100*mag(magSf - nbrMagSf)/avSf
                    << "% -- possible face ordering problem." << endl
                    << "patch:" << name()
                    << " my area:" << magSf
                    << " neighbour area:" << nbrMagSf
                    << " matching tolerance:" << coupledPolyPatch::matchTol
                     << endl
                    << "Mesh face:" << start()+facei
                    << " vertices:"
                    << UIndirectList<point>(points, operator[](facei))()
                    << endl
                    << "Neighbour face:" << start()+nbrFacei
                    << " vertices:"
                    << UIndirectList<point>(points, operator[](nbrFacei))()
                    << endl
                    << "Rerun with cyclic debug flag set"
                    << " for more information." << exit(FatalError);
            }
            else
            {
                half0Normals[facei] /= magSf;
                half1Normals[facei] /= nbrMagSf;
            }
        }


        // See if transformation is prescribed
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        switch (transform_)
        {
            case ROTATIONAL:
            {
                // Specified single rotation tensor.

                // Get best fitting face and its opposite number
                label face0 = getConsistentRotationFace(half0Ctrs);
                label face1 = face0;

                vector n0 =
                    (
                        (half0Ctrs[face0] - rotationCentre_)
                      ^ rotationAxis_
                    );
                vector n1 =
                    (
                        (half1Ctrs[face1] - rotationCentre_)
                      ^ -rotationAxis_
                    );
                n0 /= mag(n0) + VSMALL;
                n1 /= mag(n1) + VSMALL;

                if (debug)
                {
                    Pout<< "cyclicPolyPatch::calcTransforms :"
                        << " Specified rotation :"
                        << " n0:" << n0 << " n1:" << n1 << endl;
                }

                // Calculate transformation tensors from face0,1 only.
                // Note: can use tight tolerance now.
                calcTransformTensors
                (
                    pointField(1, half0Ctrs[face0]),
                    pointField(1, half1Ctrs[face1]),
                    vectorField(1, n0),
                    vectorField(1, n1),
                    scalarField(1, half0Tols[face0]),
                    1E-4
                );

                break;
            }

            default:
            {
                // Calculate transformation tensors from all faces.
                calcTransformTensors
                (
                    half0Ctrs,
                    half1Ctrs,
                    half0Normals,
                    half1Normals,
                    half0Tols
                );

                break;
            }
        }
    }
}


// Get geometric zones of patch by looking at normals.
// Method 1: any edge with sharpish angle is edge between two halves.
//           (this will handle e.g. wedge geometries).
//           Also two fully disconnected regions will be handled this way.
// Method 2: sort faces into two halves based on face normal.
bool Foam::cyclicPolyPatch::getGeometricHalves
(
    const primitivePatch& pp,
    labelList& half0ToPatch,
    labelList& half1ToPatch
) const
{
    // Calculate normals
    const vectorField& faceNormals = pp.faceNormals();

    // Find edges with sharp angles.
    boolList regionEdge(pp.nEdges(), false);

    const labelListList& edgeFaces = pp.edgeFaces();

    label nRegionEdges = 0;

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        // Check manifold edges for sharp angle.
        // (Non-manifold already handled by patchZones)
        if (eFaces.size() == 2)
        {
            if ((faceNormals[eFaces[0]] & faceNormals[eFaces[1]])< featureCos_)
            {
                regionEdge[edgeI] = true;

                nRegionEdges++;
            }
        }
    }


    // For every face determine zone it is connected to (without crossing
    // any regionEdge)
    patchZones ppZones(pp, regionEdge);

    if (debug)
    {
        Pout<< "cyclicPolyPatch::getGeometricHalves : "
            << "Found " << nRegionEdges << " edges on patch " << name()
            << " where the cos of the angle between two connected faces"
            << " was less than " << featureCos_ << nl
            << "Patch divided by these and by single sides edges into "
            << ppZones.nZones() << " parts." << endl;


        // Dumping zones to obj files.

        labelList nZoneFaces(ppZones.nZones());

        for (label zoneI = 0; zoneI < ppZones.nZones(); zoneI++)
        {
            OFstream stream
            (
                boundaryMesh().mesh().time().path()
               /name()+"_zone_"+Foam::name(zoneI)+".obj"
            );
            Pout<< "cyclicPolyPatch::getGeometricHalves : Writing zone "
                << zoneI << " face centres to OBJ file " << stream.name()
                << endl;

            labelList zoneFaces(findIndices(ppZones, zoneI));

            forAll(zoneFaces, i)
            {
                writeOBJ(stream, pp[zoneFaces[i]].centre(pp.points()));
            }

            nZoneFaces[zoneI] = zoneFaces.size();
        }
    }


    if (ppZones.nZones() == 2)
    {
        half0ToPatch = findIndices(ppZones, 0);
        half1ToPatch = findIndices(ppZones, 1);
    }
    else
    {
        if (debug)
        {
            Pout<< "cyclicPolyPatch::getGeometricHalves :"
                << " falling back to face-normal comparison" << endl;
        }
        label n0Faces = 0;
        half0ToPatch.setSize(pp.size());

        label n1Faces = 0;
        half1ToPatch.setSize(pp.size());

        // Compare to face 0 normal.
        forAll(faceNormals, faceI)
        {
            if ((faceNormals[faceI] & faceNormals[0]) > 0)
            {
                half0ToPatch[n0Faces++] = faceI;
            }
            else
            {
                half1ToPatch[n1Faces++] = faceI;
            }
        }
        half0ToPatch.setSize(n0Faces);
        half1ToPatch.setSize(n1Faces);

        if (debug)
        {
            Pout<< "cyclicPolyPatch::getGeometricHalves :"
                << " Number of faces per zone:("
                << n0Faces << ' ' << n1Faces << ')' << endl;
        }
    }

    if (half0ToPatch.size() != half1ToPatch.size())
    {
        fileName casePath(boundaryMesh().mesh().time().path());

        // Dump halves
        {
            fileName nm0(casePath/name()+"_half0_faces.obj");
            Pout<< "cyclicPolyPatch::getGeometricHalves : Writing half0"
                << " faces to OBJ file " << nm0 << endl;
            writeOBJ(nm0, UIndirectList<face>(pp, half0ToPatch)(), pp.points());

            fileName nm1(casePath/name()+"_half1_faces.obj");
            Pout<< "cyclicPolyPatch::getGeometricHalves : Writing half1"
                << " faces to OBJ file " << nm1 << endl;
            writeOBJ(nm1, UIndirectList<face>(pp, half1ToPatch)(), pp.points());
        }

        // Dump face centres
        {
            OFstream str0(casePath/name()+"_half0.obj");
            Pout<< "cyclicPolyPatch::getGeometricHalves : Writing half0"
                << " face centres to OBJ file " << str0.name() << endl;

            forAll(half0ToPatch, i)
            {
                writeOBJ(str0, pp[half0ToPatch[i]].centre(pp.points()));
            }

            OFstream str1(casePath/name()+"_half1.obj");
            Pout<< "cyclicPolyPatch::getGeometricHalves : Writing half1"
                << " face centres to OBJ file " << str1.name() << endl;
            forAll(half1ToPatch, i)
            {
                writeOBJ(str1, pp[half1ToPatch[i]].centre(pp.points()));
            }
        }

        SeriousErrorIn
        (
            "cyclicPolyPatch::getGeometricHalves"
            "(const primitivePatch&, labelList&, labelList&) const"
        )   << "Patch " << name() << " gets decomposed in two zones of"
            << "inequal size: " << half0ToPatch.size()
            << " and " << half1ToPatch.size() << endl
            << "This means that the patch is either not two separate regions"
            << " or one region where the angle between the different regions"
            << " is not sufficiently sharp." << endl
            << "Please adapt the featureCos setting." << endl
            << "Continuing with incorrect face ordering from now on!" << endl;

        return false;
    }
    else
    {
        return true;
    }
}


// Given a split of faces into left and right half calculate the centres
// and anchor points. Transform the left points so they align with the
// right ones.
void Foam::cyclicPolyPatch::getCentresAndAnchors
(
    const primitivePatch& pp,
    const faceList& half0Faces,
    const faceList& half1Faces,

    pointField& ppPoints,
    pointField& half0Ctrs,
    pointField& half1Ctrs,
    pointField& anchors0,
    scalarField& tols
) const
{
    // Get geometric data on both halves.
    half0Ctrs = calcFaceCentres(half0Faces, pp.points());
    anchors0 = getAnchorPoints(half0Faces, pp.points());
    half1Ctrs = calcFaceCentres(half1Faces, pp.points());

    switch (transform_)
    {
        case ROTATIONAL:
        {
            label face0 = getConsistentRotationFace(half0Ctrs);
            label face1 = getConsistentRotationFace(half1Ctrs);

            vector n0 = ((half0Ctrs[face0] - rotationCentre_) ^ rotationAxis_);
            vector n1 = ((half1Ctrs[face1] - rotationCentre_) ^ -rotationAxis_);
            n0 /= mag(n0) + VSMALL;
            n1 /= mag(n1) + VSMALL;

            if (debug)
            {
                Pout<< "cyclicPolyPatch::getCentresAndAnchors :"
                    << " Specified rotation :"
                    << " n0:" << n0 << " n1:" << n1 << endl;
            }

            // Rotation (around origin)
            const tensor reverseT(rotationTensor(n0, -n1));

            // Rotation
            forAll(half0Ctrs, faceI)
            {
                half0Ctrs[faceI] = Foam::transform(reverseT, half0Ctrs[faceI]);
                anchors0[faceI] = Foam::transform(reverseT, anchors0[faceI]);
            }

            ppPoints = Foam::transform(reverseT, pp.points());

            break;
        }
        //- Problem: usually specified translation is not accurate enough
        //- to get proper match so keep automatic determination over here.
        //case TRANSLATIONAL:
        //{
        //    // Transform 0 points.
        //
        //    if (debug)
        //    {
        //        Pout<< "cyclicPolyPatch::getCentresAndAnchors :"
        //            << "Specified translation : " << separationVector_ << endl;
        //    }
        //
        //    half0Ctrs += separationVector_;
        //    anchors0 += separationVector_;
        //    break;
        //}
        default:
        {
            // Assumes that cyclic is planar. This is also the initial
            // condition for patches without faces.

            // Determine the face with max area on both halves. These
            // two faces are used to determine the transformation tensors
            label max0I = findMaxArea(pp.points(), half0Faces);
            vector n0 = half0Faces[max0I].normal(pp.points());
            n0 /= mag(n0) + VSMALL;

            label max1I = findMaxArea(pp.points(), half1Faces);
            vector n1 = half1Faces[max1I].normal(pp.points());
            n1 /= mag(n1) + VSMALL;

            if (mag(n0 & n1) < 1-coupledPolyPatch::matchTol)
            {
                if (debug)
                {
                    Pout<< "cyclicPolyPatch::getCentresAndAnchors :"
                        << " Detected rotation :"
                        << " n0:" << n0 << " n1:" << n1 << endl;
                }

                // Rotation (around origin)
                const tensor reverseT(rotationTensor(n0, -n1));

                // Rotation
                forAll(half0Ctrs, faceI)
                {
                    half0Ctrs[faceI] = Foam::transform
                    (
                        reverseT,
                        half0Ctrs[faceI]
                    );
                    anchors0[faceI] = Foam::transform
                    (
                        reverseT,
                        anchors0[faceI]
                    );
                }
                ppPoints = Foam::transform(reverseT, pp.points());
            }
            else
            {
                // Parallel translation. Get average of all used points.

                primitiveFacePatch half0(half0Faces, pp.points());
                const pointField& half0Pts = half0.localPoints();
                const point ctr0(sum(half0Pts)/half0Pts.size());

                primitiveFacePatch half1(half1Faces, pp.points());
                const pointField& half1Pts = half1.localPoints();
                const point ctr1(sum(half1Pts)/half1Pts.size());

                if (debug)
                {
                    Pout<< "cyclicPolyPatch::getCentresAndAnchors :"
                        << " Detected translation :"
                        << " n0:" << n0 << " n1:" << n1
                        << " ctr0:" << ctr0 << " ctr1:" << ctr1 << endl;
                }

                half0Ctrs += ctr1 - ctr0;
                anchors0 += ctr1 - ctr0;
                ppPoints = pp.points() + ctr1 - ctr0;
            }
            break;
        }
    }


    // Calculate typical distance per face
    tols = calcFaceTol(half1Faces, pp.points(), half1Ctrs);
}


// Calculates faceMap and rotation. Returns true if all ok.
bool Foam::cyclicPolyPatch::matchAnchors
(
    const bool report,
    const primitivePatch& pp,
    const labelList& half0ToPatch,
    const pointField& anchors0,

    const labelList& half1ToPatch,
    const faceList& half1Faces,
    const labelList& from1To0,

    const scalarField& tols,

    labelList& faceMap,
    labelList& rotation
) const
{
    // Set faceMap such that half0 faces get first and corresponding half1
    // faces last.

    forAll(half0ToPatch, half0FaceI)
    {
        // Label in original patch
        label patchFaceI = half0ToPatch[half0FaceI];

        faceMap[patchFaceI] = half0FaceI;

        // No rotation
        rotation[patchFaceI] = 0;
    }

    bool fullMatch = true;

    forAll(from1To0, half1FaceI)
    {
        label patchFaceI = half1ToPatch[half1FaceI];

        // This face has to match the corresponding one on half0.
        label half0FaceI = from1To0[half1FaceI];

        label newFaceI = half0FaceI + pp.size()/2;

        faceMap[patchFaceI] = newFaceI;

        // Rotate patchFaceI such that its f[0] aligns with that of
        // the corresponding face
        // (which after shuffling will be at position half0FaceI)

        const point& wantedAnchor = anchors0[half0FaceI];

        rotation[newFaceI] = getRotation
        (
            pp.points(),
            half1Faces[half1FaceI],
            wantedAnchor,
            tols[half1FaceI]
        );

        if (rotation[newFaceI] == -1)
        {
            fullMatch = false;

            if (report)
            {
                const face& f = half1Faces[half1FaceI];
                SeriousErrorIn
                (
                    "cyclicPolyPatch::matchAnchors(..)"
                )   << "Patch:" << name() << " : "
                    << "Cannot find point on face " << f
                    << " with vertices:"
                    << UIndirectList<point>(pp.points(), f)()
                    << " that matches point " << wantedAnchor
                    << " when matching the halves of cyclic patch " << name()
                    << endl
                    << "Continuing with incorrect face ordering from now on!"
                    << endl;
            }
        }
    }
    return fullMatch;
}


Foam::label Foam::cyclicPolyPatch::getConsistentRotationFace
(
    const pointField& faceCentres
) const
{
    const scalarField magRadSqr =
        magSqr((faceCentres - rotationCentre_) ^ rotationAxis_);
    scalarField axisLen = (faceCentres - rotationCentre_) & rotationAxis_;
    axisLen = axisLen - min(axisLen);
    const scalarField magLenSqr = magRadSqr + axisLen*axisLen;

    label rotFace = -1;
    scalar maxMagLenSqr = -GREAT;
    scalar maxMagRadSqr = -GREAT;
    forAll(faceCentres, i)
    {
        if (magLenSqr[i] >= maxMagLenSqr)
        {
            if (magRadSqr[i] > maxMagRadSqr)
            {
                rotFace = i;
                maxMagLenSqr = magLenSqr[i];
                maxMagRadSqr = magRadSqr[i];
            }
        }
    }

    if (debug)
    {
        Info<< "getConsistentRotationFace(const pointField&)" << nl
            << "    rotFace = " << rotFace << nl
            << "    point =  " << faceCentres[rotFace] << endl;
    }

    return rotFace;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, size, start, index, bm),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL),
    featureCos_(0.9),
    transform_(UNKNOWN),
    rotationAxis_(vector::zero),
    rotationCentre_(point::zero),
    separationVector_(vector::zero)
{}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL),
    featureCos_(0.9),
    transform_(UNKNOWN),
    rotationAxis_(vector::zero),
    rotationCentre_(point::zero),
    separationVector_(vector::zero)
{
    if (dict.found("neighbourPatch"))
    {
        FatalIOErrorIn
        (
            "cyclicPolyPatch::cyclicPolyPatch\n"
            "(\n"
            "    const word& name,\n"
            "    const dictionary& dict,\n"
            "    const label index,\n"
            "    const polyBoundaryMesh& bm\n"
            ")",
            dict
        )   << "Found \"neighbourPatch\" entry when reading cyclic patch "
            << name << endl
            << "Is this mesh already with split cyclics?" << endl
            << "If so run a newer version that supports it"
            << ", if not comment out the \"neighbourPatch\" entry and re-run"
            << exit(FatalIOError);
    }

    dict.readIfPresent("featureCos", featureCos_);

    if (dict.found("transform"))
    {
        transform_ = transformTypeNames.read(dict.lookup("transform"));
        switch (transform_)
        {
            case ROTATIONAL:
            {
                dict.lookup("rotationAxis") >> rotationAxis_;
                dict.lookup("rotationCentre") >> rotationCentre_;
                break;
            }
            case TRANSLATIONAL:
            {
                dict.lookup("separationVector") >> separationVector_;
                break;
            }
            default:
            {
                // no additional info required
            }
        }
    }
}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const cyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL),
    featureCos_(pp.featureCos_),
    transform_(pp.transform_),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    separationVector_(pp.separationVector_)
{}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const cyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL),
    featureCos_(pp.featureCos_),
    transform_(pp.transform_),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    separationVector_(pp.separationVector_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicPolyPatch::~cyclicPolyPatch()
{
    deleteDemandDrivenData(coupledPointsPtr_);
    deleteDemandDrivenData(coupledEdgesPtr_);
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cyclicPolyPatch::initGeometry()
{
    polyPatch::initGeometry();
}

void Foam::cyclicPolyPatch::calcGeometry()
{
    polyPatch::calcGeometry();
    calcTransforms();
}

void Foam::cyclicPolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::initMovePoints(p);
}

void Foam::cyclicPolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);
    calcTransforms();
}

void Foam::cyclicPolyPatch::initUpdateMesh()
{
    polyPatch::initUpdateMesh();
}

void Foam::cyclicPolyPatch::updateMesh()
{
    polyPatch::updateMesh();
    deleteDemandDrivenData(coupledPointsPtr_);
    deleteDemandDrivenData(coupledEdgesPtr_);
}


const Foam::edgeList& Foam::cyclicPolyPatch::coupledPoints() const
{
    if (!coupledPointsPtr_)
    {
        // Look at cyclic patch as two halves, A and B.
        // Now all we know is that relative face index in halfA is same
        // as coupled face in halfB and also that the 0th vertex
        // corresponds.

        // From halfA point to halfB or -1.
        labelList coupledPoint(nPoints(), -1);

        for (label patchFaceA = 0; patchFaceA < size()/2; patchFaceA++)
        {
            const face& fA = localFaces()[patchFaceA];

            forAll(fA, indexA)
            {
                label patchPointA = fA[indexA];

                if (coupledPoint[patchPointA] == -1)
                {
                    const face& fB = localFaces()[patchFaceA + size()/2];

                    label indexB = (fB.size() - indexA) % fB.size();

                    // Filter out points on wedge axis
                    if (patchPointA != fB[indexB])
                    {
                        coupledPoint[patchPointA] = fB[indexB];
                    }
                }
            }
        }

        coupledPointsPtr_ = new edgeList(nPoints());
        edgeList& connected = *coupledPointsPtr_;

        // Extract coupled points.
        label connectedI = 0;

        forAll(coupledPoint, i)
        {
            if (coupledPoint[i] != -1)
            {
                connected[connectedI++] = edge(i, coupledPoint[i]);
            }
        }

        connected.setSize(connectedI);

        if (debug)
        {
            OFstream str
            (
                boundaryMesh().mesh().time().path()
               /"coupledPoints.obj"
            );
            label vertI = 0;

            Pout<< "Writing file " << str.name() << " with coordinates of "
                << "coupled points" << endl;

            forAll(connected, i)
            {
                const point& a = points()[meshPoints()[connected[i][0]]];
                const point& b = points()[meshPoints()[connected[i][1]]];

                str<< "v " << a.x() << ' ' << a.y() << ' ' << a.z() << nl;
                str<< "v " << b.x() << ' ' << b.y() << ' ' << b.z() << nl;
                vertI += 2;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }

        // Remove any addressing calculated for the coupled edges calculation
        const_cast<primitivePatch&>(static_cast<const primitivePatch&>(*this))
            .clearOut();
    }
    return *coupledPointsPtr_;
}


const Foam::edgeList& Foam::cyclicPolyPatch::coupledEdges() const
{
    if (!coupledEdgesPtr_)
    {
        // Build map from points on halfA to points on halfB.
        const edgeList& pointCouples = coupledPoints();

        Map<label> aToB(2*pointCouples.size());

        forAll(pointCouples, i)
        {
            const edge& e = pointCouples[i];

            aToB.insert(e[0], e[1]);
        }

        // Map from edge on half A to points (in halfB indices)
        EdgeMap<label> edgeMap(nEdges());

        for (label patchFaceA = 0; patchFaceA < size()/2; patchFaceA++)
        {
            const labelList& fEdges = faceEdges()[patchFaceA];

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                const edge& e = edges()[edgeI];

                // Convert edge end points to corresponding points on halfB
                // side.
                Map<label>::const_iterator fnd0 = aToB.find(e[0]);
                if (fnd0 != aToB.end())
                {
                    Map<label>::const_iterator fnd1 = aToB.find(e[1]);
                    if (fnd1 != aToB.end())
                    {
                        edgeMap.insert(edge(fnd0(), fnd1()), edgeI);
                    }
                }
            }
        }

        coupledEdgesPtr_ = new edgeList(nEdges()/2);
        edgeList& coupledEdges = *coupledEdgesPtr_;
        label coupleI = 0;

        for (label patchFaceB = size()/2; patchFaceB < size(); patchFaceB++)
        {
            const labelList& fEdges = faceEdges()[patchFaceB];

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                const edge& e = edges()[edgeI];

                // Look up halfA edge from HashTable.
                EdgeMap<label>::iterator iter = edgeMap.find(e);

                if (iter != edgeMap.end())
                {
                    label halfAEdgeI = iter();

                    // Store correspondence. Filter out edges on wedge axis.
                    if (halfAEdgeI != edgeI)
                    {
                        coupledEdges[coupleI++] = edge(halfAEdgeI, edgeI);
                    }

                    // Remove so we build unique list
                    edgeMap.erase(iter);
                }
            }
        }
        coupledEdges.setSize(coupleI);


        // Some checks

        forAll(coupledEdges, i)
        {
            const edge& e = coupledEdges[i];

            if (e[0] == e[1] || e[0] < 0 || e[1] < 0)
            {
                FatalErrorIn("cyclicPolyPatch::coupledEdges() const")
                    << "Problem : at position " << i
                    << " illegal couple:" << e
                    << abort(FatalError);
            }
        }

        if (debug)
        {
            OFstream str
            (
                boundaryMesh().mesh().time().path()
               /"coupledEdges.obj"
            );
            label vertI = 0;

            Pout<< "Writing file " << str.name() << " with centres of "
                << "coupled edges" << endl;

            forAll(coupledEdges, i)
            {
                const edge& e = coupledEdges[i];

                const point& a = edges()[e[0]].centre(localPoints());
                const point& b = edges()[e[1]].centre(localPoints());

                str<< "v " << a.x() << ' ' << a.y() << ' ' << a.z() << nl;
                str<< "v " << b.x() << ' ' << b.y() << ' ' << b.z() << nl;
                vertI += 2;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }

        // Remove any addressing calculated for the coupled edges calculation
        const_cast<primitivePatch&>(static_cast<const primitivePatch&>(*this))
            .clearOut();
    }
    return *coupledEdgesPtr_;
}


void Foam::cyclicPolyPatch::initOrder(const primitivePatch& pp) const
{}


//  Return new ordering. Ordering is -faceMap: for every face index
//  the new face -rotation:for every new face the clockwise shift
//  of the original face. Return false if nothing changes (faceMap
//  is identity, rotation is 0)
bool Foam::cyclicPolyPatch::order
(
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    if (pp.empty())
    {
        // No faces, nothing to change.
        return false;
    }

    if (pp.size()&1)
    {
        FatalErrorIn("cyclicPolyPatch::order(..)")
            << "Size of cyclic " << name() << " should be a multiple of 2"
            << ". It is " << pp.size() << abort(FatalError);
    }

    label halfSize = pp.size()/2;

    // Supplied primitivePatch already with new points.
    // Cyclics are limited to one transformation tensor
    // currently anyway (i.e. straight plane) so should not be too big a
    // problem.


    // Indices of faces on half0
    labelList half0ToPatch;
    // Indices of faces on half1
    labelList half1ToPatch;


    // 1. Test if already correctly oriented by starting from trivial ordering.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    half0ToPatch = identity(halfSize);
    half1ToPatch = half0ToPatch + halfSize;

    // Get faces
    faceList half0Faces(UIndirectList<face>(pp, half0ToPatch));
    faceList half1Faces(UIndirectList<face>(pp, half1ToPatch));

    // Get geometric quantities
    pointField half0Ctrs, half1Ctrs, anchors0, ppPoints;
    scalarField tols;
    getCentresAndAnchors
    (
        pp,
        half0Faces,
        half1Faces,

        ppPoints,
        half0Ctrs,
        half1Ctrs,
        anchors0,
        tols
    );

    // Geometric match of face centre vectors
    labelList from1To0;
    bool matchedAll = matchPoints
    (
        half1Ctrs,
        half0Ctrs,
        tols,
        false,
        from1To0
    );

    if (debug)
    {
        Pout<< "cyclicPolyPatch::order : test if already ordered:"
            << matchedAll << endl;

        // Dump halves
        fileName nm0("match1_"+name()+"_half0_faces.obj");
        Pout<< "cyclicPolyPatch::order : Writing half0"
            << " faces to OBJ file " << nm0 << endl;
        writeOBJ(nm0, half0Faces, ppPoints);

        fileName nm1("match1_"+name()+"_half1_faces.obj");
        Pout<< "cyclicPolyPatch::order : Writing half1"
            << " faces to OBJ file " << nm1 << endl;
        writeOBJ(nm1, half1Faces, pp.points());

        OFstream ccStr
        (
            boundaryMesh().mesh().time().path()
           /"match1_"+ name() + "_faceCentres.obj"
        );
        Pout<< "cyclicPolyPatch::order : "
            << "Dumping currently found cyclic match as lines between"
            << " corresponding face centres to file " << ccStr.name()
            << endl;

        // Recalculate untransformed face centres
        //pointField rawHalf0Ctrs = calcFaceCentres(half0Faces, pp.points());
        label vertI = 0;

        forAll(half1Ctrs, i)
        {
            //if (from1To0[i] != -1)
            {
                // Write edge between c1 and c0
                //const point& c0 = rawHalf0Ctrs[from1To0[i]];
                //const point& c0 = half0Ctrs[from1To0[i]];
                const point& c0 = half0Ctrs[i];
                const point& c1 = half1Ctrs[i];
                writeOBJ(ccStr, c0, c1, vertI);
            }
        }
    }


    // 2. Ordered in pairs (so 0,1 coupled and 2,3 etc.)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (!matchedAll)
    {
        label faceI = 0;
        for (label i = 0; i < halfSize; i++)
        {
            half0ToPatch[i] = faceI++;
            half1ToPatch[i] = faceI++;
        }

        // And redo all matching
        half0Faces = UIndirectList<face>(pp, half0ToPatch);
        half1Faces = UIndirectList<face>(pp, half1ToPatch);

        getCentresAndAnchors
        (
            pp,
            half0Faces,
            half1Faces,

            ppPoints,
            half0Ctrs,
            half1Ctrs,
            anchors0,
            tols
        );

        // Geometric match of face centre vectors
        matchedAll = matchPoints
        (
            half1Ctrs,
            half0Ctrs,
            tols,
            false,
            from1To0
        );

        if (debug)
        {
            Pout<< "cyclicPolyPatch::order : test if pairwise ordered:"
                << matchedAll << endl;
            // Dump halves
            fileName nm0("match2_"+name()+"_half0_faces.obj");
            Pout<< "cyclicPolyPatch::order : Writing half0"
                << " faces to OBJ file " << nm0 << endl;
            writeOBJ(nm0, half0Faces, ppPoints);

            fileName nm1("match2_"+name()+"_half1_faces.obj");
            Pout<< "cyclicPolyPatch::order : Writing half1"
                << " faces to OBJ file " << nm1 << endl;
            writeOBJ(nm1, half1Faces, pp.points());

            OFstream ccStr
            (
                boundaryMesh().mesh().time().path()
               /"match2_"+name()+"_faceCentres.obj"
            );
            Pout<< "cyclicPolyPatch::order : "
                << "Dumping currently found cyclic match as lines between"
                << " corresponding face centres to file " << ccStr.name()
                << endl;

            // Recalculate untransformed face centres
            //pointField rawHalf0Ctrs = calcFaceCentres(half0Faces, pp.points());
            label vertI = 0;

            forAll(half1Ctrs, i)
            {
                if (from1To0[i] != -1)
                {
                    // Write edge between c1 and c0
                    //const point& c0 = rawHalf0Ctrs[from1To0[i]];
                    const point& c0 = half0Ctrs[from1To0[i]];
                    const point& c1 = half1Ctrs[i];
                    writeOBJ(ccStr, c0, c1, vertI);
                }
            }
        }
    }


    // 3. Baffles(coincident faces) converted into cyclics (e.g. jump)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (!matchedAll)
    {
        label baffleI = 0;

        forAll(pp, faceI)
        {
            const face& f = pp.localFaces()[faceI];
            const labelList& pFaces = pp.pointFaces()[f[0]];

            label matchedFaceI = -1;

            forAll(pFaces, i)
            {
                label otherFaceI = pFaces[i];

                if (otherFaceI > faceI)
                {
                    const face& otherF = pp.localFaces()[otherFaceI];

                    // Note: might pick up two similar oriented faces
                    //       (but that is illegal anyway)
                    if (f == otherF)
                    {
                        matchedFaceI = otherFaceI;
                        break;
                    }
                }
            }

            if (matchedFaceI != -1)
            {
                half0ToPatch[baffleI] = faceI;
                half1ToPatch[baffleI] = matchedFaceI;
                baffleI++;
            }
        }

        if (baffleI == halfSize)
        {
            // And redo all matching
            half0Faces = UIndirectList<face>(pp, half0ToPatch);
            half1Faces = UIndirectList<face>(pp, half1ToPatch);

            getCentresAndAnchors
            (
                pp,
                half0Faces,
                half1Faces,

                ppPoints,
                half0Ctrs,
                half1Ctrs,
                anchors0,
                tols
            );

            // Geometric match of face centre vectors
            matchedAll = matchPoints
            (
                half1Ctrs,
                half0Ctrs,
                tols,
                false,
                from1To0
            );

            if (debug)
            {
                Pout<< "cyclicPolyPatch::order : test if baffles:"
                    << matchedAll << endl;
                // Dump halves
                fileName nm0("match3_"+name()+"_half0_faces.obj");
                Pout<< "cyclicPolyPatch::order : Writing half0"
                    << " faces to OBJ file " << nm0 << endl;
                writeOBJ(nm0, half0Faces, ppPoints);

                fileName nm1("match3_"+name()+"_half1_faces.obj");
                Pout<< "cyclicPolyPatch::order : Writing half1"
                    << " faces to OBJ file " << nm1 << endl;
                writeOBJ(nm1, half1Faces, pp.points());

                OFstream ccStr
                (
                    boundaryMesh().mesh().time().path()
                   /"match3_"+ name() + "_faceCentres.obj"
                );
                Pout<< "cyclicPolyPatch::order : "
                    << "Dumping currently found cyclic match as lines between"
                    << " corresponding face centres to file " << ccStr.name()
                    << endl;

                // Recalculate untransformed face centres
                //pointField rawHalf0Ctrs = calcFaceCentres(half0Faces, pp.points());
                label vertI = 0;

                forAll(half1Ctrs, i)
                {
                    if (from1To0[i] != -1)
                    {
                        // Write edge between c1 and c0
                        //const point& c0 = rawHalf0Ctrs[from1To0[i]];
                        const point& c0 = half0Ctrs[from1To0[i]];
                        const point& c1 = half1Ctrs[i];
                        writeOBJ(ccStr, c0, c1, vertI);
                    }
                }
            }
        }
    }


    // 4. Automatic geometric ordering
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (!matchedAll)
    {
        // Split faces according to feature angle or topology
        bool okSplit = getGeometricHalves(pp, half0ToPatch, half1ToPatch);

        if (!okSplit)
        {
            // Did not split into two equal parts.
            return false;
        }

        // And redo all matching
        half0Faces = UIndirectList<face>(pp, half0ToPatch);
        half1Faces = UIndirectList<face>(pp, half1ToPatch);

        getCentresAndAnchors
        (
            pp,
            half0Faces,
            half1Faces,

            ppPoints,
            half0Ctrs,
            half1Ctrs,
            anchors0,
            tols
        );

        // Geometric match of face centre vectors
        matchedAll = matchPoints
        (
            half1Ctrs,
            half0Ctrs,
            tols,
            false,
            from1To0
        );

        if (debug)
        {
            Pout<< "cyclicPolyPatch::order : automatic ordering result:"
                << matchedAll << endl;
            // Dump halves
            fileName nm0("match4_"+name()+"_half0_faces.obj");
            Pout<< "cyclicPolyPatch::order : Writing half0"
                << " faces to OBJ file " << nm0 << endl;
            writeOBJ(nm0, half0Faces, ppPoints);

            fileName nm1("match4_"+name()+"_half1_faces.obj");
            Pout<< "cyclicPolyPatch::order : Writing half1"
                << " faces to OBJ file " << nm1 << endl;
            writeOBJ(nm1, half1Faces, pp.points());

            OFstream ccStr
            (
                boundaryMesh().mesh().time().path()
               /"match4_"+ name() + "_faceCentres.obj"
            );
            Pout<< "cyclicPolyPatch::order : "
                << "Dumping currently found cyclic match as lines between"
                << " corresponding face centres to file " << ccStr.name()
                << endl;

            // Recalculate untransformed face centres
            //pointField rawHalf0Ctrs = calcFaceCentres(half0Faces, pp.points());
            label vertI = 0;

            forAll(half1Ctrs, i)
            {
                if (from1To0[i] != -1)
                {
                    // Write edge between c1 and c0
                    //const point& c0 = rawHalf0Ctrs[from1To0[i]];
                    const point& c0 = half0Ctrs[from1To0[i]];
                    const point& c1 = half1Ctrs[i];
                    writeOBJ(ccStr, c0, c1, vertI);
                }
            }
        }
    }


    if (!matchedAll || debug)
    {
        // Dump halves
        fileName nm0(name()+"_half0_faces.obj");
        Pout<< "cyclicPolyPatch::order : Writing half0"
            << " faces to OBJ file " << nm0 << endl;
        writeOBJ(nm0, half0Faces, pp.points());

        fileName nm1(name()+"_half1_faces.obj");
        Pout<< "cyclicPolyPatch::order : Writing half1"
            << " faces to OBJ file " << nm1 << endl;
        writeOBJ(nm1, half1Faces, pp.points());

        OFstream ccStr
        (
            boundaryMesh().mesh().time().path()
           /name() + "_faceCentres.obj"
        );
        Pout<< "cyclicPolyPatch::order : "
            << "Dumping currently found cyclic match as lines between"
            << " corresponding face centres to file " << ccStr.name()
            << endl;

        // Recalculate untransformed face centres
        //pointField rawHalf0Ctrs = calcFaceCentres(half0Faces, pp.points());
        label vertI = 0;

        forAll(half1Ctrs, i)
        {
            if (from1To0[i] != -1)
            {
                // Write edge between c1 and c0
                //const point& c0 = rawHalf0Ctrs[from1To0[i]];
                const point& c0 = half0Ctrs[from1To0[i]];
                const point& c1 = half1Ctrs[i];
                writeOBJ(ccStr, c0, c1, vertI);
            }
        }
    }


    if (!matchedAll)
    {
        SeriousErrorIn
        (
            "cyclicPolyPatch::order"
            "(const primitivePatch&, labelList&, labelList&) const"
        )   << "Patch:" << name() << " : "
            << "Cannot match vectors to faces on both sides of patch" << endl
            << "    Perhaps your faces do not match?"
            << " The obj files written contain the current match." << endl
            << "    Continuing with incorrect face ordering from now on!"
            << endl;

            return false;
    }


    // Set faceMap such that half0 faces get first and corresponding half1
    // faces last.
    matchAnchors
    (
        true,                   // report if anchor matching error
        pp,
        half0ToPatch,
        anchors0,
        half1ToPatch,
        half1Faces,
        from1To0,
        tols,
        faceMap,
        rotation
    );

    // Return false if no change neccesary, true otherwise.

    forAll(faceMap, faceI)
    {
        if (faceMap[faceI] != faceI || rotation[faceI] != 0)
        {
            return true;
        }
    }

    return false;
}


void Foam::cyclicPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    os.writeKeyword("featureCos") << featureCos_ << token::END_STATEMENT << nl;
    switch (transform_)
    {
        case ROTATIONAL:
        {
            os.writeKeyword("transform") << transformTypeNames[ROTATIONAL]
                << token::END_STATEMENT << nl;
            os.writeKeyword("rotationAxis") << rotationAxis_
                << token::END_STATEMENT << nl;
            os.writeKeyword("rotationCentre") << rotationCentre_
                << token::END_STATEMENT << nl;
            break;
        }
        case TRANSLATIONAL:
        {
            os.writeKeyword("transform") << transformTypeNames[TRANSLATIONAL]
                << token::END_STATEMENT << nl;
            os.writeKeyword("separationVector") << separationVector_
                << token::END_STATEMENT << nl;
            break;
        }
        default:
        {
            // no additional info to write
        }
    }
}


// ************************************************************************* //
