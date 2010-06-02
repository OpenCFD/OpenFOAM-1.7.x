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

InClass
    vtkPV3Foam

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3FoamPoints_H
#define vtkPV3FoamPoints_H

// VTK includes
#include "vtkPoints.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline void vtkPV3FoamInsertNextPoint
(
    vtkPoints *points,
    const Foam::point& p
)
{
    points->InsertNextPoint(p.x(), p.y(), p.z());
}

#if 0
// this should be faster, but didn't get it working ...
inline void vtkPV3FoamSetPoint
(
    vtkPoints *points,
    const Foam::label id,
    const Foam::point& p
)
{
    points->SetPoint(id, p.x(), p.y(), p.z());
}


// Convert Foam mesh vertices to VTK
inline vtkPoints* vtkPV3FoamVTKPoints(const Foam::pointField& points)
{
    vtkPoints *vtkpoints = vtkPoints::New();
    vtkpoints->SetNumberOfPoints(points.size());
    forAll(points, i)
    {
        const Foam::point& p = points[i];
        vtkpoints->SetPoint(i, p.x(), p.y(), p.z());
    }

    return vtkpoints;
}

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
