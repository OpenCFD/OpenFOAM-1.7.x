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

#include "MeshedSurface.H"
#include "boundBox.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::writeStats(Ostream& os) const
{
    os  << "points      : " << this->points().size() << nl;
    if (MeshedSurface<Face>::isTri())
    {
        os << "triangles   : " << this->size() << nl;
    }
    else
    {
        label nTri = 0;
        label nQuad = 0;
        forAll(*this, i)
        {
            const label n = this->operator[](i).size();

            if (n == 3)
            {
                nTri++;
            }
            else if (n == 4)
            {
                nQuad++;
            }
        }

        os  << "faces       : " << this->size()
            << "  (tri:" << nTri << " quad:" << nQuad
            << " poly:" << (this->size() - nTri - nQuad ) << ")" << nl;
    }

    os  << "boundingBox : " << boundBox(this->points()) << endl;
}


// ************************************************************************* //
