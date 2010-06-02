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

Class
    Foam::procLduInterface

Description
    Foam::procLduInterface

SourceFiles
    procLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef procLduInterface_H
#define procLduInterface_H

#include "labelList.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lduInterfaceField;

/*---------------------------------------------------------------------------*\
                           Class procLduInterface Declaration
\*---------------------------------------------------------------------------*/

class procLduInterface
{
    // Private data

        labelList faceCells_;
        scalarField coeffs_;
        label myProcNo_;
        label neighbProcNo_;


    // Private member functions

        //- Dissallow construction as copy
        procLduInterface(const procLduInterface&);


public:

    friend class LUscalarMatrix;


    // Constructors

        procLduInterface
        (
            const lduInterfaceField& interface,
            const scalarField& coeffs
        );

        procLduInterface(Istream& is);

        autoPtr<procLduInterface> clone()
        {
            notImplemented("procLduInterface::clone()");
            return autoPtr<procLduInterface>(NULL);
        }

        static autoPtr<procLduInterface> New(Istream& is)
        {
            return autoPtr<procLduInterface>(new procLduInterface(is));
        }
            

    // Ostream operator

        friend Ostream& operator<<(Ostream&, const procLduInterface&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
