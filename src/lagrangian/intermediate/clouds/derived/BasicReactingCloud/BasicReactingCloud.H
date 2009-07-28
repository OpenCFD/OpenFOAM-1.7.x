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

Class
    Foam::BasicReactingCloud

Description
    Reacting cloud templated on the type of carrier phase thermodynamics

SourceFiles
    BasicReactingCloud.C

\*---------------------------------------------------------------------------*/

#ifndef BasicReactingCloud_H
#define BasicReactingCloud_H

#include "ReactingCloud.H"
#include "BasicReactingParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class ThermoType>
class BasicReactingCloud;

/*---------------------------------------------------------------------------*\
                        Class BasicReactingCloud Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class BasicReactingCloud
:
    public ReactingCloud<BasicReactingParcel<ThermoType> >
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BasicReactingCloud(const BasicReactingCloud&);

        //- Disallow default bitwise assignment
        void operator=(const BasicReactingCloud&);


public:

    //- Runtime type information
    TypeName("BasicReactingCloud");


    // Constructors

        //- Construct given carrier gas fields
        BasicReactingCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            basicThermo& thermo
        );


    //- Destructor
    ~BasicReactingCloud();


    // Member Functions

        //- Write fields
        virtual void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BasicReactingCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
