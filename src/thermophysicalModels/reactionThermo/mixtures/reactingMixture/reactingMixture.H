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
    Foam::reactingMixture

Description
    Foam::reactingMixture

SourceFiles
    reactingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef reactingMixture_H
#define reactingMixture_H

#include "chemistryReader.H"
#include "multiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class reactingMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class reactingMixture
:
    public autoPtr<chemistryReader<ThermoType> >,
    public multiComponentMixture<ThermoType>,
    public PtrList<Reaction<ThermoType> >
{
    // Private member functions

        //- Disallow default bitwise copy construct
        reactingMixture(const reactingMixture&);

        //- Disallow default bitwise assignment
        void operator=(const reactingMixture&);


public:

    //- The type of thermo package this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary and mesh
        reactingMixture(const dictionary&, const fvMesh&);


    //- Destructor
    virtual ~reactingMixture()
    {}


    // Member functions

        //- Read dictionary
        void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "reactingMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
