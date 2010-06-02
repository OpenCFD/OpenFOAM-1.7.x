/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::staticPressure

Description
    Converts kinematic pressure to static pressure, from the name of the
    pressure field, and density, i.e.

        p_static = density*p_kinematic

SourceFiles
    staticPressure.C
    IOstaticPressure.H

\*---------------------------------------------------------------------------*/

#ifndef staticPressure_H
#define staticPressure_H

#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class staticPressure Declaration
\*---------------------------------------------------------------------------*/

class staticPressure
{
    // Private data

        //- Name of this set of staticPressure objects
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Name of pressure field, default is "p"
        word pName_;

        //- Density value
        scalar rho_;


    // Private Member Functions

        //- Return true if the pressure field corresponds to kinematic pressure
        bool isKinematicPressure();

        //- Disallow default bitwise copy construct
        staticPressure(const staticPressure&);

        //- Disallow default bitwise assignment
        void operator=(const staticPressure&);


public:

    //- Runtime type information
    TypeName("staticPressure");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        staticPressure
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    // Destructor

        virtual ~staticPressure();


    // Member Functions

        //- Return name of the set of staticPressure
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the staticPressure data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Calculate the staticPressure and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
