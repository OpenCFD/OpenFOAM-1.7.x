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
    Foam::forces

Description
    Calculates the forces and moments by integrating the pressure and
    skin-friction forces over a given list of patches.

    Member function calcForcesMoment()calculates and returns the forces and
    moments.

    Member function forces::write() calls calcForcesMoment() and writes the
    forces and moments into the file \<timeDir\>/forces.dat

SourceFiles
    forces.C
    IOforces.H

\*---------------------------------------------------------------------------*/

#ifndef forces_H
#define forces_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class forces Declaration
\*---------------------------------------------------------------------------*/

class forces
{
public:

    // Tuple which holds the pressure (.first()) and viscous (.second) forces
    typedef Tuple2<vector, vector> pressureViscous;

    // Tuple which holds the forces (.first()) and moment (.second)
    // pressure/viscous forces Tuples.
    typedef Tuple2<pressureViscous, pressureViscous> forcesMoments;

    //- Sum operation class to accumulate the pressure, viscous forces and moments
    class sumOp
    {
    public:

        forcesMoments operator()
        (
            const forcesMoments& fm1,
            const forcesMoments& fm2
        ) const
        {
            return forcesMoments
            (
                pressureViscous
                (
                    fm1.first().first() + fm2.first().first(),
                    fm1.first().second() + fm2.first().second()
                ),
                pressureViscous
                (
                    fm1.second().first() + fm2.second().first(),
                    fm1.second().second() + fm2.second().second()
                )
            );
        }
    };


protected:

    // Private data

        //- Name of this set of forces,
        //  Also used as the name of the probes directory.
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        // Read from dictionary

            //- Patches to integrate forces over
            labelHashSet patchSet_;

            //- Name of pressure field
            word pName_;

            //- Name of velocity field
            word UName_;

            //- Name of density field (optional)
            word rhoName_;

            //- Is the force density being supplied directly?
            Switch directForceDensity_;

            //- The name of the force density (fD) field
            word fDName_;

            //- Reference density needed for incompressible calculations
            scalar rhoRef_;

            //- Centre of rotation
            vector CofR_;


        //- Forces/moment file ptr
        autoPtr<OFstream> forcesFilePtr_;


    // Private Member Functions

        //- If the forces file has not been created create it
        void makeFile();

        //- Return the effective viscous stress (laminar + turbulent).
        tmp<volSymmTensorField> devRhoReff() const;

        //- Return rho if rhoName is specified otherwise rhoRef
        tmp<volScalarField> rho() const;

        //- Return rhoRef if the pressure field is dynamic, i.e. p/rho
        //  otherwise return 1
        scalar rho(const volScalarField& p) const;

        //- Disallow default bitwise copy construct
        forces(const forces&);

        //- Disallow default bitwise assignment
        void operator=(const forces&);

        //- Output file header information
        virtual void writeFileHeader();


public:

    //- Runtime type information
    TypeName("forces");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        forces
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~forces();


    // Member Functions

        //- Return name of the set of forces
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the forces data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Write the forces
        virtual void write();

        //- Calculate and return forces and moment
        virtual forcesMoments calcForcesMoment() const;

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
