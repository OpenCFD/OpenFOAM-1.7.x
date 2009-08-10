/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
    Foam::faceZonesIntegration

Description
    Integrates surfaceScalarFields on faceZones

SourceFiles
    faceZonesIntegration.C
    IOfaceZonesIntegration.H

\*---------------------------------------------------------------------------*/

#ifndef faceZonesIntegration_H
#define faceZonesIntegration_H

#include "fvCFD.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "HashPtrTable.H"
#include "OFstream.H"
#include "Switch.H"
#include "pointFieldFwd.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class faceZonesIntegration Declaration
\*---------------------------------------------------------------------------*/

class faceZonesIntegration
{

protected:

    // Private data

        //- Name of this set of face zone integration,
        //  Also used as the name of the probes directory.
        word name_;

        const objectRegistry& obr_;

        // Read from dictionary

            //- on/off switch
            bool active_;

            //- Switch to send output to Info as well as to file
            Switch log_;

            //- faceZones to integrate over
            wordList faceZonesSet_;

            //- Names of the surface fields
            wordList fItems_;


        //- Current open files
        HashPtrTable<OFstream> faceZonesIntegrationFilePtr_;



    // Private Member Functions

        //- If the integration file has not been created create it
        void makeFile();

        scalar calcFaceZonesIntegral
        (
            const surfaceScalarField& fD,
            const faceZone& fz
        ) const;


        //- Disallow default bitwise copy construct
        faceZonesIntegration(const faceZonesIntegration&);

        //- Disallow default bitwise assignment
        void operator=(const faceZonesIntegration&);

        //- Output file header information
        virtual void writeFileHeader();


public:

    //- Runtime type information
    TypeName("faceZonesIntegration");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        faceZonesIntegration
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    // Destructor

        virtual ~faceZonesIntegration();


    // Member Functions

        //- Return name of the set of zones
        virtual const word& name() const
        {
            return name_;
        };

        //- Read the zone integration data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Write the integration
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
