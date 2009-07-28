/*---------------------------------------------------------------------------* \
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
    Foam::IOOutputFilter

Description
    IOdictionary wrapper around OutputFilter to allow them to read from
    their associated dictionaries.

SourceFiles
    IOOutputFilter.C

\*---------------------------------------------------------------------------*/

#ifndef IOOutputFilter_H
#define IOOutputFilter_H

#include "IOdictionary.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class IOOutputFilter Declaration
\*---------------------------------------------------------------------------*/

template<class OutputFilter>
class IOOutputFilter
:
    public IOdictionary,
    public OutputFilter
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        IOOutputFilter(const IOOutputFilter&);
        void operator=(const IOOutputFilter&);


public:

    // Constructors

        //- Construct for given objectRegistry and dictionary
        //  Allow dictionary to be optional
        //  Allow the possibility to load fields from files
        IOOutputFilter
        (
            const word& outputFilterName,
            const objectRegistry&,
            const fileName& dictName = OutputFilter::typeName() + "Dict",
            const IOobject::readOption rOpt = IOobject::MUST_READ,
            const bool loadFromFile = false
        );


    //- Destructor
    virtual ~IOOutputFilter();


    // Member Functions

        //- Return name
        virtual const word& name() const
        {
            return IOdictionary::name();
        }

        //- Read the probes
        virtual bool read();

        //- Sample and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm)
        {
            read();
            OutputFilter::updateMesh(mpm);
        }

        //- Update for changes of mesh
        virtual void movePoints(const pointField& points)
        {
            read();
            OutputFilter::movePoints(points);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IOOutputFilter.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
