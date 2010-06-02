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
    Foam::fvSchemes

Description
    Selector class for finite volume differencing schemes.
    fvMesh is derived from fvShemes so that all fields have access to the
    fvSchemes from the mesh reference they hold.

SourceFiles
    fvSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef fvSchemes_H
#define fvSchemes_H

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fvSchemes Declaration
\*---------------------------------------------------------------------------*/

class fvSchemes
:
    public IOdictionary
{
    // Private data

        dictionary ddtSchemes_;
        ITstream defaultDdtScheme_;

        dictionary d2dt2Schemes_;
        ITstream defaultD2dt2Scheme_;

        dictionary interpolationSchemes_;
        ITstream defaultInterpolationScheme_;

        dictionary divSchemes_;
        ITstream defaultDivScheme_;

        dictionary gradSchemes_;
        ITstream defaultGradScheme_;

        dictionary snGradSchemes_;
        ITstream defaultSnGradScheme_;

        dictionary laplacianSchemes_;
        ITstream defaultLaplacianScheme_;

        dictionary fluxRequired_;
        bool defaultFluxRequired_;

        dictionary cacheFields_;


    // Private Member Functions

        //- Clear the dictionaries and streams before reading
        void clear();

        //- Disallow default bitwise copy construct
        fvSchemes(const fvSchemes&);

        //- Disallow default bitwise assignment
        void operator=(const fvSchemes&);


public:

    //- Debug switch
    static int debug;


    // Constructors

        //- Construct for objectRegistry
        fvSchemes(const objectRegistry& obr);


    // Member Functions

        // Access

            const dictionary& schemesDict() const;

            ITstream& ddtScheme(const word& name) const;

            ITstream& d2dt2Scheme(const word& name) const;

            ITstream& interpolationScheme(const word& name) const;

            ITstream& divScheme(const word& name) const;

            ITstream& gradScheme(const word& name) const;

            ITstream& snGradScheme(const word& name) const;

            ITstream& laplacianScheme(const word& name) const;

            bool fluxRequired(const word& name) const;

            bool cache(const word& name) const;


        // Read

            //- Read the fvSchemes
            bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
