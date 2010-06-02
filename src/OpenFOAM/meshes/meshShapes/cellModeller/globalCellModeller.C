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

Description
    cellModeller global initializations

\*---------------------------------------------------------------------------*/

#include "cellModeller.H"
#include "OSspecific.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * Static data * * * * * * * * * * * * * * * * //


// PtrList of models
Foam::PtrList<Foam::cellModel> Foam::cellModeller::models_
(
    IFstream(findEtcFile("cellModels", true))()
);

// List of model pointers
Foam::List<Foam::cellModel*> Foam::cellModeller::modelPtrs_;

// HashTable of model pointers
Foam::HashTable<const Foam::cellModel*> Foam::cellModeller::modelDictionary_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Construct a dummy cellModeller which reads the models and fills
// the above tables
cellModeller globalCellModeller_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
