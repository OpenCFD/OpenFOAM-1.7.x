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

\*---------------------------------------------------------------------------*/

#include "fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeFvPatchField(fvPatchTypeField)                                    \
                                                                              \
defineNamedTemplateTypeNameAndDebug(fvPatchTypeField, 0);                     \
template<>                                                                    \
int fvPatchTypeField::disallowGenericFvPatchField                             \
(                                                                             \
    debug::debugSwitch("disallowGenericFvPatchField", 0)                      \
);                                                                            \
defineTemplateRunTimeSelectionTable(fvPatchTypeField, patch);                 \
defineTemplateRunTimeSelectionTable(fvPatchTypeField, patchMapper);           \
defineTemplateRunTimeSelectionTable(fvPatchTypeField, dictionary);

makeFvPatchField(fvPatchScalarField)
makeFvPatchField(fvPatchVectorField)
makeFvPatchField(fvPatchSphericalTensorField)
makeFvPatchField(fvPatchSymmTensorField)
makeFvPatchField(fvPatchTensorField)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
