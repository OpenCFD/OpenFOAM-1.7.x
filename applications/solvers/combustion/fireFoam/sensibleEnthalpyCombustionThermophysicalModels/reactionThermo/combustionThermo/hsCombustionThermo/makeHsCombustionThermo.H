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

InClass
    Foam::hsCombustionThermo

Description

\*---------------------------------------------------------------------------*/

#ifndef makeHsCombustionThermo_H
#define makeHsCombustionThermo_H

#include "addToRunTimeSelectionTable.H"
#include "basicSensiblePsiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeHsCombustionThermo(CThermo,MixtureThermo,Mixture,Transport,Thermo,EqnOfState) \
                                                                              \
typedef MixtureThermo                                                         \
    <Mixture<Transport<specieThermo<Thermo<EqnOfState> > > > >                \
    MixtureThermo##Mixture##Transport##Thermo##EqnOfState;                    \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    MixtureThermo##Mixture##Transport##Thermo##EqnOfState,                    \
    #MixtureThermo                                                            \
        "<"#Mixture"<"#Transport"<specieThermo<"#Thermo"<"#EqnOfState">>>>>", \
    0                                                                         \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    basicSensiblePsiThermo,                                                   \
    MixtureThermo##Mixture##Transport##Thermo##EqnOfState,                    \
    fvMesh                                                                    \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    CThermo,                                                                  \
    MixtureThermo##Mixture##Transport##Thermo##EqnOfState,                    \
    fvMesh                                                                    \
)


#define makeHsCombustionMixtureThermo(CThermo,MixtureThermo,Mixture,ThermoPhys) \
                                                                              \
typedef MixtureThermo<Mixture<ThermoPhys> >                                   \
    MixtureThermo##Mixture##ThermoPhys;                                       \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    MixtureThermo##Mixture##ThermoPhys,                                       \
    #MixtureThermo"<"#Mixture"<"#ThermoPhys">>",                              \
    0                                                                         \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    basicSensiblePsiThermo,                                                   \
    MixtureThermo##Mixture##ThermoPhys,                                       \
    fvMesh                                                                    \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    CThermo,                                                                  \
    MixtureThermo##Mixture##ThermoPhys,                                       \
    fvMesh                                                                    \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
