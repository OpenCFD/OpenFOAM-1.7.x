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
    Foam::PhaseChangeModel

Description
    Templated phase change model class

SourceFiles
    PhaseChangeModel.C
    NewPhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef PhaseChangeModel_H
#define PhaseChangeModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class PhaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PhaseChangeModel
{
public:

    // Public enumerations

        //- Enthalpy transfer type
        enum enthalpyTransferType
        {
            etLatentHeat,
            etEnthalpyDifference
        };

        //- Name representations of enthalpy transfer types
        static const Foam::wordList enthalpyTransferTypeNames;


protected:

    // Protected data

        //- The cloud dictionary
        const dictionary& dict_;

        //- Reference to the owner cloud class
        CloudType& owner_;

        //- The coefficient dictionary
        const dictionary coeffDict_;

        //- Enthalpy transfer type enumeration
        enthalpyTransferType enthalpyTransfer_;


    // Protected member functions

        //- Convert word to enthalpy transfer type
        enthalpyTransferType wordToEnthalpyTransfer(const word& etName) const;

        //- Sherwood number
        scalar Sh() const;


public:

    //- Runtime type information
    TypeName("PhaseChangeModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PhaseChangeModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        PhaseChangeModel(CloudType& owner);

        //- Construct from dictionary
        PhaseChangeModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~PhaseChangeModel();


    //- Selector
    static autoPtr<PhaseChangeModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Access

        //- Return the owner cloud object
        const CloudType& owner() const;

        //- Return the cloud dictionary
        const dictionary& dict() const;

        //- Return the coefficient dictionary
        const dictionary& coeffDict() const;

        //- Return the enthalpy transfer type enumeration
        const enthalpyTransferType& enthalpyTransfer() const;


    // Member Functions

        //- Flag to indicate whether model activates phase change model
        virtual bool active() const = 0;

        //- Update model
        virtual void calculate
        (
            const scalar dt,
            const label cellI,
            const scalar Re,
            const scalar d,
            const scalar nu,
            const scalar T,
            const scalar Ts,
            const scalar pc,
            scalarField& dMassPC
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePhaseChangeModel(CloudType)                                       \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(PhaseChangeModel<CloudType>, 0);      \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        PhaseChangeModel<CloudType>,                                          \
        dictionary                                                            \
    );


#define makePhaseChangeModelThermoType(SS, CloudType, ParcelType, ThermoType) \
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        SS<CloudType<ParcelType<ThermoType> > >,                              \
        0                                                                     \
    );                                                                        \
                                                                              \
    PhaseChangeModel<CloudType<ParcelType<ThermoType> > >::                   \
        adddictionaryConstructorToTable                                       \
            <SS<CloudType<ParcelType<ThermoType> > > >                        \
            add##SS##CloudType##ParcelType##ThermoType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PhaseChangeModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
