/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::PostProcessingModel

Description
    Templated post-processing model class

SourceFiles
    PostProcessingModel.C
    NewPostProcessingModel.C

\*---------------------------------------------------------------------------*/

#ifndef PostProcessingModel_H
#define PostProcessingModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PostProcessingModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PostProcessingModel
{
    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        //- Reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;


protected:

    // Protected member functions

        //- Write post-processing info
        virtual void write() = 0;


public:

    //- Runtime type information
    TypeName("PostProcessingModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PostProcessingModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        PostProcessingModel(CloudType& owner);

        //- Construct from dictionary
        PostProcessingModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~PostProcessingModel();


    //- Selector
    static autoPtr<PostProcessingModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        // Access

            //- Return the owner cloud dictionary
            inline const dictionary& dict() const;

            //- Return const access the owner cloud object
            inline const CloudType& owner() const;

            //- Return non-const access the owner cloud object for manipulation
            inline CloudType& owner();

            //- Return the coefficients dictionary
            inline const dictionary& coeffDict() const;


        // Evaluation

            //- Main post-processing function
            virtual void post();

            //- Flag to indicate whether model activates post-processing model
            virtual bool active() const = 0;

            //- Gather post-processing info on patch
            virtual void postPatch
            (
                const typename CloudType::parcelType& p,
                const label patchI
            ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePostProcessingModel(CloudType)                                    \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(PostProcessingModel<CloudType>, 0);   \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        PostProcessingModel<CloudType>,                                       \
        dictionary                                                            \
    );


#define makePostProcessingModelType(SS, CloudType, ParcelType)                \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType<ParcelType> >, 0);       \
                                                                              \
    PostProcessingModel<CloudType<ParcelType> >::                             \
        adddictionaryConstructorToTable<SS<CloudType<ParcelType> > >          \
            add##SS##CloudType##ParcelType##ConstructorToTable_;


#define makePostProcessingModelThermoType(SS, CloudType, ParcelType, ThermoType)\
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        SS<CloudType<ParcelType<ThermoType> > >,                              \
        0                                                                     \
    );                                                                        \
                                                                              \
    PostProcessingModel<CloudType<ParcelType<ThermoType> > >::                \
        adddictionaryConstructorToTable                                       \
            <SS<CloudType<ParcelType<ThermoType> > > >                        \
            add##SS##CloudType##ParcelType##ThermoType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PostProcessingModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PostProcessingModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
