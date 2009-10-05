/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "greyMeanAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(greyMeanAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            greyMeanAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyMeanAbsorptionEmission::greyMeanAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.subDict(typeName + "Coeffs"))),
    speciesNames_(0),
    specieIndex_(0),
    lookUpTable_
    (
        fileName(coeffsDict_.lookup("lookUpTableFileName")),
        mesh.time().constant(),
        mesh
    ),
    thermo_(mesh.lookupObject<basicThermo>("thermophysicalProperties")),
    EhrrCoeff_(readScalar(coeffsDict_.lookup("EhrrCoeff"))),
    Yj_(nSpecies_)
{
    label nFunc = 0;
    const dictionary& functionDicts = dict.subDict(typeName + "Coeffs");

    forAllConstIter(dictionary, functionDicts, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }
        const word& key = iter().keyword();
        speciesNames_.insert(key, nFunc);
        const dictionary& dict = iter().dict();
        coeffs_[nFunc].initialise(dict);
        nFunc++;
    }

    // Check that all the species on the dictionary are present in the
    // look-up table and save the corresponding indices of the look-up table

    label j = 0;
    forAllConstIter(HashTable<label>, speciesNames_, iter)
    {
        if (mesh.foundObject<volScalarField>("ft"))
        {
            if (lookUpTable_.found(iter.key()))
            {
                label index = lookUpTable_.findFieldIndex(iter.key());

                Info<< "specie: " << iter.key() << " found on look-up table "
                    << " with index: " << index << endl;

                specieIndex_[iter()] = index;
            }
            else if (mesh.foundObject<volScalarField>(iter.key()))
            {
                volScalarField& Y =
                    const_cast<volScalarField&>
                    (
                        mesh.lookupObject<volScalarField>(iter.key())
                    );
                Yj_.set(j, &Y);
                specieIndex_[iter()] = 0;
                j++;
                Info<< "specie: " << iter.key() << " is being solved" << endl;
            }
            else
            {
                FatalErrorIn
                (
                    "Foam::radiation::greyMeanAbsorptionEmission(const"
                    "dictionary& dict, const fvMesh& mesh)"
                )   << "specie: " << iter.key()
                    << " is neither in look-up table: "
                    << lookUpTable_.tableName()
                    << " nor is being solved" << nl
                    << exit(FatalError);
            }
        }
        else
        {
            FatalErrorIn
            (
                "Foam::radiation::greyMeanAbsorptionEmission(const"
                "dictionary& dict, const fvMesh& mesh)"
            )   << "specie ft is not present " << nl
                << exit(FatalError);

        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::greyMeanAbsorptionEmission::~greyMeanAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanAbsorptionEmission::aCont(const label bandI) const
{
    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();
    const volScalarField& ft = mesh_.lookupObject<volScalarField>("ft");

    label nSpecies = speciesNames_.size();

    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0.0)
        )
    );

    scalarField& a = ta().internalField();

    forAll(a, i)
    {
        const List<scalar>& species = lookUpTable_.lookUp(ft[i]);

        for (label n=0; n<nSpecies; n++)
        {
            label l = 0;
            scalar Yipi = 0;
            if (specieIndex_[n] != 0)
            {
                //moles x pressure [atm]
                Yipi = species[specieIndex_[n]]*p[i]*9.869231e-6;
            }
            else
            {
                // mass fraction
                Yipi = Yj_[l][i];
                l++;
            }

            const absorptionCoeffs::coeffArray& b = coeffs_[n].coeffs(T[i]);

            scalar Ti = T[i];
            // negative temperature exponents
            if (coeffs_[n].invTemp())
            {
                Ti = 1./T[i];
            }
            a[i] +=
                Yipi
               *(
                    ((((b[5]*Ti + b[4])*Ti + b[3])*Ti + b[2])*Ti + b[1])*Ti
                  + b[0]
                );
        }
    }
    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanAbsorptionEmission::eCont(const label bandI) const
{
    tmp<volScalarField> e
    (
        new volScalarField
        (
            IOobject
            (
                "e",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("e", dimless/dimLength, 0.0)
        )
    );

    return e;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanAbsorptionEmission::ECont(const label bandI) const
{
    tmp<volScalarField> E
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    if (mesh_.foundObject<volScalarField>("dQ"))
    {
        const volScalarField& dQ =
            mesh_.lookupObject<volScalarField>("dQ");
        E().internalField() = EhrrCoeff_*dQ;
    }

    return E;
}


// ************************************************************************* //
