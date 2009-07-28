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

\*---------------------------------------------------------------------------*/

#include "sampledSets.H"
#include "volFields.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class Type>
Foam::sampledSets::volFieldSampler<Type>::volFieldSampler
(
    const word& interpolationScheme,
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const PtrList<sampledSet>& samplers
)
:
    List<Field<Type> >(samplers.size()),
    name_(field.name())
{
    autoPtr<interpolation<Type> > interpolator
    (
        interpolation<Type>::New(interpolationScheme, field)
    );

    forAll(samplers, seti)
    {
        Field<Type>& values = this->operator[](seti);
        const sampledSet& samples = samplers[seti];

        values.setSize(samples.size());
        forAll(samples, samplei)
        {
            const point& samplePt = samples[samplei];
            label celli = samples.cells()[samplei];
            label facei = samples.faces()[samplei];

            values[samplei] = interpolator().interpolate
            (
                samplePt,
                celli,
                facei
            );
        }
    }
}


template <class Type>
Foam::sampledSets::volFieldSampler<Type>::volFieldSampler
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const PtrList<sampledSet>& samplers
)
:
    List<Field<Type> >(samplers.size()),
    name_(field.name())
{
    forAll(samplers, seti)
    {
        Field<Type>& values = this->operator[](seti);
        const sampledSet& samples = samplers[seti];

        values.setSize(samples.size());
        forAll(samples, samplei)
        {
            values[samplei] = field[samples.cells()[samplei]];
        }
    }
}


template <class Type>
Foam::sampledSets::volFieldSampler<Type>::volFieldSampler
(
    const List<Field<Type> >& values,
    const word& name
)
:
    List<Field<Type> >(values),
    name_(name)
{}


template<class Type>
Foam::label Foam::sampledSets::grep
(
    fieldGroup<Type>& fieldList,
    const wordList& fieldTypes
) const
{
    fieldList.setSize(fieldNames_.size());
    label nFields = 0;

    forAll(fieldNames_, fieldi)
    {
        if
        (
            fieldTypes[fieldi]
         == GeometricField<Type, fvPatchField, volMesh>::typeName
        )
        {
            fieldList[nFields] = fieldNames_[fieldi];
            nFields++;
        }
    }

    fieldList.setSize(nFields);

    return nFields;
}


template<class Type>
void Foam::sampledSets::writeSampleFile
(
    const coordSet& masterSampleSet,
    const PtrList<volFieldSampler<Type> >& masterFields,
    const label seti,
    const fileName& timeDir,
    const writer<Type>& formatter
)
{
    wordList valueSetNames(masterFields.size());
    List<const Field<Type>*> valueSets(masterFields.size());

    forAll(masterFields, fieldi)
    {
        valueSetNames[fieldi] = masterFields[fieldi].name();
        valueSets[fieldi] = &masterFields[fieldi][seti];
    }

    fileName fName
    (
        timeDir/formatter.getFileName(masterSampleSet, valueSetNames)
    );

    formatter.write
    (
        masterSampleSet,
        valueSetNames,
        valueSets,
        OFstream(fName)()
    );
}


template<class T>
void Foam::sampledSets::combineSampledValues
(
    const PtrList<volFieldSampler<T> >& sampledFields,
    const labelListList& indexSets,
    PtrList<volFieldSampler<T> >& masterFields
)
{
    forAll(sampledFields, fieldi)
    {
        List<Field<T> > masterValues(indexSets.size());

        forAll(indexSets, seti)
        {
            // Collect data from all processors
            List<Field<T> > gatheredData(Pstream::nProcs());
            gatheredData[Pstream::myProcNo()] = sampledFields[fieldi][seti];
            Pstream::gatherList(gatheredData);

            if (Pstream::master())
            {
                Field<T> allData
                (
                    ListListOps::combine<Field<T> >
                    (
                        gatheredData,
                        Foam::accessOp<Field<T> >()
                    )
                );

                masterValues[seti] = UIndirectList<T>
                (
                    allData,
                    indexSets[seti]
                )();
            }
        }

        masterFields.set
        (
            fieldi,
            new volFieldSampler<T>
            (
                masterValues,
                sampledFields[fieldi].name()
            )
        );
    }
}


template<class Type>
void Foam::sampledSets::sampleAndWrite
(
    fieldGroup<Type>& fields
)
{
    if (fields.size())
    {
        bool interpolate = interpolationScheme_ != "cell";

        // Create or use existing writer
        if (fields.formatter.empty())
        {
            fields.formatter = writer<Type>::New(writeFormat_);
        }

        // Storage for interpolated values
        PtrList<volFieldSampler<Type> > sampledFields(fields.size());

        forAll(fields, fieldi)
        {
            if (Pstream::master() && verbose_)
            {
                Pout<< "sampledSets::sampleAndWrite: "
                    << fields[fieldi] << endl;
            }

            if (loadFromFiles_)
            {
                GeometricField<Type, fvPatchField, volMesh> vf
                (
                    IOobject
                    (
                        fields[fieldi],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                );

                if (interpolate)
                {
                    sampledFields.set
                    (
                        fieldi,
                        new volFieldSampler<Type>
                        (
                            interpolationScheme_,
                            vf,
                            *this
                        )
                    );
                }
                else
                {
                    sampledFields.set
                    (
                        fieldi,
                        new volFieldSampler<Type>(vf, *this)
                    );
                }
            }
            else
            {
                if (interpolate)
                {
                    sampledFields.set
                    (
                        fieldi,
                        new volFieldSampler<Type>
                        (
                            interpolationScheme_,
                            mesh_.lookupObject
                            <GeometricField<Type, fvPatchField, volMesh> >
                            (fields[fieldi]),
                            *this
                        )
                    );
                }
                else
                {
                    sampledFields.set
                    (
                        fieldi,
                        new volFieldSampler<Type>
                        (
                            mesh_.lookupObject
                            <GeometricField<Type, fvPatchField, volMesh> >
                            (fields[fieldi]),
                            *this
                        )
                    );
                }
            }
        }

        // Combine sampled fields from processors.
        // Note: only master results are valid

        PtrList<volFieldSampler<Type> > masterFields(sampledFields.size());
        combineSampledValues(sampledFields, indexSets_, masterFields);

        if (Pstream::master())
        {
            forAll(masterSampledSets_, seti)
            {
                writeSampleFile
                (
                    masterSampledSets_[seti],
                    masterFields,
                    seti,
                    outputPath_/mesh_.time().timeName(),
                    fields.formatter()
                );
            }
        }
    }
}


// ************************************************************************* //
