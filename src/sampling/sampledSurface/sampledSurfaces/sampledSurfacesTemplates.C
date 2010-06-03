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

#include "sampledSurfaces.H"
#include "volFields.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::sampledSurfaces::sampleAndWrite
(
    const GeometricField<Type, fvPatchField, volMesh>& vField,
    const surfaceWriter<Type>& formatter
)
{
    // interpolator for this field
    autoPtr< interpolation<Type> > interpolator;

    const word& fieldName   = vField.name();
    const fileName outputDir = outputPath_/vField.time().timeName();

    forAll(*this, surfI)
    {
        const sampledSurface& s = operator[](surfI);

        Field<Type> values;

        if (s.interpolate())
        {
            if (interpolator.empty())
            {
                interpolator = interpolation<Type>::New
                (
                    interpolationScheme_,
                    vField
                );
            }

            values = s.interpolate(interpolator());
        }
        else
        {
            values = s.sample(vField);
        }

        if (Pstream::parRun())
        {
            // Collect values from all processors
            List<Field<Type> > gatheredValues(Pstream::nProcs());
            gatheredValues[Pstream::myProcNo()] = values;
            Pstream::gatherList(gatheredValues);

            if (Pstream::master())
            {
                // Combine values into single field
                Field<Type> allValues
                (
                    ListListOps::combine<Field<Type> >
                    (
                        gatheredValues,
                        accessOp<Field<Type> >()
                    )
                );

                // Renumber (point data) to correspond to merged points
                if (mergeList_[surfI].pointsMap.size() == allValues.size())
                {
                    inplaceReorder(mergeList_[surfI].pointsMap, allValues);
                    allValues.setSize(mergeList_[surfI].points.size());
                }

                // Write to time directory under outputPath_
                // skip surface without faces (eg, a failed cut-plane)
                if (mergeList_[surfI].faces.size())
                {
                    formatter.write
                    (
                        outputDir,
                        s.name(),
                        mergeList_[surfI].points,
                        mergeList_[surfI].faces,
                        fieldName,
                        allValues
                    );
                }
            }
        }
        else
        {
            // Write to time directory under outputPath_
            // skip surface without faces (eg, a failed cut-plane)
            if (s.faces().size())
            {
                formatter.write
                (
                    outputDir,
                    s.name(),
                    s.points(),
                    s.faces(),
                    fieldName,
                    values
                );
            }
        }
    }
}


template<class Type>
void Foam::sampledSurfaces::sampleAndWrite
(
    fieldGroup<Type>& fields
)
{
    if (fields.size())
    {
        // create or use existing surfaceWriter
        if (fields.formatter.empty())
        {
            fields.formatter = surfaceWriter<Type>::New(writeFormat_);
        }

        forAll(fields, fieldI)
        {
            if (Pstream::master() && verbose_)
            {
                Pout<< "sampleAndWrite: " << fields[fieldI] << endl;
            }

            if (loadFromFiles_)
            {
                sampleAndWrite
                (
                    GeometricField<Type, fvPatchField, volMesh>
                    (
                        IOobject
                        (
                            fields[fieldI],
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE,
                            false
                        ),
                        mesh_
                    ),
                    fields.formatter()
                );
            }
            else
            {
                objectRegistry::const_iterator iter =
                    mesh_.find(fields[fieldI]);

                if
                (
                    iter != mesh_.objectRegistry::end()
                 && iter()->type()
                 == GeometricField<Type, fvPatchField, volMesh>::typeName
                )
                {
                   sampleAndWrite
                   (
                       mesh_.lookupObject
                       <GeometricField<Type, fvPatchField, volMesh> >
                       (
                           fields[fieldI]
                       ),
                       fields.formatter()
                   );
                }
            }
        }
    }
}


// ************************************************************************* //
