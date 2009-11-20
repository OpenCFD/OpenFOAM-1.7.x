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

#include "OutputFilterFunctionObject.H"
#include "IOOutputFilter.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

template<class OutputFilter>
void Foam::OutputFilterFunctionObject<OutputFilter>::readDict()
{
    dict_.readIfPresent("region", regionName_);
    dict_.readIfPresent("dictionary", dictName_);
    dict_.readIfPresent("enabled", enabled_);
    dict_.readIfPresent("storeFilter", storeFilter_);
}

template<class OutputFilter>
void Foam::OutputFilterFunctionObject<OutputFilter>::allocateFilter()
{
    if (dictName_.size())
    {
        ptr_.reset
        (
            new IOOutputFilter<OutputFilter>
            (
                name(),
                time_.lookupObject<objectRegistry>(regionName_),
                dictName_
            )
        );
    }
    else
    {
        ptr_.reset
        (
            new OutputFilter
            (
                name(),
                time_.lookupObject<objectRegistry>(regionName_),
                dict_
            )
        );
    }
}

template<class OutputFilter>
void Foam::OutputFilterFunctionObject<OutputFilter>::destroyFilter()
{
    ptr_.reset();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class OutputFilter>
Foam::OutputFilterFunctionObject<OutputFilter>::OutputFilterFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    time_(t),
    dict_(dict),
    regionName_(polyMesh::defaultRegion),
    dictName_(),
    enabled_(true),
    storeFilter_(true),
    outputControl_(t, dict)
{
    readDict();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class OutputFilter>
void Foam::OutputFilterFunctionObject<OutputFilter>::on()
{
    enabled_ = true;
}


template<class OutputFilter>
void Foam::OutputFilterFunctionObject<OutputFilter>::off()
{
    enabled_ = false;
}


template<class OutputFilter>
bool Foam::OutputFilterFunctionObject<OutputFilter>::start()
{
    readDict();

    if (enabled_&&storeFilter_)
    {
        allocateFilter();
    }

    return true;
}


template<class OutputFilter>
bool Foam::OutputFilterFunctionObject<OutputFilter>::execute()
{
    if (enabled_)
    {
        if (!storeFilter_)
        {
            allocateFilter();
        }

        ptr_->execute();

        if (enabled_ && outputControl_.output())
        {
            ptr_->write();
        }

        if (!storeFilter_)
        {
            destroyFilter();
        }
    }

    return true;
}


template<class OutputFilter>
bool Foam::OutputFilterFunctionObject<OutputFilter>::end()
{
    if (enabled_)
    {
        if (!storeFilter_)
        {
            allocateFilter();
        }

        ptr_->end();

        if (enabled_ && outputControl_.output())
        {
            ptr_->write();
        }

        if (!storeFilter_)
        {
            destroyFilter();
        }
    }

    return true;
}


template<class OutputFilter>
bool Foam::OutputFilterFunctionObject<OutputFilter>::read
(
    const dictionary& dict
)
{
    if (dict != dict_)
    {
        dict_ = dict;
        outputControl_.read(dict);

        return start();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
