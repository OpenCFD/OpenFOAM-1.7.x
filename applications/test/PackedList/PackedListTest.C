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

Application

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "uLabel.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "IFstream.H"
#include "PackedBoolList.H"
#include <climits>


using namespace Foam;

template<unsigned nBits>
inline void reportInfo()
{
    unsigned offset = PackedList<nBits>::packing();

    unsigned useSHL = ((1u << (nBits * offset)) - 1);
    unsigned useSHR = (~0u >> (sizeof(unsigned)*CHAR_BIT - nBits * offset));

    Info<< nl
        << "PackedList<" << nBits << ">" << nl
        << " max_value: " << PackedList<nBits>::max_value() << nl
        << " packing: " << PackedList<nBits>::packing() << nl
        << " utilization: " << (nBits * offset) << nl;

    Info<< " Masking:" << nl
        << "  shift << " << unsigned(nBits * offset) << nl
        << "  shift >> " << unsigned((sizeof(unsigned)*CHAR_BIT) - nBits * offset)
        << nl;

    hex(Info);
    Info<< "   maskLower: " << PackedList<nBits>::maskLower(PackedList<nBits>::packing())
        << nl
        << "      useSHL: " << useSHL << nl
        << "      useSHR: " << useSHR << nl;

    if (useSHL != useSHR)
    {
        Info<< "WARNING:  different results for SHL and SHR" << nl;
    }

    Info<< nl;
    dec(Info);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.insert("file .. fileN");

    argList::validOptions.insert("mask", "");
    argList::validOptions.insert("count", "");
    argList::validOptions.insert("info", "");

    argList args(argc, argv, false, true);


    if (args.optionFound("mask"))
    {
        Info<< "bit width: " << unsigned(sizeof(unsigned)*CHAR_BIT) << endl;
        reportInfo<1>();
        reportInfo<2>();
        reportInfo<3>();
        reportInfo<4>();
        reportInfo<5>();
        reportInfo<6>();
        reportInfo<7>();
        reportInfo<8>();
        reportInfo<9>();
        reportInfo<10>();
        reportInfo<11>();
        reportInfo<12>();
        reportInfo<13>();
        reportInfo<14>();
        reportInfo<15>();
        reportInfo<16>();
        reportInfo<17>();
        reportInfo<18>();
        reportInfo<19>();
        reportInfo<20>();
        reportInfo<21>();
        reportInfo<22>();
        reportInfo<23>();
        reportInfo<24>();
        reportInfo<25>();
        reportInfo<26>();
        reportInfo<27>();
        reportInfo<28>();
        reportInfo<29>();
        reportInfo<30>();
        reportInfo<31>();

        return 0;
    }
    else if (args.additionalArgs().empty())
    {
        args.printUsage();
    }


    forAll(args.additionalArgs(), argI)
    {
        const string& srcFile = args.additionalArgs()[argI];
        Info<< nl << "reading " << srcFile << nl;

        IFstream ifs(srcFile);
        List<label> rawLst(ifs);

        PackedBoolList packLst(rawLst);

        Info<< "size: " << packLst.size() << nl;

        if (args.optionFound("count"))
        {
            unsigned int rawCount = 0;
            forAll(rawLst, elemI)
            {
                if (rawLst[elemI])
                {
                    rawCount++;
                }
            }
            Info<< "raw count: " << rawCount << nl
                << "packed count: " << packLst.count() << nl;
        }

        if (args.optionFound("info"))
        {
            packLst.print(Info);
        }

        Info<< nl;
        IOobject::writeDivider(Info);
    }

    return 0;
}

// ************************************************************************* //
