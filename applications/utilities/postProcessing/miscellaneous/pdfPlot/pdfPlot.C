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

Description
    Generates an .obj file to plot a probability distribution function

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pdf.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createFields.H"

    OFstream objFileNew("s.m");

    Random rndGen(label(0));

    autoPtr<pdf> p
    (
        pdf::New(pdfDictionary, rndGen)
    );

    scalar xMin = p->minValue();

    scalar xMax = p->maxValue();

    for(label i=0;i<nSamples;i++)
    {
        scalar ps = p->sample();
        label n = label((ps-xMin)*nIntervals/(xMax-xMin));
        //Info << "p[" << i << "] = " << ps << ", n = " << n << endl;
        samples[n]++;
    }

    for(label i=0;i<nIntervals;i++)
    {
        scalar x = xMin + i*(xMax-xMin)/(nIntervals-1);
        objFileNew << x << " \t" << samples[i] << endl;
    }
    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
