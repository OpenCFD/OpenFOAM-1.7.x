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

\*---------------------------------------------------------------------------*/

if (writeVacf)
{
    OFstream vacfFile(runTime.path()/"vacf");

    if (!vacf.writeAveraged(vacfFile))
    {
        FatalErrorIn(args.executable())
            << "Failed writing to "
            << vacfFile.name()
            << abort(FatalError);
    }
}

Info<< "Diffusion coefficient = "
    << vacf.integral() << endl;

if (writePacf)
{
    OFstream pacfFile(runTime.path()/"pacf");

    if (!pacf.writeAveraged(pacfFile))
    {
        FatalErrorIn(args.executable())
            << "Failed writing to "
            << pacfFile.name()
            << abort(FatalError);
    }
}

Info<< "Viscosity = "
    << pacf.integral()/averageTemperature/moleculeCloud::kb/meshVolume
    << endl;

if (writeHFacf)
{
    OFstream hfacfFile
    (
        runTime.path()/ + "hfacf"
    );

    if (!hfacf.writeAveraged(hfacfFile))
    {
        FatalErrorIn(args.executable())
            << "Failed writing to "
            << hfacfFile.name()
            << abort(FatalError);
    }
}

Info<< "Thermal conductivity = "
    << hfacf.integral()
        /averageTemperature
        /averageTemperature
        /moleculeCloud::kb
        / meshVolume
    << endl;
