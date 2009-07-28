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

Application
    boundaryFoam

Description
    Steady-state solver for 1D turbulent flow, typically to generate boundary
    layer conditions at an inlet, for use in a simulation.

    Boundary layer code to calculate the U, k and epsilon distributions.
    Used to create inlet boundary conditions for experimental comparisons
    for which U and k have not been measured.
    Turbulence model is runtime selectable.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "wallFvPatch.H"
#include "makeGraph.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        fvVectorMatrix divR = turbulence->divDevReff(U);
        divR.source() = flowMask & divR.source();

        fvVectorMatrix UEqn
        (
            divR == gradP
        );

        UEqn.relax();

        UEqn.solve();


        // Correct driving force for a constant mass flow rate

        dimensionedVector UbarStar = flowMask & U.weightedAverage(mesh.V());

        U += (Ubar - UbarStar);
        gradP += (Ubar - UbarStar)/(1.0/UEqn.A())().weightedAverage(mesh.V());

        label id = y.size() - 1;

        scalar wallShearStress =
            flowDirection & turbulence->R()()[id] & wallNormal;

        scalar yplusWall
//            = ::sqrt(mag(wallShearStress))*y[id]/laminarTransport.nu()()[id];
            = ::sqrt(mag(wallShearStress))*y[id]/turbulence->nuEff()()[id];

        Info<< "Uncorrected Ubar = " << (flowDirection & UbarStar.value())<< tab
            << "pressure gradient = " << (flowDirection & gradP.value()) << tab
            << "min y+ = " << yplusWall << endl;


        turbulence->correct();


        if (runTime.outputTime())
        {
            volSymmTensorField R
            (
                IOobject
                (
                    "R",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                turbulence->R()
            );

            runTime.write();

            const word& gFormat = runTime.graphFormat();

            makeGraph(y, flowDirection & U, "Uf", gFormat);

            makeGraph(y, laminarTransport.nu(), gFormat);

            makeGraph(y, turbulence->k(), gFormat);
            makeGraph(y, turbulence->epsilon(), gFormat);

            //makeGraph(y, flowDirection & R & flowDirection, "Rff", gFormat);
            //makeGraph(y, wallNormal & R & wallNormal, "Rww", gFormat);
            //makeGraph(y, flowDirection & R & wallNormal, "Rfw", gFormat);

            //makeGraph(y, sqrt(R.component(tensor::XX)), "u", gFormat);
            //makeGraph(y, sqrt(R.component(tensor::YY)), "v", gFormat);
            //makeGraph(y, sqrt(R.component(tensor::ZZ)), "w", gFormat);
            makeGraph(y, R.component(tensor::XY), "uv", gFormat);

            makeGraph(y, mag(fvc::grad(U)), "gammaDot", gFormat);
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
