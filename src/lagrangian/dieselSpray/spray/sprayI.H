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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Time& spray::runTime() const
{
    return runTime_;
}


inline const fvMesh& spray::mesh() const
{
    return mesh_;
}


inline const volVectorField& spray::U() const
{
    return U_;
}


inline const volScalarField& spray::rho() const
{
    return rho_;
}


inline const volScalarField& spray::p() const
{
    return p_;
}


inline const volScalarField& spray::T() const
{
    return T_;
}


inline PtrList<injector>& spray::injectors()
{
    return injectors_;
}


inline const PtrList<injector>& spray::injectors() const
{
    return injectors_;
}


inline const atomizationModel& spray::atomization() const
{
    return atomization_;
}


inline const breakupModel& spray::breakup() const
{
    return breakupModel_;
}


inline const collisionModel& spray::collisions() const
{
    return collisionModel_;
}


inline const dispersionModel& spray::dispersion() const
{
    return dispersionModel_;
}


inline const dragModel& spray::drag() const
{
    return drag_;
}


inline const evaporationModel& spray::evaporation() const
{
    return evaporation_;
}


inline const heatTransferModel& spray::heatTransfer() const
{
    return heatTransfer_;
}


inline const injectorModel& spray::injection() const
{
    return injectorModel_;
}


inline const wallModel& spray::wall() const
{
    return wall_;
}


inline tmp<volVectorField> spray::momentumSource() const
{
    tmp<volVectorField> tsource
    (
        new volVectorField
        (
            IOobject
            (
                "sms",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "zero",
                dimensionSet(1, -2, -2, 0, 0),
                vector::zero
            )
        )
    );

    tsource().internalField() = sms_/runTime_.deltaT().value()/mesh_.V();

    return tsource;
}


inline tmp<volScalarField> spray::evaporationSource(const label si) const
{
    tmp<volScalarField> tsource
    (
        new volScalarField
        (
            IOobject
            (
                "srhos",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.0)
        )
    );

    if (isLiquidFuel_[si])
    {
        label fi = gasToLiquidIndex_[si];
        tsource().internalField() = srhos_[fi]/runTime_.deltaT().value()/mesh_.V();
    }
    else
    {
        scalarField s(mesh_.nCells(), 0.0);
        tsource().internalField() = s;
    }

    return tsource;
}


inline tmp<volScalarField> spray::heatTransferSource() const
{
    tmp<volScalarField> tsource
    (
        new volScalarField
        (
            IOobject
            (
                "shs",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );

    tsource().internalField() = shs_/runTime_.deltaT().value()/mesh_.V();

    return tsource;
}


inline Random& spray::rndGen()
{
    return rndGen_;
}


inline label spray::subCycles() const
{
    return subCycles_;
}


inline const vector& spray::g() const
{
    return g_;
}


inline const liquidMixture& spray::fuels() const
{
    return fuels_;
}


inline const PtrList<gasThermoPhysics>& spray::gasProperties() const
{
    return gasProperties_;
}


inline const basicMultiComponentMixture& spray::composition() const
{
    return composition_;
}


inline const List<label>& spray::liquidToGasIndex() const
{
    return liquidToGasIndex_;
}


inline const List<label>& spray::gasToLiquidIndex() const
{
    return gasToLiquidIndex_;
}


inline const List<bool>& spray::isLiquidFuel() const
{
    return isLiquidFuel_;
}


inline const bool& spray::twoD() const
{
    return twoD_;
}


inline const vector& spray::axisOfSymmetry() const
{
    return axisOfSymmetry_;
}


inline const vector& spray::axisOfWedge() const
{
    return axisOfWedge_;
}


inline const vector& spray::axisOfWedgeNormal() const
{
    return axisOfWedgeNormal_;
}


inline const scalar& spray::angleOfWedge() const
{
    return angleOfWedge_;
}


inline const interpolation<vector>& spray::UInterpolator() const
{
    return UInterpolator_;
}


inline const interpolation<scalar>& spray::rhoInterpolator() const
{
    return rhoInterpolator_;
}


inline const interpolation<scalar>& spray::pInterpolator() const
{
    return pInterpolator_;
}


inline const interpolation<scalar>& spray::TInterpolator() const
{
    return TInterpolator_;
}


inline vectorField& spray::sms()
{
    return sms_;
}


inline const vectorField& spray::sms() const
{
    return sms_;
}


inline scalarField& spray::shs()
{
    return shs_;
}


inline const scalarField& spray::shs() const
{
    return shs_;
}


inline PtrList<scalarField>& spray::srhos()
{
    return srhos_;
}


inline const PtrList<scalarField>& spray::srhos() const
{
    return srhos_;
}


inline const scalar& spray::ambientPressure() const
{
    return ambientPressure_;
}


inline const scalar& spray::ambientTemperature() const
{
    return ambientTemperature_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
