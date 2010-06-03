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

Class
    Foam::parcel

Description
    A parcel of particles

\*---------------------------------------------------------------------------*/

#ifndef parcel_H
#define parcel_H

#include "Particle.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class spray;

/*---------------------------------------------------------------------------*\
                           Class parcel Declaration
\*---------------------------------------------------------------------------*/

class parcel
:
    public Particle<parcel>
{
    // Private member data

        // Reference to the names of the liquid components
            List<word> liquidComponents_;

        // Defining data (read and written to field files)

            //- Diameter of droplets in parcel
            scalar d_;

            //- Temperature of droplets in parcel
            scalar T_;

            //- Total parcel mass
            scalar m_;

            //- Spherical deviation of droplets in parcel
            scalar y_;

            //- Rate of change of spherical deviation of droplets in parcel
            scalar yDot_;

            //- Characteristic time of droplets in parcel
            scalar ct_;

            //- Stripped mass of droplets in parcel
            scalar ms_;

            //- Time spent in turbulent eddy
            scalar tTurb_;

            //- Part of liquid core (1-fully liquid, 0-droplet)
            scalar liquidCore_;

            //- injected from injector
            //  Should really be a label, but is scalar due to
            //  post-processing reasons
            scalar injector_;

            //- Velocity of parcel
            vector U_;

            //- Turbulent velocity fluctuation
            vector Uturb_;

            //- Normal n_ and axis describe the 2D plane
            //  in which the particle moves
            vector n_;

            //- Liquid components molar fractions
            scalarField X_;

        // Derived state information (not read or written)

            //- Momentum relaxation time of droplets in parcel
            scalar tMom_;


    // Private member functions

        //- Set the relaxation times
        void setRelaxationTimes
        (
            label celli,
            scalar& tauMomentum,
            scalarField& tauEvaporation,
            scalar& tauHeatTransfer,
            scalarField& tauBoiling,
            const spray& sprayDatabase,
            const scalar rho,
            const vector& Up,
            const scalar temperature,
            const scalar pressure,
            const scalarField& Yf,
            const scalarField& m0,
            const scalar dt
        );


        void updateParcelProperties
        (
            const scalar dt,
            spray& sprayData,
            const label celli,
            const label facei
        );


public:

    friend class Cloud<parcel>;


    // Constructors

        //- Construct from components
        parcel
        (
            const Cloud<parcel>& cloud,
            const vector& position,
            const label celli,
            const vector& n,
            const scalar d,
            const scalar T,
            const scalar m,
            const scalar y,
            const scalar yDot,
            const scalar ct,
            const scalar ms,
            const scalar tTurb,
            const scalar liquidCore,
            const scalar injector,
            const vector& U,
            const vector& Uturb,
            const scalarField& X,
            const List<word>& liquidNames
        );

        //- Construct from Istream reading field values if required
        parcel
        (
            const Cloud<parcel>& cloud,
            Istream& is,
            bool readFields = true
        );


    // Member Functions

        // Access

            //- Return the names of the liquid components
            inline const List<word>& liquidNames() const;

            //- Return the names of the liquid fuel components - identical with liquidNames
            inline const List<word>& fuelNames() const;

            //- Return diameter of droplets in parcel
            inline scalar d() const;

            //- Return diameter of droplets in parcel
            inline scalar& d();

            //- Return temperature of droplets in parcel
            inline scalar T() const;

            //- Return temperature of droplets in parcel
            inline scalar& T();

            //- Return total parcel mass
            inline scalar m() const;

            //- Return total parcel mass
            inline scalar& m();

            //- Return spherical deviation of droplets in parcel
            inline scalar dev() const;

            //- Return spherical deviation of droplets in parcel
            inline scalar& dev();

            //- Return rate of change of spherical deviation of
            //  droplets in parcel
            inline scalar ddev() const;

            //- Return rate of change of spherical deviation of
            //  droplets in parcel
            inline scalar& ddev();

            //- Return characteristic time of droplets in parcel
            inline scalar ct() const;

            //- Return characteristic time of droplets in parcel
            inline scalar& ct();

            //- Return stripped mass of droplets in parcel
            inline scalar& ms();

            //- Return stripped mass of droplets in parcel
            inline scalar ms() const;

            //- Return time spent in turbulent eddy
            inline scalar& tTurb();

            //- Return time spent in turbulent eddy
            inline scalar tTurb() const;

            //- Return part of liquid liquidCore
            inline scalar& liquidCore();

            //- Return part of liquid liquidCore
            inline scalar liquidCore() const;

            //- Return the injector from which is injected
            inline scalar& injector();

            //- Return the injector from which is injected
            inline scalar injector() const;

            //- Return velocity of parcel
            inline const vector& U() const;

            //- Return velocity of parcel
            inline vector& U();

            //- Return turbulent velocity fluctuation
            inline const vector& Uturb() const;

            //- Return turbulent velocity fluctuation
            inline vector& Uturb();

            //- Return the normal used for 2D purposes
            inline const vector& n() const;

            //- Return the normal used for 2D purposes
            inline vector& n();

            //- Return the liquid components molar fractions
            inline const scalarField& X() const;

            //- Return the liquid components molar fractions
            inline scalarField& X();

            //- Return the momentum relaxation time of droplets in parcel
            inline scalar& tMom();

            //- Return the momentum relaxation time of droplets in parcel
            inline scalar tMom() const;


        // Derived information

            //- Return statistical number of drops in parcel
            scalar N(const scalar rho) const;

            //- Return relative velocity between given vector and parcel
            inline vector Urel(const vector&) const;


        // Dimensionless Numbers

            //- Reynolds number based on rho an dynamic viscosity
            scalar Re
            (
                const scalar rho,
                const vector& U,
                const scalar mu
            ) const;

            //- Reynolds number based on kinematic viscosity
            scalar Re
            (
                const vector& U,
                const scalar nu
            ) const;

            //- Weber number
            scalar We
            (
                const vector& U,
                const scalar rho,
                const scalar sigma
            ) const;

            //- Schmidt number based on dynamic viscosity and rho
            scalar Sc
            (
                const scalar mu,
                const scalar rho,
                const scalar massDiffusion
            ) const;

            //- Schmidt number based on kinematic viscosity
            scalar Sc
            (
                const scalar nu,
                const scalar massDiffusion
            ) const;

            //- Prandtl number
            scalar Pr
            (
                const scalar cp,
                const scalar mu,
                const scalar kappa
            ) const;

            //- Volume of one droplet in the parcel
            scalar Vd() const;

            //- Volume of all droplets in parcel
            scalar V(const scalar rho) const;


        // Parcel operations

            bool move(spray& sprayData);

            //- Transform the position and physical properties of the particle
            //  according to the given transformation tensor
            void transformProperties(const tensor& T);

            //- Transform the position and physical properties of the particle
            //  according to the given separation vector
            void transformProperties(const vector& separation);

            //- fix the 2D plane normal,
            // when particle hits a face it is slightly perturbed
            // towards the face centre and n_ will no longer be valid
            inline void correctNormal(const vector& sym);


    // I/O

            static void readFields(Cloud<parcel>& c);

            static void writeFields(const Cloud<parcel>& c);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const parcel&);
};


template<>
inline bool contiguous<parcel>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "parcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
