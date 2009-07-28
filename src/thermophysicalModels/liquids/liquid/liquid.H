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

Class
    Foam::liquid

Description
    The thermophysical properties of a liquid

SourceFiles
    liquid.C

\*---------------------------------------------------------------------------*/

#ifndef liquid_H
#define liquid_H

#include "scalar.H"
#include "IOstreams.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class liquid Declaration
\*---------------------------------------------------------------------------*/

class liquid
{
    // Private data

        //- Molecular weight [kg/kmol]
        scalar W_;

        //- Critical temperature [K]
        scalar Tc_;

        //- Critical pressure [Pa]
        scalar Pc_;

        //- Critical volume [m^3/mol]
        scalar Vc_;

        //- Critical compressibility factor []
        scalar Zc_;

        //- Triple point temperature [K]
        scalar Tt_;

        //- Triple point pressure [Pa]
        scalar Pt_;

        //- Normal boiling temperature [K]
        scalar Tb_;

        //- Dipole moment []
        scalar dipm_;

        //- Pitzer's accentric factor []
        scalar omega_;

        //- Solubility parameter [(J/m^3)^0.5]
        scalar delta_;


public:

    TypeName("liquid");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            liquid,
            ,
            (),
            ()
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            liquid,
            Istream,
            (Istream& is),
            (is)
        );


    // Constructors

        //- Construct from components
        liquid
        (
            scalar W,
            scalar Tc,
            scalar Pc,
            scalar Vc,
            scalar Zc,
            scalar Tt,
            scalar Pt,
            scalar Tb,
            scalar dipm,
            scalar omega,
            scalar delta
        )
        :
            W_(W),
            Tc_(Tc),
            Pc_(Pc),
            Vc_(Vc),
            Zc_(Zc),
            Tt_(Tt),
            Pt_(Pt),
            Tb_(Tb),
            dipm_(dipm),
            omega_(omega),
            delta_(delta)
        {}

        //- Construct from Istream
        liquid(Istream& is)
        :
            W_(readScalar(is)),
            Tc_(readScalar(is)),
            Pc_(readScalar(is)),
            Vc_(readScalar(is)),
            dipm_(readScalar(is)),
            omega_(readScalar(is)),
            delta_(readScalar(is))
        {}

        //- Return a pointer to a new liquid created from input
        static autoPtr<liquid> New(Istream& is);


    //- Destructor
    virtual ~liquid()
    {}


    // Member Functions

        // Physical constants which define the specie

            //- Molecular weight [kg/kmol]
            inline scalar W() const;

            //- Critical temperature [K]
            inline scalar Tc() const;

            //- Critical pressure [Pa]
            inline scalar Pc() const;

            //- Critical volume [m^3/mol]
            inline scalar Vc() const;

            //- Critical compressibilty factor
            inline scalar Zc() const;

            //- Triple point temperature [K]
            inline scalar Tt() const;

            //- Triple point pressure [Pa]
            inline scalar Pt() const;

            //- Normal boiling temperature [K]
            inline scalar Tb() const;

            //- Dipole moment []
            inline scalar dipm() const;

            //- Pitzer's ascentric factor []
            inline scalar omega() const;

            //- Solubility parameter [(J/m^3)^(1/2)]
            inline scalar delta() const;


        // Physical property pure virtual functions

            //- Liquid rho [kg/m^3]
            virtual scalar rho(scalar p, scalar T) const = 0;

            //- Vapour pressure [Pa]
            virtual scalar pv(scalar p, scalar T) const = 0;

            //- Heat of vapourisation [J/kg]
            virtual scalar hl(scalar p, scalar T) const = 0;

            //- Liquid heat capacity [J/(kg K)]
            virtual scalar cp(scalar p, scalar T) const = 0;

            //- Liquid h [J/kg]
            virtual scalar h(scalar p, scalar T) const = 0;

            //- Ideal gas heat capacity [J/(kg K)]
            virtual scalar cpg(scalar p, scalar T) const = 0;

            //- Liquid viscosity [Pa s]
            virtual scalar mu(scalar p, scalar T) const = 0;

            //- Vapour viscosity [Pa s]
            virtual scalar mug(scalar p, scalar T) const = 0;

            //- Liquid thermal conductivity  [W/(m K)]
            virtual scalar K(scalar p, scalar T) const = 0;

            //- Vapour thermal conductivity  [W/(m K)]
            virtual scalar Kg(scalar p, scalar T) const = 0;

            //- Surface tension [N/m]
            virtual scalar sigma(scalar p, scalar T) const = 0;

            //- Vapour diffussivity [m2/s]
            virtual scalar D(scalar p, scalar T) const = 0;


    // I-O

        //- Write the function coefficients
        virtual void writeData(Ostream& os) const
        {
            os  << W_ << token::SPACE
                << Tc_ << token::SPACE
                << Pc_ << token::SPACE
                << Vc_ << token::SPACE
                << dipm_ << token::SPACE
                << omega_<< token::SPACE
                << delta_;
        }

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const liquid& l)
        {
            l.writeData(os);
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "liquidI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
