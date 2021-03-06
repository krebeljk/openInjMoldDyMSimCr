/*---------------------------------------------------------------------------*\
modified from: OpenFOAM-3.0.0/src/thermophysicalModels/specie/transport/sutherland/sutherlandTransport.H
 =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
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
    Foam::taitCr

Description
    Perfect gas equation of state.

SourceFiles
    taitCrI.H
    taitCr.C

\*---------------------------------------------------------------------------*/

#ifndef taitCr_H
#define taitCr_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class taitCr;

template<class Specie>
inline taitCr<Specie> operator+
(
    const taitCr<Specie>&,
    const taitCr<Specie>&
);

template<class Specie>
inline taitCr<Specie> operator-
(
    const taitCr<Specie>&,
    const taitCr<Specie>&
);

template<class Specie>
inline taitCr<Specie> operator*
(
    const scalar,
    const taitCr<Specie>&
);

template<class Specie>
inline taitCr<Specie> operator==
(
    const taitCr<Specie>&,
    const taitCr<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const taitCr<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class taitCr Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class taitCr
:
    public Specie
{
    // Private data

        //- Fluid constant
        scalar b1m_;
        scalar b2m_;
        scalar b3m_;
        scalar b4m_;
        scalar b1s_;
        scalar b2s_;
        scalar b3s_;
        scalar b4s_;
        scalar b5_;
        scalar b6_;
        scalar b10_;
        static scalar C_;// definicija v .C // univerzalna konstanta v Taitovi enacbi

public:

    // Constructors

        //- Construct from components
        inline taitCr
        (
            const Specie& sp,
            const scalar b1m,
            const scalar b2m,
            const scalar b3m,
            const scalar b4m,
            const scalar b1s,
            const scalar b2s,
            const scalar b3s,
            const scalar b4s,
            const scalar b5,
            const scalar b6,
            const scalar b10
        );

        //- Construct from Istream
        taitCr(Istream&);

        //- Construct from dictionary
        taitCr(const dictionary& dict);

        //- Construct as named copy
        inline taitCr(const word& name, const taitCr&);

        //- Construct and return a clone
        inline autoPtr<taitCr> clone() const;

        // Selector from Istream
        inline static autoPtr<taitCr> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<taitCr> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "taitCr<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return fluid constant [J/(kg K)]
            inline scalar R() const;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return density from crystallinity [kg/m^3]
            inline scalar rho(scalar p, scalar T, scalar cr) const;

            //- Return entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T, scalar cr) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (cp - cv) [J/(kmol K]
            inline scalar cpMcv(scalar p, scalar T) const;

            //- Return coef. of vol. therm. expansion [1/K]
            inline scalar cTExp(scalar p, scalar T) const;
            inline scalar cTExp(scalar p, scalar T, scalar cr) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const taitCr&);
        inline void operator-=(const taitCr&);

        inline void operator*=(const scalar);


    // Friend operators

        friend taitCr operator+ <Specie>
        (
            const taitCr&,
            const taitCr&
        );

        friend taitCr operator- <Specie>
        (
            const taitCr&,
            const taitCr&
        );

        friend taitCr operator* <Specie>
        (
            const scalar s,
            const taitCr&
        );

        friend taitCr operator== <Specie>
        (
            const taitCr&,
            const taitCr&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const taitCr&
        );

    // Izracun Entropije
private:

            //- Solidification interpolation to continuity
            inline scalar interpolation
            (
                 scalar vS
                ,scalar vM
                ,scalar p
                ,scalar T
            ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "taitCrI.H"

#ifdef NoRepository
#   include "taitCr.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
