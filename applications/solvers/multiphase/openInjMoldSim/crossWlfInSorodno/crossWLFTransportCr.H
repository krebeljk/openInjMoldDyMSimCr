/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::crossWLFTransportCr

Description
    Transport package using crossWLF's formula.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    \f[
        \mu = A_s \frac{\sqrt{T}}{1 + T_s / T}
    \f]

SourceFiles
    crossWLFTransportCrI.H
    crossWLFTransportCr.C

\*---------------------------------------------------------------------------*/

#ifndef crossWLFTransportCr_H
#define crossWLFTransportCr_H
#include "interpolation2DTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class crossWLFTransportCr;

template<class Thermo>
inline crossWLFTransportCr<Thermo> operator+
(
    const crossWLFTransportCr<Thermo>&,
    const crossWLFTransportCr<Thermo>&
);

template<class Thermo>
inline crossWLFTransportCr<Thermo> operator-
(
    const crossWLFTransportCr<Thermo>&,
    const crossWLFTransportCr<Thermo>&
);

template<class Thermo>
inline crossWLFTransportCr<Thermo> operator*
(
    const scalar,
    const crossWLFTransportCr<Thermo>&
);

template<class Thermo>
inline crossWLFTransportCr<Thermo> operator==
(
    const crossWLFTransportCr<Thermo>&,
    const crossWLFTransportCr<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const crossWLFTransportCr<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class crossWLFTransportCr Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class crossWLFTransportCr
:
    public Thermo
{
    // Private data

        // crossWLF's coefficients
        scalar n_, Tau_, D1_, D2_, D3_, D3s_, A1_, A2_;
        scalar etaMin_, etaMax_, TnoFlow_, deltaTempInterp_;
        // Titomanlio 1997 crystalline solidification params
        scalar crF_, crH_, crM_, crUlt_;

	// thermal conductivity [W/(Km)]
        interpolation2DTable<scalar> kappa_;


    // Private Member Functions

        //- Calculate the crossWLF coefficients
        //  given two viscosities and temperatures
        inline void calcCoeffs
        (
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );


public:

    // Constructors

        //- Construct from components
        inline crossWLFTransportCr
        (
            const Thermo& t,
            const scalar n_,
            const scalar Tau_,
            const scalar D1_,
            const scalar D2_,
            const scalar D3_,
            const scalar D3s_,
            const scalar A1_,
            const scalar A2_,
            const interpolation2DTable<scalar> kappa_,
            const scalar etaMin_,
            const scalar etaMax_,
            const scalar TnoFlow_,
            const scalar deltaTempInterp_,
            const scalar crF_,
            const scalar crH_,
            const scalar crM_,
            const scalar crUlt_
        );

        //- Construct from two viscosities
        inline crossWLFTransportCr
        (
            const Thermo& t,
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Construct as named copy
        inline crossWLFTransportCr(const word&, const crossWLFTransportCr&);

        //- Construct from Istream
        crossWLFTransportCr(Istream&);

        //- Construct from dictionary
        crossWLFTransportCr(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<crossWLFTransportCr> clone() const;

        // Selector from Istream
        inline static autoPtr<crossWLFTransportCr> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<crossWLFTransportCr> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "crossWLFCr<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T, const scalar strig) const;

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T, const scalar strig, const scalar cr) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar cr) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline crossWLFTransportCr& operator=(const crossWLFTransportCr&);

        inline void operator+=(const crossWLFTransportCr&);

        inline void operator-=(const crossWLFTransportCr&);

        inline void operator*=(const scalar);


    // Friend operators

        friend crossWLFTransportCr operator+ <Thermo>
        (
            const crossWLFTransportCr&,
            const crossWLFTransportCr&
        );

        friend crossWLFTransportCr operator- <Thermo>
        (
            const crossWLFTransportCr&,
            const crossWLFTransportCr&
        );

        friend crossWLFTransportCr operator* <Thermo>
        (
            const scalar,
            const crossWLFTransportCr&
        );

        friend crossWLFTransportCr operator== <Thermo>
        (
            const crossWLFTransportCr&,
            const crossWLFTransportCr&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const crossWLFTransportCr&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "crossWLFTransportCrI.H"

#ifdef NoRepository
#   include "crossWLFTransportCr.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
