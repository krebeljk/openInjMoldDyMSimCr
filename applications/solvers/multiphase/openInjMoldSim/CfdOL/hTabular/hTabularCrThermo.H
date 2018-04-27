/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::hTabularCrThermo

Description
    Thermodynamics package templated on the equation of state, using tabular
    functions for cp, h and s

    h and s derived from cp

SourceFiles
    hTabularCrThermoI.H
    hTabularCrThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hTabularCrThermo_H
#define hTabularCrThermo_H

#include "scalar.H"
#include "interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState>
class hTabularCrThermo;

template<class EquationOfState>
inline hTabularCrThermo<EquationOfState> operator+
(
    const hTabularCrThermo<EquationOfState>&,
    const hTabularCrThermo<EquationOfState>&
);

template<class EquationOfState>
inline hTabularCrThermo<EquationOfState> operator-
(
    const hTabularCrThermo<EquationOfState>&,
    const hTabularCrThermo<EquationOfState>&
);

template<class EquationOfState>
inline hTabularCrThermo<EquationOfState> operator*
(
    const scalar,
    const hTabularCrThermo<EquationOfState>&
);

template<class EquationOfState>
inline hTabularCrThermo<EquationOfState> operator==
(
    const hTabularCrThermo<EquationOfState>&,
    const hTabularCrThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hTabularCrThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                      Class hTabularCrThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hTabularCrThermo
:
    public EquationOfState
{
    // Private data

        //- Heat of formation
        //  Note: input in [J/kg], but internally uses [J/kmol]
        scalar Hf_;
	
	interpolation2DTable<scalar> cpTable;
	interpolation2DTable<scalar> hfusTable;//latent heat


    // Private Member Functions

        //- Construct from components
        inline hTabularCrThermo
        (
            const EquationOfState& pt,
            const scalar Hf
        );


public:

    // Constructors

        //- Construct from Istream
        hTabularCrThermo(Istream& is);

        //- Construct from dictionary
        hTabularCrThermo(const dictionary& dict);

        //- Construct as copy
        inline hTabularCrThermo(const hTabularCrThermo&);

        //- Construct as a named copy
        inline hTabularCrThermo(const word&, const hTabularCrThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hTabularCr<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kmol K)]
            inline scalar cp(const scalar p, const scalar T) const;

            //- Heat of fusion [J/kg]
            inline scalar Hfus(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kmol]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kmol]
            inline scalar hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;
	    
            inline scalar s(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline hTabularCrThermo& operator=(const hTabularCrThermo&);
        inline void operator+=(const hTabularCrThermo&);
        inline void operator-=(const hTabularCrThermo&);
        inline void operator*=(const scalar);


    // Friend operators

        friend hTabularCrThermo operator+ <EquationOfState>
        (
            const hTabularCrThermo&,
            const hTabularCrThermo&
        );

        friend hTabularCrThermo operator- <EquationOfState>
        (
            const hTabularCrThermo&,
            const hTabularCrThermo&
        );

        friend hTabularCrThermo operator* <EquationOfState>
        (
            const scalar,
            const hTabularCrThermo&
        );

        friend hTabularCrThermo operator== <EquationOfState>
        (
            const hTabularCrThermo&,
            const hTabularCrThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hTabularCrThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hTabularCrThermoI.H"

#ifdef NoRepository
#   include "hTabularCrThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
