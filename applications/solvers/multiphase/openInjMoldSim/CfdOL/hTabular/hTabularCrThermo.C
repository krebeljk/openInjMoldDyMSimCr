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

\*---------------------------------------------------------------------------*/

#include "hTabularCrThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::hTabularCrThermo<EquationOfState>::hTabularCrThermo
(
    Istream& is
)
:
    EquationOfState(is),
    Hf_(readScalar(is))
{
    Hf_ *= this->W();
    cpTable = interpolation2DTable<scalar>("constant/cpTableCr");
    cpTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);
    hfusTable = interpolation2DTable<scalar>("constant/hfusTableCr");
    hfusTable .outOfBounds(interpolation2DTable<scalar>::CLAMP);
}


template<class EquationOfState>
Foam::hTabularCrThermo<EquationOfState>::hTabularCrThermo
(
    const dictionary& dict
)
:
    EquationOfState(dict),
    Hf_(readScalar(dict.subDict("thermodynamics").lookup("Hf")))
{
    Hf_ *= this->W();
    cpTable = interpolation2DTable<scalar>("constant/cpTableCr");
    cpTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);
    hfusTable = interpolation2DTable<scalar>("constant/hfusTableCr");
    hfusTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::hTabularCrThermo<EquationOfState>::write
(
    Ostream& os
) const
{
    EquationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("Hf", Hf_/this->W());
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hTabularCrThermo<EquationOfState>& pt
)
{
    os  << static_cast<const EquationOfState&>(pt) << tab
        << pt.Hf_/pt.W() << tab;

    os.check
    (
        "operator<<"
        "("
            "Ostream&, "
            "const hTabularCrThermo<EquationOfState>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
