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

 scalar n_, Tau_, D1_, D2_, D3_, D3s_, A1_, A2_;

\*---------------------------------------------------------------------------*/

#include "crossWLFTransportCr.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::crossWLFTransportCr<Thermo>::crossWLFTransportCr(Istream& is)
:
    Thermo(is),
    n_(readScalar(is)),
    Tau_(readScalar(is)),
    D1_(readScalar(is)),
    D2_(readScalar(is)),
    D3_(readScalar(is)),
    D3s_(readScalar(is)),
    A1_(readScalar(is)),
    A2_(readScalar(is)),
    etaMin_(readScalar(is)),
    etaMax_(readScalar(is)),
    TnoFlow_(readScalar(is)),
    deltaTempInterp_(readScalar(is)),
    crF_(readScalar(is)),
    crH_(readScalar(is)),
    crM_(readScalar(is)),
    crUlt_(readScalar(is))
{
    is.check("crossWLFTransportCr<Thermo>::crossWLFTransportCr(Istream&)");
    kappa_ = interpolation2DTable<scalar>("constant/kappaTableCr");
    kappa_.outOfBounds(interpolation2DTable<scalar>::CLAMP);
}


template<class Thermo>
Foam::crossWLFTransportCr<Thermo>::crossWLFTransportCr(const dictionary& dict)
:
    Thermo(dict),
    n_(readScalar(dict.subDict("transport").lookup("n"))),
    Tau_(readScalar(dict.subDict("transport").lookup("Tau"))),
    D1_(readScalar(dict.subDict("transport").lookup("D1"))),
    D2_(readScalar(dict.subDict("transport").lookup("D2"))),
    D3_(readScalar(dict.subDict("transport").lookup("D3"))),
    D3s_(dict.subDict("transport").lookupOrDefault<scalar>("D3s"
    , readScalar(dict.subDict("transport").lookup("D3")))),
    A1_(readScalar(dict.subDict("transport").lookup("A1"))),
    A2_(readScalar(dict.subDict("transport").lookup("A2"))),
    etaMin_(readScalar(dict.subDict("transport").lookup("etaMin"))),
    etaMax_(readScalar(dict.subDict("transport").lookup("etaMax"))),
    TnoFlow_(readScalar(dict.subDict("transport").lookup("TnoFlow"))),
    deltaTempInterp_(dict.subDict("transport").lookupOrDefault<scalar>("deltaTempInterp", 5.0)),
    crF_(dict.subDict("transport").lookupOrDefault<scalar>("crF", 1e3)),
    crH_(dict.subDict("transport").lookupOrDefault<scalar>("crH", 0.2)),
    crM_(dict.subDict("transport").lookupOrDefault<scalar>("crM", 2.0)),
    crUlt_(dict.subDict("transport").lookupOrDefault<scalar>("crUlt", 0.66))
{
    kappa_ = interpolation2DTable<scalar>("constant/kappaTableCr");
    kappa_.outOfBounds(interpolation2DTable<scalar>::CLAMP);
    Info << "CrossWLF:" << endl;
    Info << "n                : " << n_              << endl;
    Info << "Tau              : " << Tau_            << endl;
    Info << "D1               : " << D1_             << endl;
    Info << "D2               : " << D2_             << endl;
    Info << "D3               : " << D3_             << endl;
    Info << "D3s              : " << D3s_            << endl;
    Info << "A1               : " << A1_             << endl;
    Info << "A2               : " << A2_             << endl;
    Info << "Tab. Therm. Cond.  "                    << endl;
    Info << "etaMin           : " << etaMin_         << endl;
    Info << "etaMax           : " << etaMax_         << endl;
    Info << "TnoFlow          : " << TnoFlow_        << endl;
    Info << "deltaTempInterp  : " << deltaTempInterp_<< endl;
    Info << "crF              : " << crF_            << endl;
    Info << "crH              : " << crH_            << endl;
    Info << "crM              : " << crM_            << endl;
    Info << "crUlt            : " << crUlt_          << endl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::crossWLFTransportCr<Thermo>::write(Ostream& os) const
{
    os  << this->specie::name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("n", n_);
    dict.add("Tau", Tau_);
    dict.add("D1", D1_);
    dict.add("D2", D2_);
    dict.add("D3", D3_);
    dict.add("D3s", D3s_);
    dict.add("A1", A1_);
    dict.add("A2", A2_);
    dict.add("etaMin", etaMin_);
    dict.add("etaMax", etaMax_);
    dict.add("TnoFlow", TnoFlow_);
    dict.add("deltaTempInterp", deltaTempInterp_);
    dict.add("crF", crF_);
    dict.add("crH", crH_);
    dict.add("crM", crM_);
    dict.add("crUlt", crUlt_);
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const crossWLFTransportCr<Thermo>& st
)
{
    os << static_cast<const Thermo&>(st)
    << tab << st.n_
    << tab << st.Tau_
    << tab << st.D1_
    << tab << st.D2_
    << tab << st.D3_
    << tab << st.D3s_
    << tab << st.A1_
    << tab << st.A2_
    << tab << st.etaMin_
    << tab << st.etaMax_
    << tab << st.TnoFlow_
    << tab << st.deltaTempInterp_ << endl
    << tab << st.crF_ << endl
    << tab << st.crH_ << endl
    << tab << st.crM_ << endl
    << tab << st.crUlt_ << endl;

    os.check
    (
        "Ostream& operator<<(Ostream&, const crossWLFTransportCr<Thermo>&)"
    );

    return os;
}


// ************************************************************************* //
