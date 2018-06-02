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
    Foam::mojCrHeRhoThermo

Description
    Energy for a mixture based on density

SourceFiles
    mojCrHeRhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef mojCrHeRhoThermo_H
#define mojCrHeRhoThermo_H

#include "mojRhoThermo.H"
#include "mojHeThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class mojCrHeRhoThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class mojCrHeRhoThermo
:
    public mojHeThermo<BasicPsiThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        mojCrHeRhoThermo(const mojCrHeRhoThermo<BasicPsiThermo, MixtureType>&);


public:

    //- Runtime type information
    TypeName("mojCrHeRhoThermo");


    // Constructors

        //- Construct from mesh and phase name
        mojCrHeRhoThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    //- Destructor
    virtual ~mojCrHeRhoThermo();


    // Member functions

        //- Update properties
        virtual void correct();

        //- Heat capacity at constant pressure for patch [J/kg/K]
        virtual tmp<scalarField> Cp
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi
        ) const;

        //- Heat capacity at constant pressure [J/kg/K]
        virtual tmp<volScalarField> Cp() const;

        //- Latent heat [J/kg]
        virtual tmp<scalarField> Hfus
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi
        ) const;

        //- Latent heat [J/kg]
        virtual tmp<volScalarField> Hfus() const;

        //- thermal conductivity [J/s/K/m]
        virtual tmp<scalarField> kappa
        (
            const scalarField& p,
            const scalarField& cr,
            const label patchi
        ) const;

        //- thermal conductivity [J/s/K/m]
        virtual tmp<volScalarField> kappa() const;

        //- Vol. therm. exp. [1/K]
        virtual tmp<scalarField> cTExp
        (
            const scalarField& p,
            const scalarField& T,
            const scalarField& cr,
            const label patchi
        ) const;


        //- Vol. therm. exp. [1/K]
        virtual tmp<volScalarField> cTExp() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mojCrHeRhoThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
