/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "mojRhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mojRhoThermo, 0);
    defineRunTimeSelectionTable(mojRhoThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mojRhoThermo::mojRhoThermo(const fvMesh& mesh, const word& phaseName)
:
    mojFluidThermo(mesh, phaseName),
    rho_
    (
        IOobject
        (
            phasePropertyName("thermo:rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )
    ,

    cr_
    (
        IOobject
        (
            phasePropertyName("thermo:cr"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("thermo:initcr", dimensionSet(0,0,0,0,0,0,0), 0.0)
    )
    ,

    crV_
    (
        IOobject
        (
            phasePropertyName("thermo:crV"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("thermo:initcrV", dimensionSet(0,0,0,0,0,0,0), 0.0)
    )

{}


Foam::mojRhoThermo::mojRhoThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    mojFluidThermo(mesh, dict, phaseName),
    rho_
    (
        IOobject
        (
            phasePropertyName("thermo:rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )

    ,
    cr_
    (
        IOobject
        (
            phasePropertyName("thermo:cr"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("thermo:initcr", dimensionSet(0,0,0,0,0,0,0), 0.0)
    )
    ,
    crV_
    (
        IOobject
        (
            phasePropertyName("thermo:crV"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("thermo:initcrV", dimensionSet(0,0,0,0,0,0,0), 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mojRhoThermo> Foam::mojRhoThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return mojBasicThermo::New<mojRhoThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mojRhoThermo::~mojRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::mojRhoThermo::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::mojRhoThermo::rho(const label patchi) const
{
    return rho_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::mojRhoThermo::rho()
{
    return rho_;
}


Foam::volScalarField& Foam::mojRhoThermo::cr()
{
    return cr_;
}
Foam::volScalarField& Foam::mojRhoThermo::crV()
{
    return crV_;
}

const Foam::volScalarField& Foam::mojRhoThermo::psi() const
{
    return psi_;
}

Foam::tmp<Foam::volScalarField> Foam::mojRhoThermo::Hfus() const
{
    const fvMesh& mesh = this->cr_.mesh();

    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Hfus",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
        dimensionedScalar("thermo:initHfus", dimEnergy/dimMass, 0.0)
        )
    );
    return tCp;
}

Foam::tmp<Foam::volScalarField> Foam::mojRhoThermo::mu() const
{
    return mu_;
}


Foam::tmp<Foam::scalarField> Foam::mojRhoThermo::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}


// ************************************************************************* //
