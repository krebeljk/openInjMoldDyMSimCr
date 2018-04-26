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

#include "mojCrHeRhoThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::mojCrHeRhoThermo<BasicPsiThermo, MixtureType>::calculate()
{
    //const scalarField& hCells = this->he().internalField(); //Kristjan: governed by TEqn
    const scalarField& pCells = this->p_.internalField();
    const scalarField& strigCells = this->strig_.internalField();

    scalarField& TCells = this->T_.internalField();
    scalarField& psiCells = this->psi_.internalField();
    scalarField& rhoCells = this->rho_.internalField();
    scalarField& muCells = this->mu_.internalField();
    scalarField& alphaCells = this->alpha_.internalField();
    scalarField& crCells = this->cr_.internalField();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        //TCells[celli] = mixture_.THE //Kristjan: governed by TEqn
        //(
        //    hCells[celli],
        //    pCells[celli],
        //    TCells[celli]
        //);

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli], crCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli], crCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli], strigCells[celli]);

        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pstrig = this->strig_.boundaryField()[patchi];
        fvPatchScalarField& pcr = this->cr_.boundaryField()[patchi];
        //fvPatchScalarField& pU = strig.boundaryField()[patchi]; // tukaj boundary

        fvPatchScalarField& ppsi = this->psi_.boundaryField()[patchi];
        fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];

        //fvPatchScalarField& ph = this->he().boundaryField()[patchi];//Kristjan: governed by TEqn


        fvPatchScalarField& pmu = this->mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                //ph[facei] = mixture_.HE(pp[facei], pT[facei]);//Kristjan: governed by TEqn


                ppsi[facei] = mixture_.psi(pp[facei], pT[facei], pcr[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei], pcr[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei], pstrig[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                //pT[facei] = mixture_.THE(ph[facei], pp[facei], pT[facei]);//Kristjan: governed by TEqn

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei], pcr[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei], pcr[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei], pstrig[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::mojCrHeRhoThermo<BasicPsiThermo, MixtureType>::mojCrHeRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    mojHeThermo<BasicPsiThermo, MixtureType>(mesh, phaseName)
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::mojCrHeRhoThermo<BasicPsiThermo, MixtureType>::~mojCrHeRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::mojCrHeRhoThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering mojCrHeRhoThermo<MixtureType>::correct()" << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "exiting mojCrHeRhoThermo<MixtureType>::correct()" << endl;
    }
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::mojHeThermo<BasicThermo, MixtureType>::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));
    scalarField& cp = tCp();

    forAll(T, facei)
    {
        cp[facei] =
            this->patchFaceMixture(patchi, facei).Cp(p[facei], T[facei]);
    }

    return tCp;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::mojHeThermo<BasicThermo, MixtureType>::Cp() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& cp = tCp();

    forAll(this->T_, celli)
    {
        cp[celli] =
            this->cellMixture(celli).Cp(this->p_[celli], this->T_[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pCp = cp.boundaryField()[patchi];

        forAll(pT, facei)
        {
            pCp[facei] =
                this->patchFaceMixture(patchi, facei).Cp(pp[facei], pT[facei]);
        }
    }

    return tCp;
}
// ************************************************************************* //
