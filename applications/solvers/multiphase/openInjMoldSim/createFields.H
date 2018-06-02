Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField shrRate
(
    IOobject
    (
        "shrRate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt(2.0*symm(fvc::grad(U))&&symm(fvc::grad(U)))
);

#include "createPhi.H"

Info<< "Constructing mojTwoPhaseMixtureThermo\n" << endl;
mojTwoPhaseMixtureThermo mixture(mesh);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

Info<< "Reading thermophysical properties\n" << endl;

volScalarField& p = mixture.p();
volScalarField& T = mixture.T();
volScalarField& strig = mixture.thermo1().strig();
volScalarField& crPoly = mixture.thermo1().cr();//reference to rel cry field of poly
volScalarField& rho1 = mixture.thermo1().rho();
const volScalarField& psi1 = mixture.thermo1().psi();
volScalarField& rho2 = mixture.thermo2().rho();
const volScalarField& psi2 = mixture.thermo2().psi();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);


dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    mixture
);

mesh.setFluxRequired(p_rgh.name());

//kristjan - za HTC

volScalarField mojKappaOut
(
    IOobject
    (
        "mojKappaOut",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mixture.kappa()
);

//kristjan - pregled viskoznosti

volScalarField visc
(
    IOobject
    (
        "visc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*mixture.thermo1().mu() + alpha2*mixture.thermo2().mu()
);

//kristjan - residuals p_rgh

volScalarField p_rgh_resid
(
    IOobject
    (
        "p_rgh_resid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("initOne", dimensionSet(0,0,-1,0,0,0,0), 1.0)
);

//kristjan - Calc Celsius
volScalarField Tc
(
    IOobject
    (
        "Tc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("initTc", dimensionSet(0,0,0,1,0,0,0), 0.0)
);
dimensionedScalar twoSevenThree("twoSevenThree", dimensionSet(0,0,0,1,0,0,0), 273.15);
Tc = T - twoSevenThree;

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha1 solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField dgdt
(
    pos(alpha2)*fvc::div(phi)/max(alpha2, scalar(0.0001))
);

// Construct interface from alpha1 distribution
interfaceProperties interface(alpha1, U, mixture);

// Construct compressible turbulence model
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New(rho, U, rhoPhi, mixture)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));



Info<< "Reading solidification properties\n" << endl;
    const word sldDictName("solidificationProperties");

    // Create and input-output object - this holds the path to the dict and its name
    IOobject sldDictIO
    (
        sldDictName, // name of the file
        mesh.time().constant(), // path to where the file is
        mesh, // reference to the mesh needed by the constructor
        IOobject::READ_IF_PRESENT // indicate that reading this dictionary is optional
    );

    // Check the if the dictionary is present and follows the OF format
    if (!sldDictIO.headerOk())
    {
        Info << "No constant/solidificationProperties dictionary found." << endl;
        Info << "Flow cessation modeled as visc(TnoFlow)=etaMax." << nl << nl << endl;
    }

    // Get access to the dictionary
    dictionary sldDict;
    // Initialise the dictionary object
    sldDict = IOdictionary(sldDictIO);

    // This template method needs to know the type of the variable and can provide
    // a default value if the entry is not found in the dictionary
    dimensionedScalar shrMod("shearModulus", dimensionSet(1,-1,-2,0,0,0,0),
        sldDict.lookupOrDefault<scalar>("shearModulus", 900e6)
    );
    dimensionedScalar shrRateLimEl("shrRateLimEl", dimensionSet(0,0,-1,0,0,0,0),
        sldDict.lookupOrDefault<scalar>("shrRateLimEl", 1.0)
    );
    dimensionedScalar viscLimEl("viscLimEl", dimensionSet(1,-1,-1,0,0,0,0),
        sldDict.lookupOrDefault<scalar>("viscLimEl", 1e6)
    );

    if (sldDictIO.headerOk())
    {
        // Summarise what's been read and print in the console
        Info << nl << "Solidification properties: " << nl
             << "Shear modulus: " << shrMod << nl
             << "Highest shear rate with elasticity: " << shrRateLimEl << nl
             << "Lowest viscosity with elasticity: " << viscLimEl << nl
             << nl << endl;
    }


    volSymmTensorField elSigDev
    (
        IOobject
        (
            "elSigDev",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        )
        , mesh
        , dimensionedSymmTensor("zero", dimensionSet(1,-1,-2,0,0,0,0), symmTensor::zero)
        , "zeroGradient"
    );

//Crystallinity
//Zuidema 2001
Info<< "Reading crystallization properties\n" << endl;
    const word cryDictName("crystallizationProperties");


    // Create and input-output object - this holds the path to the dict and its name
    IOobject cryDictIO
    (
        cryDictName, // name of the file
        mesh.time().constant(), // path to where the file is
        mesh, // reference to the mesh needed by the constructor
        IOobject::READ_IF_PRESENT // indicate that reading this dictionary is optional
    );

    // Check the if the dictionary is present and follows the OF format
    if (!cryDictIO.headerOk())
    {
        Info << "No constant/crystallizationProperties dictionary found." << endl;
        Info << "Defaulting to Zuidema 2001 parameters of crystallization." << nl << nl << endl;
    }

    // Get access to the dictionary
    dictionary cryDict;
    // Initialise the dictionary object
    cryDict = IOdictionary(cryDictIO);

    //Growth rate - params
    dimensionedScalar zuGmax("zuGmax", dimensionSet(0,1,-1,0,0,0,0),
        cryDict.lookupOrDefault<scalar>("Gmax", 8.1e-6)
    ); //m/s
    dimensionedScalar zuTref("zuTref", dimensionSet(0,0,0,1,0,0,0),
        cryDict.lookupOrDefault<scalar>("Tref", 356.8)
    ); //K
    dimensionedScalar zua("zua", dimensionSet(-1,1,2,1,0,0,0),
        cryDict.lookupOrDefault<scalar>("a", 5e-7)
    ); //K/Pa
    dimensionedScalar zuPref("zuPref", dimensionSet(1,-1,-2,0,0,0,0),
        cryDict.lookupOrDefault<scalar>("pRef", 1e5)
    ); //Pa
    dimensionedScalar zuBeta("zuBeta", dimensionSet(0,0,0,2,0,0,0),
        cryDict.lookupOrDefault<scalar>("beta", 1126.9)
    ); //K2
    //Nucleation - params
    dimensionedScalar zuN1("zuN1", dimensionSet(0,-3,0,-1,0,0,0),
        cryDict.lookupOrDefault<scalar>("N1", -2.6087e13)
    ); //K^-1 m^-3
    dimensionedScalar zuN2("zuN2", dimensionSet(0,-3,0,0,0,0,0),
        cryDict.lookupOrDefault<scalar>("N2", 6.5783e15)
    ); //m^-3

    if (cryDictIO.headerOk())
    {
        // Summarise what's been read and print in the console
        Info << nl << "Crystallization properties: " << nl
             << "Gmax: " << zuGmax << nl
             << "Tref: " << zuTref << nl
             << "a: " << zua << nl
             << "pRef: " << zuPref << nl
             << "beta: " << zuBeta << nl
             << "N1: " << zuN1 << nl
             << "N2: " << zuN2 << nl
             << nl << endl;
    }

//Crystallization state fields
volScalarField zuFi2
(
    IOobject
    (
        "zuFi2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
    ,mesh
    ,dimensionedScalar("zuFi2init", dimensionSet(0,-2,0,0,0,0,0), 0.0)
    ,"zeroGradient"
);
volScalarField zuFi1
(
    IOobject
    (
        "zuFi1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
    ,mesh
    ,dimensionedScalar("zuFi1init", dimensionSet(0,-1,0,0,0,0,0), 0.0)
    ,"zeroGradient"
);
volScalarField zuFi0
(
    IOobject
    (
        "zuFi0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
    ,mesh
    ,dimensionedScalar("zuFi0init", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ,"zeroGradient"
);
volScalarField zuCr
(
    IOobject
    (
        "zuCr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    )
    ,mesh
    ,dimensionedScalar("zuCr", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ,"zeroGradient"
);


Info<< "Reading crystPropDict properties\n" << endl;
    const word cryDictName("crystPropDict");

    // Create and input-output object - this holds the path to the dict and its name
    IOobject cryDictIO
    (
        cryDictName, // name of the file
        mesh.time().constant(), // path to where the file is
        mesh, // reference to the mesh needed by the constructor
        IOobject::READ_IF_PRESENT // indicate that reading this dictionary is optional
    );

    // Check the if the dictionary is present and follows the OF format
    if (!cryDictIO.headerOk())
    {
        Info << "No constant/crystPropDict dictionary found." << endl;
        Info << "Ultimate crystallinity default value." << endl;
    }

    // Get access to the dictionary
    dictionary cryDict;
    // Initialise the dictionary object
    cryDict = IOdictionary(cryDictIO);

    // This template method needs to know the type of the variable and can provide
    // a default value if the entry is not found in the dictionary
    dimensionedScalar ultCryst("ultCryst", dimensionSet(0,0,0,0,0,0,0),
        cryDict.lookupOrDefault<scalar>("ultCryst", 0.66)
    );

    //if (cryDictIO.headerOk())
    // {
        // Summarise what's been read and print in the console
        Info << "Ultimate crystallinity: " << ultCryst << nl
             << nl << endl;
    //}
