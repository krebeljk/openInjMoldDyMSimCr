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

#include "createPhi.H"

Info<< "Constructing mojTwoPhaseMixtureThermo\n" << endl;
mojTwoPhaseMixtureThermo mixture(mesh);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

Info<< "Reading thermophysical properties\n" << endl;

volScalarField& p = mixture.p();
volScalarField& T = mixture.T();
volScalarField& shrRate = mixture.shrRate();
volScalarField& crPolyV = mixture.thermo1().crV();//reference to rel cry field of poly
volScalarField& crPolyM = mixture.thermo1().cr();//reference to rel cry field of poly
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

//kristjan - Calc Celsius
dimensionedScalar twoSevenThree("twoSevenThree", dimensionSet(0,0,0,1,0,0,0), 273.15);
volScalarField Tc
(
    IOobject
    (
        "Tc",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    T - twoSevenThree
);

//kristjan - viscosity

volScalarField visc
(
    IOobject
    (
        "visc",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*mixture.thermo1().mu() + alpha2*mixture.thermo2().mu()
);

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


volScalarField dgdtw
(
    "dgdtw",
    dgdt
);
dgdtw.write();

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
    dimensionedScalar viscLimEl("viscLimEl", dimensionSet(1,-1,-1,0,0,0,0),
        sldDict.lookupOrDefault<scalar>("viscLimEl", 1e6)
    );

volScalarField solid
(
    IOobject
    (
        "solid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("initSolid", dimless, 0.0)
);
    if (sldDictIO.headerOk())
    {
        // Summarise what's been read and print in the console
        Info << nl << "Solidification properties: " << nl
             << "Shear modulus: " << shrMod << nl
             << "Lowest viscosity with elasticity: " << viscLimEl << nl
             << nl << endl;

        volScalarField sigSmStop1 = (visc-viscLimEl)/(5e-2*viscLimEl);
        sigSmStop1 = min(max(sigSmStop1 ,0.0),1.0);
        solid = sigSmStop1;
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
        Info << "Defaulting to HDPE parameters of crystallization." << nl << nl << endl;
    }

    // Get access to the dictionary
    dictionary cryDict;
    // Initialise the dictionary object
    cryDict = IOdictionary(cryDictIO);

    //Growth rate - params
    dimensionedScalar zuG0("zuG0", dimensionSet(0,1,-1,0,0,0,0),
        cryDict.lookupOrDefault<scalar>("G0", 1e3)
    ); //m/s
    dimensionedScalar zua("zua", dimensionSet(-1,1,2,1,0,0,0),
        cryDict.lookupOrDefault<scalar>("a", 1.543e-7)
    ); //K/Pa
    dimensionedScalar zuPref("zuPref", dimensionSet(1,-1,-2,0,0,0,0),
        cryDict.lookupOrDefault<scalar>("pRef", 1e5)
    ); //Pa
    dimensionedScalar zuED("zuED", dimensionSet(1,2,-2,0,-1,0,0),
        cryDict.lookupOrDefault<scalar>("ED", 29.3e3)
    ); //J/(mol K)
    dimensionedScalar zuC3("zuC3", dimensionSet(0,0,0,1,0,0,0),
        cryDict.lookupOrDefault<scalar>("C3", 265)
    ); //K
    dimensionedScalar zuTm0G("zuTm0G", dimensionSet(0,0,0,1,0,0,0),
        cryDict.lookupOrDefault<scalar>("Tm0G", 419)
    ); //K
    dimensionedScalar zuGretard("zuGretard", dimensionSet(0,0,0,0,0,0,0),
        cryDict.lookupOrDefault<scalar>("Gretard", 0.0)
    ); //K
    dimensionedScalar zuR("zuR", dimensionSet(1,2,-2,-1,-1,0,0), 8.3144598); //J/(mol K) gas constant (should be taken form OpenFOAM)
    dimensionedScalar zuDelTSafe("zuDelTSafe", dimensionSet(0,0,0,1,0,0,0), 5); //K ... divide with this instead of 0

    //Nucleation - params
    dimensionedScalar zuN0("zuN0", dimensionSet(0,-3,0,0,0,0,0),
        cryDict.lookupOrDefault<scalar>("N0", 1e8)
    ); //m^-3
    dimensionedScalar zuTm0N("zuTm0N", dimensionSet(0,0,0,1,0,0,0),
        cryDict.lookupOrDefault<scalar>("Tm0N", 419)
    ); //K
    dimensionedScalar zuDelT0("zuDelT0", dimensionSet(0,0,0,1,0,0,0),
        cryDict.lookupOrDefault<scalar>("delT0", 2.5)
    ); //K
    dimensionedScalar zuTNStop("zuTNStop", dimensionSet(0,0,0,1,0,0,0),
        cryDict.lookupOrDefault<scalar>("TNStop", 333.15) // stop nucleation
    ); //K

    //Ultimate crystallinity
    dimensionedScalar ultCryst("ultCryst", dimensionSet(0,0,0,0,0,0,0),
        cryDict.lookupOrDefault<scalar>("ultCryst", 0.66)
    );//[1]

    if (cryDictIO.headerOk())
    {
        // Summarise what's been read and print in the console
        Info << nl << "Crystallization properties: " << nl
             << "G0: " << zuG0 << nl
             << "a: " << zua << nl
             << "pRef: " << zuPref << nl
             << "ED: " << zuED << nl
             << "C3: " << zuC3 << nl
             << "Tm0G: " << zuTm0G << nl
             << "Gretard: " << zuGretard << nl
             << "N0: " << zuN0 << nl
             << "Tm0N: " << zuTm0N << nl
             << "delT0: " << zuDelT0 << nl
             << "TNStop: " << zuTNStop << nl
             << "Ultimate crystallinity: " << ultCryst << nl
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
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
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
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
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
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    )
    ,mesh
    ,dimensionedScalar("zuFi0init", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ,"zeroGradient"
);
volScalarField zuCr // crystalline volume fraction
(
    IOobject
    (
        "zuCr",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    )
    ,mesh
    ,dimensionedScalar("zuCr", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ,"zeroGradient"
);

volScalarField zuCrM // crystalline mass fraction
(
    IOobject
    (
        "zuCrM",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    )
    ,mesh
    ,dimensionedScalar("zuCrM", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ,"zeroGradient"
);

crPolyV = zuCr; //transfer vol-crystallinity to thermo
mixture.correct(); //thermo calcs its mass-cr, rho, mu, kappa etc.

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


// pressure preventing material collapse bellow pMin
volScalarField pAux
(
    IOobject
    (
        "pAux",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeropAux", dimPressure, 0.0)
);
