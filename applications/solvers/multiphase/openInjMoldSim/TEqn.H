{
    fvScalarMatrix TEqn
    (
         mixture.Cp()
       *(
            fvm::ddt(rho, T)
          + fvm::div(rhoPhi, T)
        )
      - fvm::laplacian(mixture.kappa(), T)
      - mixture.mu() * pow(strig,2)///viscous heating
      - T * mixture.cTExp()
        * (
            fvc::ddt(p)
          + (U & fvc::grad(p))
          )// Dp/Dt
      - rho1*alpha1
       *ultCryst*mixture.thermo1().Hfus()//Latent heat
        * (
            fvc::ddt(zuCr)
          + (U & fvc::grad(zuCr))
          )// Dcr/Dt
    );

    TEqn.relax();
    TEqn.solve();

    mixture.correct();

    Info<< "min(T) " << min(T).value() << endl;
}
