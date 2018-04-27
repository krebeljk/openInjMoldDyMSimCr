{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - 1.0/mixture.Cp()*fvm::laplacian(mixture.kappa(), T)
      - mixture.mu() * pow(strig,2)/mixture.Cp()//viscous heating
      - T * mixture.cTExp()
        * (
            fvc::ddt(p)
          + fvc::div(phi, p)
          )// Dp/Dt
          /mixture.Cp()
      - rho * mixture.thermo1().Hfus()//Latent heat
        * (
            fvc::ddt(zuCr)
          + fvc::div(phi, zuCr)
          )// Dcr/Dt
          /mixture.Cp()
    );

    TEqn.relax();
    TEqn.solve();

    mixture.correct();

    Info<< "min(T) " << min(T).value() << endl;
}
