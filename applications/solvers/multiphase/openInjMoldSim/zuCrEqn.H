{

    volScalarField zuGro
    (
       zuGmax * 
       exp( -2.0 *
           sqr(
              T - zuTref - zua * (p - zuPref)
           )/zuBeta
       ) 
    );
    volScalarField zuNuc
    (
        zuN1 * (T - twoSevenThree) + zuN2
    );

    fvScalarMatrix zuFi2eq
    (
        fvm::ddt(zuFi2)
      + fvm::div(phi, zuFi2)
      - 8.0 * constant::mathematical::pi
        * zuGro * zuNuc
    );

    zuFi2eq.relax();
    zuFi2eq.solve();

    fvScalarMatrix zuFi1eq
    (
        fvm::ddt(zuFi1)
      + fvm::div(phi, zuFi1)
      - zuFi2 * zuGro
    );

    zuFi1eq.relax();
    zuFi1eq.solve();

    fvScalarMatrix zuFi0eq
    (
        fvm::ddt(zuFi0)
      + fvm::div(phi, zuFi0)
      - zuFi1 * zuGro
    );

    zuFi0eq.relax();
    zuFi0eq.solve();
}
