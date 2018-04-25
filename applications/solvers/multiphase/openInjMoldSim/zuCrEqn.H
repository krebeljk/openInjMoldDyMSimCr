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

    //Info<< "min(T) " << min(T).value() << endl;
}
