{

    // growth rate
    volScalarField zuGro
    (
       zuGmax * 
       exp( -2.0 *
           sqr(
              T - zuTref - zua * (p - zuPref)
           )/zuBeta
       ) 
    );
    // nucleation
    volScalarField zuNuc
    (
        zuN1 * (T - twoSevenThree) + zuN2
    );

    // fi2
    fvScalarMatrix zuFi2eq
    (
        fvm::ddt(zuFi2)
      + fvm::div(phi, zuFi2)
      - 8.0 * constant::mathematical::pi
        * zuGro * zuNuc
        * pos(0.9999 - zuCr) //stop at cr=1
        * pos(alpha1 - 0.9) //only for poly
    );

    zuFi2eq.relax();
    zuFi2eq.solve();

    // fi1
    fvScalarMatrix zuFi1eq
    (
        fvm::ddt(zuFi1)
      + fvm::div(phi, zuFi1)
      - zuFi2 * zuGro
        * pos(0.9999 - zuCr) //stop at cr=1
        * pos(alpha1 - 0.9) //only for poly
    );

    zuFi1eq.relax();
    zuFi1eq.solve();

    // fi0
    fvScalarMatrix zuFi0eq
    (
        fvm::ddt(zuFi0)
      + fvm::div(phi, zuFi0)
      - zuFi1 * zuGro
        * pos(0.9999 - zuCr) //stop at cr=1
        * pos(alpha1 - 0.9) //only for poly
    );

    zuFi0eq.relax();
    zuFi0eq.solve();

    // zuCr
    zuCr = 1 - exp(-zuFi0);
    crPoly = zuCr; //transfer crystallinity to thermo

}
