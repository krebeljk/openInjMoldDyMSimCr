{

    // growth rate
    volScalarField zuT = T  - zua * (p - zuPref);// pressure shifted

    volScalarField zuGro
    (
       zuG0
       * exp( - zuED /(zuR * zuT)
              - zuC3 * zuTm0G/( zuT * max((zuTm0G - zuT), zuDelTSafe))
            )
    );
    // nucleation
    volScalarField zuNuc
    (
        zuN0 * exp( (zuTm0N - max(zuT, zuTNStop))/zuDelT0 )//freeze below zuTNStop
        * pos(zuTm0N - zuT) //inhibit above zuTm0N with N=0
    );

    // fi2
    fvScalarMatrix zuFi2eq
    (
        fvm::ddt(zuFi2)
      + fvm::div(phi, zuFi2)
      - fvc::div(U) * zuFi2
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
      - fvc::div(U) * zuFi1
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
      - fvc::div(U) * zuFi0
      - zuFi1 * zuGro
        * pos(0.9999 - zuCr) //stop at cr=1
        * pos(alpha1 - 0.9) //only for poly
    );

    zuFi0eq.relax();
    zuFi0eq.solve();

    // zuCr
    zuCr = 1 - exp(-zuFi0);
    crPolyV = zuCr; //transfer vol-crystallinity to thermo
    zuCrM = crPolyM; //write cr mass fraction output

}
