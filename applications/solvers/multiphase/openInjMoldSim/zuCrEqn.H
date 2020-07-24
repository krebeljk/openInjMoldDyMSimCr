{

    // growth rate
    volScalarField zuT = T  - zua * (p - zuPref);// pressure shifted

    volScalarField smoothStop1 = -(zuCr-0.999)/0.02;//ramp down final cry rate before cr=1 (that cant be reached)
    smoothStop1 = min(max(smoothStop1,0.0),1.0);
    volScalarField smoothStop2 = (alpha1-0.9)/0.03;//only for poly
    smoothStop2 = min(max(smoothStop2,0.0),1.0);

    volScalarField zuGro
    (
       zuG0
       * exp( - zuED /(zuR * zuT)
              - zuC3 * zuTm0G/( zuT * max((zuTm0G - zuT), zuDelTSafe))
              - zuGretard * zuCr //empirical growth impingement // Mandelkern discussed //Skg_03_05_E_07_09_03
            )
        * smoothStop1 //stop at cr=1
        * smoothStop2  //only for poly
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
      + fvm::SuSp(-fvc::div(phi), zuFi2)
      - 8.0 * constant::mathematical::pi
        * zuGro * zuNuc
    );

    zuFi2eq.relax();
    zuFi2eq.solve();

    // fi1
    fvScalarMatrix zuFi1eq
    (
        fvm::ddt(zuFi1)
      + fvm::div(phi, zuFi1)
      + fvm::SuSp(-fvc::div(phi), zuFi1)
      - zuFi2 * zuGro
    );

    zuFi1eq.relax();
    zuFi1eq.solve();

    // fi0
    fvScalarMatrix zuFi0eq
    (
        fvm::ddt(zuFi0)
      + fvm::div(phi, zuFi0)
      + fvm::SuSp(-fvc::div(phi), zuFi0)
      - zuFi1 * zuGro
    );

    zuFi0eq.relax();
    zuFi0eq.solve();

    // zuCr
    zuCr = 1 - exp(-zuFi0);
    crPolyV = zuCr; //transfer vol-crystallinity to thermo
    zuCrM = crPolyM; //write cr mass fraction output

}
