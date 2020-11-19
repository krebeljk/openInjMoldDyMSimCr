#include "readTimeControls.H"

correctPhi = pimple.dict().lookupOrDefault<Switch>("correctPhi", true);

checkMeshCourantNo =
    pimple.dict().lookupOrDefault<Switch>("checkMeshCourantNo", false);

scalar maxSolidCo =
    runTime.controlDict().lookupOrDefault<scalar>("maxSolidCo", 1.0); // Kristjan
scalar pAuxRlx =
    runTime.controlDict().lookupOrDefault<scalar>("pAuxRlx", 0.1); // Kristjan
