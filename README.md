# openInjMoldDyMSimCr v1.0

This is a fork of [openInjMoldSim](https://github.com/krebeljk/openInjMoldSim) that implements the Kolmogorof-Avrami-Evans model of crystallization as Schenider's equations.
The material data provided with the case approximate HDPE behavior on cooling and compare the packing pressure evolution.

## Functionality
* The dynamic mesh functionality is used to model mold deformation.
* Pressure dependence of the thermal heat coefficient is modeled with a custom boundary condition [externalWallHeatFluxTemperatureP](https://github.com/krebeljk/externalWallHeatFluxTemperatureP) (DOI:10.5281/zenodo.4308349).
* The thermal conductivity depends on pressure and crystallinity.
* Specific heat and latent heat are provided as separate temperature dependent tables.
* Latent heat is released according to relative crystallinity progress.
* Elastic deviatoric stress is developed in the solid phase.

## Contact
krebeljk()gmail.com

## Acknowledgments

* The work was supported by the [Laboratory for Numerical Modelling and Simulation - LNMS](http://lab.fs.uni-lj.si/lnms/).

### Special thanks
* **Janez Turk** - Introduced the key modifications to the original OpenFOAM library.

## License

This project is licensed under the GPU License - see the [LICENSE.md](LICENSE.md) file for details.
