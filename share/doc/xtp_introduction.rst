Introduction
############

VOTCA-XTP is a versatile software package that can perform many tasks related to
electronic-structure calculations and materials modeling. One of its main
features is the possibility to perform excited state calculations on organic
molecules with the many-body Green's functions theory (GW-BSE) in a polarizable
environment (QM/MM). Both local and charge-transfer excitations can be modelled
using this technique. Furthermore XTP provides a complete workflow to go from an
MD topology to a full excited state transport simulation.

What VOTCA-XTP Can Do
---------------------

Ground- and Excited-State Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Ground-state calculations, either using the build in DFT package or the `ORCA <https://sites.google.com/site/orcainputlibrary/home>`_ package.
* Excited-state calculations, using the many-body Green's function method and Bethe Salpeter Equation (GW-BSE).
* Ground- and excited-state calculations in a polarizable environment (QM/MM).

Charge Transport Simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Convert a GROMACS or LAMMPS topology to a system of conjugated segments. 
* Replace distorted geometries by their ground-state geometries and map multipole data to atomic expansion sites.
* Compute site and reorganization energies.
* Compute pair data, such as transfer integrals and transfer rates.
* Perform kinetic Monte Carlo simulations.


