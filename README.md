# WRPF

Program written in Python 3.9.10 to calculate response properties in different approximations of principal propagator in No--Relativistic level

* RPA

Now, only static principal propagator and, lineal and quadratic responses

LRESC and PPESC for shielding is implemented. LRESC have two version, one version is with some constants by scale factor, instead another version is without scale factor.

# Installation

Build the venv

python3.9 -m venv pyint

# Integrals

Integrals implemented using point nucleu, cartessian primitives and until i shell. The transformation to one or two--spherical also i shell. Also, only is considered real primitive functions (gaussian).

* One--body
    * Potenial 
    * Kinetic
    * Angular momentum
    * Spin dipolar
    * Fermi contact
    * Darwin
    * Mass--velocity
    * Nuclear eletric field gradient
    * Dipole length
    * Dipole velocity
    * Paramagnetic spin-orbit
    * Diamagnetic nuclear shielding tensor
    * Kinetic-energy correction to the diamagnetic contribution to nuclear shielding
    * Kinetic-energy correction to the paramagnetic spin-orbit
    * Orbital-zeeman correction to the paramagnetic spin-orbit 
    * Kinetic energy correction to the orbital zeeman
    * Spin-orbit
    * Laplacian: Dxx, Dyy, Dzz, Dxy, Dxz, and Dzz
    * External magnetic-field dependence of the spin–orbit operator integrals
    * PA²P

* Two--body
    * Electron repulsion

PPotential integrals use boys function.

# Fock Matriz 

* Fock--Matrix Implemented
    * Hartree--Fock
    * Add Massvelocity and Darwin to Hartree--Fock energies 
  
