scsi

Author: Johan Bengtsson

Self-Consistent Symplectic Integrator for charged particle beam dynamics,
based on TPSA (Truncated Power Series Algebra), aka PTC (Polymorphic Tracking
Code), originated 1994; by implementing a transparent polymorphic number object
with reference counting for FP/TPSA in C++.

The symplectic integrator for realistic modeling of magnetic lattices for
ring-based synchrotrons was initially implemented in Pascal, by the author,
with care taken for the software architecture and resulting records/modules
(-> "objects") to reflect the structure of the mathematical objects describing
the underlying beam dynamics model.


Numerical Recipes was replaced GSL by Piotr Goryl and Bartek Sulkowski,
SOLARIS, 2014.

The Lex/Yacc based lattice parser is an adaptation of glps implemented by
Lingyun Yang, ALS, 2009.

The Python API was prototyped by James Rowland, DIAMOND, 2004.

The symplectic integrator for RADIA kick maps was implemented by Laurent
Nadolski, SOLEIL, 2002.

The original Pascal library/code was machine translated to C (with p2c) by
Michael Boege, SLS, 1998.


Requirements:

   GNU autoconf and automake environments.
   GNU C/C++ and FORTRAN-95 compilers: gcc and gfortran.
   GNU Scientific Library GSL.


To install:

   mkdir git_repos
   cd git_repos
   git clone git@github.com:jbengtsson/scsi.git
   cd scsi
   ./make_scsi.sh
