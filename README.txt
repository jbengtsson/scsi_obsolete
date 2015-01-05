scsi

Self-Consistent Symplectic Integrator for charged particle beam dynamics,
based on TPSA (Truncated Power Series Algebra),
aka PTC (Polymorphic Tracking Code); originated 1994.

Author: Johan Bengtsson

The Lex/Yacc based lattice parser is based on glps by Lingyan Yang.
The Python API was prototyped by James Rowland, DIAMOND, 2004.

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
