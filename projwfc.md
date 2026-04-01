3/27/26, 10:50 PM projwfc.x: input description
Input File Description
Program: projwfc.x / PWscf / Quantum ESPRESSO (version: 7.5)
TABLE OF CONTENTS
INTRODUCTION
&PROJWFC
prefix | outdir | ngauss | degauss | Emin | Emax | DeltaE | lsym | diag_basis | pawproj |
filpdos | filproj | lwrite_overlaps | lbinary_data | kresolveddos | tdosinboxes | n_proj_boxes |
irmin | irmax | plotboxes
Notes
Format of output files
Orbital Order
Defining boxes for the Local DOS(E)
Important notices
INTRODUCTION
Purpose of projwfc.x:
projects wavefunctions onto orthogonalized atomic wavefunctions,
calculates Lowdin charges, spilling parameter, projected DOS
(separated into up and down components for LSDA). Alternatively:
computes the local DOS(E) integrated in volumes given in input
(see tdosinboxes) or k-resolved DOS (see kresolveddos).
Atomic projections are written to file "atomic_proj.xml".
Structure of the input data:
============================
&PROJWFC
...
/
Namelist: &PROJWFC
prefix CHARACTER
Default: 'pwscf'
prefix of input file produced by pw.x (wavefunctions are needed)
[Back to Top]
outdir CHARACTER
Default: value of the ESPRESSO_TMPDIR environment variable if set; current directory ('./')
otherwise
directory containing the input data, i.e. the same as in pw.x
[Back to Top]
ngauss INTEGER
Default: 0
Type of gaussian broadening:
0 ... Simple Gaussian (default)
1 ... Methfessel-Paxton of order 1
https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html 1/5

3/27/26, 10:50 PM projwfc.x: input description
-1 ... "cold smearing" (Marzari-Vanderbilt-DeVita-Payne)
-99 ... Fermi-Dirac function
[Back to Top]
degauss REAL
Default: 0.0
gaussian broadening, Ry (not eV!)
[Back to Top]
Emin, Emax REAL
Default: (band extrema)
min & max energy (eV) for DOS plot
[Back to Top]
DeltaE REAL
energy grid step (eV)
[Back to Top]
lsym LOGICAL
Default: .true.
if .true. the projections are symmetrized,
the partial density of states are computed
if .false. the projections are not symmetrized, the partial
DOS can be computed only in the k-resolved case
[Back to Top]
diag_basis LOGICAL
Default: .false.
if .false. the projections of Kohn-Sham states are
done on the orthogonalized atomic orbitals
in the global XYZ coordinate frame.
if .true. the projections of Kohn-Sham states are
done on the orthogonalized atomic orbitals
that are rotated to the basis in which the
atomic occupation matrix is diagonal
(i.e. local XYZ coordinate frame).
[Back to Top]
pawproj LOGICAL
Default: .false.
if .true. use PAW projectors and all-electron PAW basis
functions to calculate weight factors for the partial
densities of states. Following Bloechl, PRB 50, 17953 (1994),
Eq. (4 & 6), the weight factors thus approximate the real
charge within the augmentation sphere of each atom.
Only for PAW, not implemented in the noncolinear case.
[Back to Top]
filpdos CHARACTER
Default: (value of prefix variable)
prefix for output files containing PDOS(E)
[Back to Top]
https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html 2/5

3/27/26, 10:50 PM projwfc.x: input description
filproj CHARACTER
Default: (standard output)
file containing the projections
[Back to Top]
lwrite_overlaps LOGICAL
Default: .false.
if .true., the overlap matrix of the atomic orbitals
prior to orthogonalization is written to "atomic_proj.xml".
Does not work together with parallel diagonalization:
for parallel runs, use "mpirun -np N projwfc.x -nd 1 ... "
[Back to Top]
lbinary_data LOGICAL
Default: .false.
CURRENTLY DISABLED.
if .true., write atomic projections to a binary file.
[Back to Top]
kresolveddos LOGICAL
Default: .false.
if .true. the k-resolved DOS is computed: not summed over
all k-points but written as a function of the k-point index.
In this case all k-point weights are set to unity
[Back to Top]
tdosinboxes LOGICAL
Default: .false.
if .true. compute the local DOS integrated in volumes
Volumes are defined as boxes with edges parallel to the unit cell,
containing the points of the (charge density) FFT grid included within
irmin and irmax, in the three dimensions:
from irmin(j,n) to irmax(j,n) for j=1,2,3 (n=1,n_proj_boxes).
[Back to Top]
n_proj_boxes INTEGER
Default: 1
number of boxes where the local DOS is computed
[Back to Top]
irmin(i,n_proj_boxes), (i,n_proj_boxes) = (1,1) . . . (3,n_proj_boxes) INTEGER
Default: 1 for each box
first point of the given box
BEWARE: irmin is a 2D array of the form: irmin(3,n_proj_boxes)
[Back to Top]
irmax(i,n_proj_boxes), (i,n_proj_boxes) = (1,1) . . . (3,n_proj_boxes) INTEGER
Default: 0 for each box
last point of the given box;
( 0 stands for the last point in the FFT grid )
https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html 3/5

3/27/26, 10:50 PM projwfc.x: input description
BEWARE: irmax is a 2D array of the form: irmax(3,n_proj_boxes)
[Back to Top]
plotboxes LOGICAL
Default: .false.
if .true., the boxes are written in output as xsf files with
3D datagrids, valued 1.0 inside the box volume and 0 outside
(visualize them as isosurfaces with isovalue 0.5)
[Back to Top]
Notes
Format of output files
Projections are written to standard output, and also to file
filproj if given as input.
The total DOS and the sum of projected DOS are written to file
"filpdos".pdos_tot.
* The format for the collinear, spin-unpolarized case and the
non-collinear, spin-orbit case is:
E DOS(E) PDOS(E)
...
* The format for the collinear, spin-polarized case is:
E DOSup(E) DOSdw(E) PDOSup(E) PDOSdw(E)
...
* The format for the non-collinear, non spin-orbit case is:
E DOS(E) PDOSup(E) PDOSdw(E)
...
In the collinear case and the non-collinear, non spin-orbit case
projected DOS are written to file "filpdos".pdos_atm#N(X)_wfc#M(l),
where N = atom number , X = atom symbol, M = wfc number, l=s,p,d,f
(one file per atomic wavefunction found in the pseudopotential file)
* The format for the collinear, spin-unpolarized case is:
E LDOS(E) PDOS_1(E) ... PDOS_2l+1(E)
...
where LDOS = \sum m=1,2l+1 PDOS_m(E)
and PDOS_m(E) = projected DOS on atomic wfc with component m
* The format for the collinear, spin-polarized case and the
non-collinear, non spin-orbit case is as above with
two components for both LDOS(E) and PDOS_m(E)
In the non-collinear, spin-orbit case (i.e. if there is at least one
fully relativistic pseudopotential) wavefunctions are projected
onto eigenstates of the total angular-momentum.
Projected DOS are written to file "filpdos".pdos_atm#N(X)_wfc#M(l_j),
where N = atom number , X = atom symbol, M = wfc number, l=s,p,d,f
and j is the value of the total angular momentum.
In this case the format is:
E LDOS(E) PDOS_1(E) ... PDOS_2j+1(E)
...
If kresolveddos=.true., the k-point index is prepended
to the formats above, e.g. (collinear, spin-unpolarized case)
ik E DOS(E) PDOS(E)
All DOS(E) are in states/eV plotted vs E in eV
https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html 4/5

3/27/26, 10:50 PM projwfc.x: input description
Orbital Order
Order of m-components for each l in the output:
1, cos(phi), sin(phi), cos(2*phi), sin(2*phi), .., cos(l*phi), sin(l*phi)
where phi is the azimuthal angle: x=r cos(theta)cos(phi), y=r cos(theta)sin(phi)
This is determined in file upflib/ylmr2.f90 that calculates spherical harmonics.
for l=1:
1 pz (m=0)
2 px (real combination of m=+/-1 with cosine)
3 py (real combination of m=+/-1 with sine)
for l=2:
1 dz2 (m=0)
2 dzx (real combination of m=+/-1 with cosine)
3 dzy (real combination of m=+/-1 with sine)
4 dx2-y2 (real combination of m=+/-2 with cosine)
5 dxy (real combination of m=+/-2 with sine)
Defining boxes for the Local DOS(E)
Boxes are specified using the variables irmin and irmax:
FFT grid points are included from irmin(j,n) to irmax(j,n)
for j=1,2,3 and n=1,...,n_proj_boxes
irmin and irmax range from 1 to nr1 or nr2 or nr3
Values larger than nr1/2/3 or smaller than 1 are folded
to the unit cell.
If irmax<irmin FFT grid points are included from 1 to irmax
and from irmin to nr1/2/3.
Important notices
The tetrahedron method is used if
- the input data file has been produced by pw.x using the option
occupations='tetrahedra', AND
- a value for degauss is not given as input to namelist &projwfc
* Gaussian broadening is used in all other cases:
- if degauss is set to some value in namelist &PROJWFC, that value
(and the optional value for ngauss) is used
- if degauss is NOT set to any value in namelist &PROJWFC, the
value of degauss and of ngauss are read from the input data
file (they will be the same used in the pw.x calculations)
- if degauss is NOT set to any value in namelist &PROJWFC, AND
there is no value of degauss and of ngauss in the input data
file, degauss=DeltaE (in Ry) and ngauss=0 will be used
Obsolete variables, ignored:
io_choice
smoothing
This file has been created by helpdoc utility on Wed Sep 03 14:29:00 CEST 2025.
https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html 5/5