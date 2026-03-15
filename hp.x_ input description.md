

## Input File Description
Program: hp.x / PWscf / Quantum ESPRESSO (version: 7.5)
## TABLE OF CONTENTS
## INTRODUCTION
## &INPUTHP
prefix | outdir | iverbosity | max_seconds | nq1 | nq2 | nq3 | skip_equivalence_q |
determine_num_pert_only | determine_q_mesh_only | find_atpert | docc_thr | skip_type |
equiv_type | perturb_only_atom | start_q | last_q | sum_pertq | compute_hp | conv_thr_chi |
thresh_init | ethr_nscf | niter_max | alpha_mix(i) | nmix | num_neigh | lmin | rmax | dist_thr |
no_metq0
## INTRODUCTION
Input data format: { } = optional, [ ] = it depends, # = comment
Structure of the input data:
## ===============================================================================
## &INPUTHP
## ...
## /

Namelist: &INPUTHP
prefixCHARACTER
## Default:'pwscf'
Prepended to input/output filenames; must be the same
used in the calculation of unperturbed system.

[Back to Top]
outdirCHARACTER
Default:value of the ESPRESSO_TMPDIR environment variable if set;
current directory ('./') otherwise
Directory containing input, output, and scratch files;
must be the same as specified in the calculation of
the unperturbed system.

[Back to Top]
iverbosityINTEGER
## Default:1
= 1 : minimal output
= 2 : as above + symmetry matrices, final response
matrices chi0 and chi1 and their inverse matrices,
full U matrix
= 3 : as above + various detailed info about the NSCF
calculation at k and k+q
= 4 : as above + response occupation matrices at every
iteration and for every q point in the star

[Back to Top]
max_secondsREAL
## Default:1.d7
2/15/26, 12:16 PMhp.x: input description
https://www.quantum-espresso.org/Doc/INPUT_HP.html1/6

Maximum allowed run time before the job stops smoothly.

[Back to Top]
nq1, nq2, nq3INTEGER
## Default:1,1,1
Parameters of the Monkhorst-Pack grid (no offset).
Same meaning as for nk1, nk2, nk3 in the input of pw.x.

[Back to Top]
skip_equivalence_qLOGICAL
## Default:.false.
If .true. then the HP code will skip the equivalence
analysis of q points, and thus the full grid of q points
will be used. Otherwise the symmetry is used to determine
equivalent q points (star of q), and then perform
calculations only for inequivalent q points.

[Back to Top]
determine_num_pert_onlyLOGICAL
## Default:.false.
## See:find_atpert
If .true. determines the number of perturbations
(i.e. which atoms will be perturbed) and exits smoothly
without performing any calculation. For DFT+U+V, it also
determines the indices of inter-site couples.

[Back to Top]
determine_q_mesh_onlyLOGICAL
## Default:.false.
## See:perturb_only_atom
If .true. determines the number of q points
for a given perturbed atom and exits smoothly.
This keyword can be used only if perturb_only_atom
is set to .true.

[Back to Top]
find_atpertINTEGER
## Default:1
Method for searching of atoms which must be perturbed.
1 = Find how many inequivalent Hubbard atoms there are
by analyzing unperturbed occupations.
2 = Find how many Hubbard atoms to perturb based on
how many different Hubbard atomic types there are.
Warning: atoms which have the same type but which
are inequivalent by symmetry or which have different
occupations will not be distinguished in this case
(use option 1 or 3 instead).
3 = Find how many inequivalent Hubbard atoms
there are using symmetry. Atoms which have the
same type but are not equivalent by symmetry will
be distinguished in this case.
4 = Perturb all Hubbard atoms (the most expensive option)

[Back to Top]
docc_thrREAL
Default:5.D-5
2/15/26, 12:16 PMhp.x: input description
https://www.quantum-espresso.org/Doc/INPUT_HP.html2/6

Threshold for a comparison of unperturbed occupations
which is needed for the selection of atoms which must
be perturbed. Can be used only when find_atpert = 1.

[Back to Top]
skip_type(i), i=1,ntypLOGICAL
## Default:skip_type(i) = .false.
## See:equiv_type
skip_type(i), where i runs over types of atoms.
If skip_type(i)=.true. then no linear-response
calculation will be performed for the i-th atomic type:
in this case equiv_type(i) must be specified, otherwise
the HP code will stop. This option is useful if the
system has atoms of the same type but opposite spin
pollarizations (anti-ferromagnetic case).
This keyword cannot be used when find_atpert = 1.

[Back to Top]
equiv_type(i), i=1,ntypINTEGER
## Default:equiv_type(i) = 0
## See:skip_type
equiv_type(i), where i runs over types of atoms.
equiv_type(i)=j, will make type i equivalent to type j
(useful when nspin=2). Such a merging of types is done
only at the post-processing stage.
This keyword cannot be used when find_atpert = 1.

[Back to Top]
perturb_only_atom(i), i=1,ntypLOGICAL
## Default:perturb_only_atom(i) = .false.
## See:compute_hp
If perturb_only_atom(i)=.true. then only the i-th
atom will be perturbed and considered in the run.
This variable is useful when one wants to split
the whole calculation on parts.
Note: this variable has a higher priority than skip_type.

[Back to Top]
start_qINTEGER
## Default:1
See:last_q, sum_pertq
Computes only the q points from start_q to last_q.
IMPORTANT: start_q must be smaller or equal to
the total number of q points found.

[Back to Top]
last_qINTEGER
Default:number of q points
See:start_q, sum_pertq
Computes only the q points from start_q to last_q.
IMPORTANT: last_q must be smaller or equal to
the total number of q points found.

2/15/26, 12:16 PMhp.x: input description
https://www.quantum-espresso.org/Doc/INPUT_HP.html3/6

[Back to Top]
sum_pertqLOGICAL
## Default:.false.
See:start_q, last_q, perturb_only_atom
If it is set to .true. then the HP code will collect
pieces of the response occupation matrices for all
q points. This variable should be used only when
start_q, last_q and perturb_only_atom are used.

[Back to Top]
compute_hpLOGICAL
## Default:.false.
## See:perturb_only_atom
If it is set to .true. then the HP code will collect
pieces of the chi0 and chi matrices (which must have
been produced in previous runs) and then compute
Hubbard parameters. The HP code will look for files
tmp_dir/HP/prefix.chi.i.dat. Note that all files
prefix.chi.i.dat (where i runs over all perturbed
atoms) must be placed in one folder tmp_dir/HP/.
compute_hp=.true. must be used only when the
calculation was parallelized over perturbations.

[Back to Top]
conv_thr_chiREAL
Default:1.D-5
Convergence threshold for the response function chi,
which is defined as a trace of the response
occupation matrix.

[Back to Top]
thresh_initREAL
Default:1.D-14
Initial threshold for the solution of the linear
system (first iteration). Needed to converge the
bare (non-interacting) response function chi0.
The specified value will be multiplied by the
number of electrons in the system.

[Back to Top]
ethr_nscfREAL
Default:1.D-11
Threshold for the convergence of eigenvalues during
the iterative diagonalization of the Hamiltonian in
the non-self-consistent-field (NSCF) calculation at
k and k+q points. Note, this quantity is NOT extensive.

[Back to Top]
niter_maxINTEGER
## Default:100
Maximum number of iterations in the iterative
solution of the linear-response Kohn-Sham equations.

[Back to Top]
alpha_mix(i)REAL
2/15/26, 12:16 PMhp.x: input description
https://www.quantum-espresso.org/Doc/INPUT_HP.html4/6

## Default:alpha_mix(1)=0.3
Mixing parameter (for the i-th iteration) for updating
the response SCF potential using the modified Broyden
method. See: D.D. Johnson, PRB 38, 12807 (1988).

[Back to Top]
nmixINTEGER
## Default:4
Number of iterations used in potential mixing
using the modified Broyden method. See:
D.D. Johnson, PRB 38, 12807 (1988).

[Back to Top]
num_neighINTEGER
## Default:6
Number of nearest neighbors of every Hubbard atom which
will be considered when writting Hubbard V parameters to
the file parameters.out, which can be used in the
subsequent DFT+U+V calculation. This keyword is used only
for DFT+U+V (post-processing stage).

[Back to Top]
lminINTEGER
## Default:2
Minimum value of the orbital quantum number of the Hubbard
atoms starting from which (and up to the maximum l in the
system) Hubbard V will be written to the file parameters.out.
lmin refers to the orbital quantum number of the atom
corresponding to the first site-index in Hubbard_V(:,:,:).
This keyword is used only for DFT+U+V and only
in the post-processing stage. Example: lmin=1 corresponds to
writing to file V between e.g. oxygen (with p states) and its
neighbors, and including V between transition metals (with d
states) and their neighbors. Instead, when lmin=2 only the
latter will be written to parameters.out.

[Back to Top]
rmaxREAL
Default:100.D0
Maximum distance (in Bohr) between two atoms to search
neighbors (used only at the postprocessing step for
DFT+U+V). This keyword is useful when there
are e.g. defects in the system.

[Back to Top]
dist_thrREAL
Default:6.D-4
Threshold (in Bohr) for comparing inter-atomic distances
when reconstructing the missing elements of the response
susceptibility in the post-processing step.

[Back to Top]
no_metq0LOGICAL
## Default:.false.
If .true. the metallic response term at q=0 is ignored
(i.e. the last term in Eq. (22) in PRB 103, 045141 (2021)).
This is useful for magnetic insulators to avoid the divergence
2/15/26, 12:16 PMhp.x: input description
https://www.quantum-espresso.org/Doc/INPUT_HP.html5/6

of the calculation.

[Back to Top]
This file has been created by helpdoc utility on Wed Sep 03 14:25:41 CEST 2025.
2/15/26, 12:16 PMhp.x: input description
https://www.quantum-espresso.org/Doc/INPUT_HP.html6/6