import os
import sys
from pathlib import Path
from ase import io
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.constraints import FixAtoms

metal = "VN"
run_root = Path(".")
run_root.mkdir(exist_ok=True)
pseudo_dir = Path("/scratch/anizami/QE_2/USPP/")
ecutwfc = 45.0
ecutrho = 8 * ecutwfc
kpts = (4, 4, 1)
pseudos = {
    "Ti": "ti_pbe_v1.4.uspp.F.UPF",
    "V": "v_pbe_v1.4.uspp.F.UPF",
    "Sc": "sc_pbe_v1.uspp.F.UPF",
    "Nb": "nb_pbe_v1.uspp.F.UPF",
    "Zr": "zr_pbe_v1.uspp.F.UPF",
    "N": "n_pbe_v1.2.uspp.F.UPF",
}

structure_path = Path("input_file/VN_relaxed.extxyz")
run_dir = run_root / "VN_relaxed"
run_dir.mkdir(parents=True, exist_ok=True)

atoms = io.read(structure_path)

# Fix atoms below Z = 10 Angstrom (bottom layer)
mask = atoms.positions[:, 2] < 10.0
fixed_atoms = FixAtoms(mask=mask)
atoms.set_constraint(fixed_atoms)

os.environ.setdefault("OMP_NUM_THREADS", "1")
profile = EspressoProfile(command="srun /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/MPI/gcc12/openmpi4/quantumespresso/7.5/bin/pw.x", pseudo_dir=str(pseudo_dir))
calc = Espresso(
    profile=profile,
    pseudopotentials=pseudos,
    input_data={
        "control": {
            "calculation": "relax",
            "prefix": "slab",
            "disk_io": "low",
            "tprnfor": True,
        },
        "system": {
            "ecutwfc": ecutwfc,
            "ecutrho": ecutrho,
            "input_dft": "rpbe",
            "vdw_corr": "dft-d3",
            "dftd3_version": 4,
            "occupations": "smearing",
            "smearing": "cold",
            "degauss": 0.005,
        },
        "electrons": {         # <-- ADD THIS BLOCK
            "conv_thr": 5.0e-6,
        },
    },
    kpts=kpts,
    pseudo_dir=str(pseudo_dir),
    directory=str(run_dir),
)
atoms.calc = calc
energy = atoms.get_potential_energy()
output_path = run_dir / "VN_relaxed.extxyz"
io.write(output_path, atoms, format="extxyz")
print("relaxation complete",
    f"structure -> {output_path}",
    f"energy -> {energy:.6f} eV",
)