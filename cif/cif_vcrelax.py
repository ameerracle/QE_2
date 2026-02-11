import os
from pathlib import Path
from ase import io
from ase.build import add_vacuum
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.optimize import QuasiNewton

metals_compounds = ['TiN','VN', 'ScN','NbN','ZrN']
metals = ['Ti','V','Sc','Nb','Zr']
run_root = Path(".")
run_root.mkdir(exist_ok=True)
pseudo_dir = Path("/scratch/anizami/QE_2/USPP/")
ecutwfc = 45.0
ecutrho = 8 * ecutwfc
kpts = (4, 4, 4)
pseudos = {
    "Ti": "ti_pbe_v1.4.uspp.F.UPF",
    "V": "v_pbe_v1.4.uspp.F.UPF",
    "Sc": "sc_pbe_v1.uspp.F.UPF",
    "Nb": "nb_pbe_v1.uspp.F.UPF",
    "Zr": "zr_pbe_v1.uspp.F.UPF",
    "N": "n_pbe_v1.2.uspp.F.UPF",
}

for metal in metals_compounds:
    structure_path = Path(f"input_file/{metal}_relaxed.extxyz")
    run_dir = run_root / (metal +'_vcrelax')
    run_dir.mkdir(parents=True, exist_ok=True)

    atoms = io.read(structure_path)

    os.environ.setdefault("OMP_NUM_THREADS", "1")
    #profile = EspressoProfile(command="srun /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/MPI/gcc12/openmpi4/quantumespresso/7.5/bin/pw.x", pseudo_dir=str(pseudo_dir))

    #profile = EspressoProfile(command="srun /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/MPI/gcc12/openmpi4/quantumespresso/7.5/bin/pw.x", pseudo_dir=str(pseudo_dir))
    calc = Espresso(
        profile=profile,
        pseudopotentials=pseudos,
        input_data={
            "control": {"calculation": "vc-relax", "prefix": "cif", 'disk_io':'low'},
            "system": {"ecutwfc": ecutwfc, "ecutrho": ecutrho},
            "input_dft": 'pbe',
            "vdw_corr": 'dft-d3',
            "dftd3_version":4,
            "occupations": "smearing",
            "smearing": "cold",
            "degauss": 0.005, # in Ry,
        },

        kpts=kpts,
        pseudo_dir=str(pseudo_dir),
        directory=str(run_dir),
    )
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    extxyz_path = run_dir / f"{metal}_relaxed.extxyz"
    io.write(extxyz_path, atoms, format="extxyz")
    print("vc-relax complete",
        f"structure -> {extxyz_path}",
        f"energy -> {energy:.6f} eV",
    )