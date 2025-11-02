import os
from ase.build import bulk
from ase import Atoms
import math

# -------------------- Setup directories --------------------
def get_output_dir(subfolder='xyz_outputs'):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, subfolder)
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

# -------------------- Supercell creation --------------------
def build_supercell(base_cell: Atoms, size, use_repeat=True):
    """Create a supercell from base_cell.
    
    Args:
        base_cell: ASE Atoms object
        size: tuple for repeat, e.g., (2,2,2)
        use_repeat: if True, use .repeat(), else use * operator
    """
    if use_repeat:
        return base_cell.repeat(size)
    else:
        return base_cell * size

# -------------------- NV center introduction --------------------
def add_nv_center(cell: Atoms, substitution_index=0, vacancy_index=1):
    """Create NV center by substituting a carbon with nitrogen and removing a neighbor."""
    nv_cell = cell.copy()
    nv_cell[substitution_index].symbol = 'N'
    nv_cell.pop(vacancy_index)
    return nv_cell

# -------------------- File writing --------------------
def write_xyz(cell: Atoms, filename: str):
    cell.write(filename)
    print(f"{os.path.basename(filename)} has {len(cell)} atoms")


def write_qe(cell: Atoms, filename: str):
    """
    Write a Quantum ESPRESSO input file for SCF calculation of the given cell.
    
    INPUT DATA REFERENCED FROM YU JIN'S 216 ATOM SUPERCELL CALCULATIONS
    - https://notebook.rcc.uchicago.edu/files/arXiv.2309.03513/Datasets/PointDefects/NV-_diamond/PBE_216_Relaxation/1E/
    """
    pseudo_dict = {'C': 'C_ONCV_PBE-1.0.upf', 'N': 'N_ONCV_PBE-1.0.upf'}

    number_of_atoms = len(cell)
    # number bands >= (N_valence_electrons/2) + extra bands for unoccupied states
    # For C (4 valence e-) and N (5 valence e-), total valence e- = 4*(nat-1) + 5*1 = 4*nat + 1
    number_of_bands = math.ceil((4 * (number_of_atoms - 1) + 5) / 2) + 25  # for NV center
    print("Number of bands set to:", number_of_bands)
    
    input_data = {
        'control': {
            'calculation': 'scf',
            'verbosity': 'high',
            'tprnfor': True,
            'tstress': True,
            'wf_collect': True,
            'prefix': 'pwscf',
            'outdir': './',
            'pseudo_dir': './',
            'forc_conv_thr': 0.00039,
            'nstep': 100
        },
        'system': {
            'ibrav': 0,
            'ecutwfc': 85,
            'tot_charge': -1,
            'nspin': 2,
            'tot_magnetization': 2.0,
            'nbnd': number_of_bands,
            'nat': len(cell),
            'ntyp': len(set(cell.get_chemical_symbols()))
        },
        'electrons': {
            'conv_thr': 1e-8,
            'diago_full_acc': True
        },
        'ions': {
            'ion_dynamics': 'bfgs',
            'ion_positions': 'default'
        }
    }

    # ASE can write espresso input
    cell.write(
        filename,
        format='espresso-in',
        pseudopotentials=pseudo_dict,
        input_data=input_data,
        kpts="gamma",  # Gamma point
        cell_dependent=True  # ensures CELL_PARAMETERS block is written
    )

    print(f"{os.path.basename(filename)} has {len(cell)} atoms")

# -------------------- Main script --------------------
def main():
    xyz_output_dir = get_output_dir('xyz_outputs')
    qe_output_dir = get_output_dir('qe_inputs')

    # Build conventional diamond cell
    print("- CREATING CONVENTIONAL DIAMOND CELL -")
    diamond = bulk('C', 'diamond', a=3.567, cubic=True)
    print(diamond)
    print("Number of atoms in the conventional cell:", len(diamond))

    # Create and write 2x2x2 supercell
    supercell_2x2x2 = build_supercell(diamond, (2, 2, 2), use_repeat=True)
    write_xyz(supercell_2x2x2, os.path.join(xyz_output_dir, 'diamond_2x2x2.xyz'))
    
    # Create and write 3x3x3 supercell
    supercell_3x3x3 = build_supercell(diamond, (3, 3, 3), use_repeat=False)
    write_xyz(supercell_3x3x3, os.path.join(xyz_output_dir, 'diamond_3x3x3.xyz'))
    
    print("\n- INTRODUCING NV CENTERS -")
    # NV center in 2x2x2
    nv_2x2x2 = add_nv_center(supercell_2x2x2)
    write_xyz(nv_2x2x2, os.path.join(xyz_output_dir, 'nv_diamond_2x2x2.xyz'))
    write_qe(nv_2x2x2, os.path.join(qe_output_dir, 'pw_nv_diamond_2x2x2.in'))

    # NV center in 3x3x3
    nv_3x3x3 = add_nv_center(supercell_3x3x3)
    write_xyz(nv_3x3x3, os.path.join(xyz_output_dir, 'nv_diamond_3x3x3.xyz'))
    write_qe(nv_3x3x3, os.path.join(qe_output_dir, 'pw_nv_diamond_3x3x3.in'))

if __name__ == "__main__":
    main()
