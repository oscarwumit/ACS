"""
This module runs initial sp screening for conformers.
"""

import yaml
import os
import argparse
import numpy as np

from ase.calculators.psi4 import Psi4

from joblib import Parallel, delayed
from acs.converter.geom import xyz_dict_to_ase_atom
from acs.common import read_yaml_file, write_yaml_file


def geom_and_calc_producer(bookkeep):
    # todo: different default values based on ts or non-ts
    charge = bookkeep['species'].get('charge', 0)
    multiplicity = bookkeep['species'].get('multiplicity', 2)
    # todo: deal with other multiplicity
    reference = 'uhf' if multiplicity == 2 else 'rhf'  # assume only 1 or 2 here, need consider other cases
    initial_screening_sp_level_of_theory = bookkeep['level_of_theory'].get('initial_screening_sp', 'wb97x-d/def2-svp')
    method, basis = initial_screening_sp_level_of_theory.split('/')
    memory = bookkeep.get('memory', '7000MB')
    num_thread = 1

    calc = set_up_psi4_calculator(charge=charge,
                                    multiplicity=multiplicity,
                                    reference=reference,
                                    method=method,
                                    basis=basis,
                                    memory=memory,
                                    num_thread=num_thread,
                                  )

    for fingerprint in bookkeep['conformer_to_screen_hash_ids']:
        ase_atom = xyz_dict_to_ase_atom(bookkeep['conformers'][fingerprint]['arc_xyz'])
        yield fingerprint, ase_atom, calc


def set_up_psi4_calculator(charge,
                           multiplicity,
                           reference,
                           method,
                           basis,
                           memory,
                           num_thread,
                           ):

    calc = Psi4(charge=charge,
                multiplicity=multiplicity,
                reference=reference,
                method=method,
                basis=basis,
                memory=memory,
                num_thread=num_thread,
                )
    return calc


def get_ase_energy(fingerprint, ase_atom, calc):
    try:
        ase_atom.calc = calc
        energy = ase_atom.get_potential_energy()
        return fingerprint, energy
    except:
        return fingerprint, np.nan


def parse_command_line_arguments(command_line_args=None):
    """
    Parse command-line arguments.

    Args:
        command_line_args: The command line arguments.

    Returns:
        The parsed command-line arguments by key words.
    """

    parser = argparse.ArgumentParser(description='Automatic Conformer Search (ACS)')
    parser.add_argument('file', metavar='FILE', type=str, nargs=1, help='a file describing the job to execute')
    args = parser.parse_args(command_line_args)
    args.file = args.file[0]

    return args


def main():
    """
    The run_screening executable function
    """
    # todo: write log
    # todo: exceptions

    # 0. Parse input
    args = parse_command_line_arguments()
    screening_input_file = args.file
    input_file_dir = os.path.abspath(os.path.dirname(args.file))
    project_info = read_yaml_file(screening_input_file)

    # 1. Run job in parallel
    result = Parallel(n_jobs=-1)(delayed(get_ase_energy)(*data) for data in geom_and_calc_producer(project_info))
    # 2. Save results
    crashing_conformer_hash_ids = list()
    for fingerprint, energy in result:
        project_info['conformers'][fingerprint]['is_crashing'] = np.isnan(energy)
        if np.isnan(energy):
            crashing_conformer_hash_ids.append(fingerprint)
        project_info['conformers'][fingerprint]['initial_screening_sp_energy'] = energy
    project_info['crashing_conformer_hash_ids'] = tuple(crashing_conformer_hash_ids)
    outfile_path = os.path.join(input_file_dir, 'initial_conf_screening_result.yml')
    write_yaml_file(path=outfile_path, content=project_info)


if __name__ == '__main__':
    main()





