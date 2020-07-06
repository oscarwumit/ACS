"""
This module analyzes cosmo-rs solvation free energy calculation results.
"""

import os
import argparse

from acs.common import read_yaml_file, write_yaml_file, read_cosmo_gsolv



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
    The analyze_cosmo_result executable function
    """
    # 0. Parse input
    # 0.1 Parse opt project info
    args = parse_command_line_arguments()
    # input_file_dir = os.path.abspath(os.path.dirname(args.file)) # todo: use this to deduce abs path if folder moved
    cosmo_project_input_file = args.file
    cosmo_project_info = read_yaml_file(cosmo_project_input_file)
    project_dir = cosmo_project_info['project_folder_path']
    cosmo_dir = os.path.join(project_dir, 'cosmo')

    # 1. Analyze cosmo result
    # 1.1. Parse Gsolv from cosmo tab file
    valid_conformer_hash_ids = cosmo_project_info['valid_conformer_hash_ids']
    for fingerprint in valid_conformer_hash_ids:
        tab_file_path = cosmo_project_info['conformers'][fingerprint]['file_path']['output']['solv_correction']
        g_solv = read_cosmo_gsolv(path=tab_file_path, use_hartree=True)
        cosmo_project_info['conformers'][fingerprint]['energy']['solv_correction'] = g_solv

        sp_after_opt = cosmo_project_info['conformers'][fingerprint]['energy'].get('sp_after_opt', None)
        # todo: deal with the case where sp_after_opt is None (maybe always parse sp result first)
        if sp_after_opt is not None:
            sp_include_solv_correction = sp_after_opt + g_solv
            cosmo_project_info['conformers'][fingerprint]['energy']['sp_include_solv_correction'] = \
                sp_include_solv_correction
        else:
            raise NotImplementedError

    # 2. Save fine opt project info (assume all calcs are done and all energies are parsed)
        full_project_info_path = os.path.join(project_dir, 'full_project_info.yml')
        write_yaml_file(path=full_project_info_path, content=cosmo_project_info)



if __name__ == '__main__':
    main()