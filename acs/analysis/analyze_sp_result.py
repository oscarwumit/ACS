"""
This module analyzes sp after opt results.
"""

import os
import argparse

from acs.common import read_yaml_file, write_yaml_file, get_e_elect_from_log
from copy import deepcopy
from acs.exceptions import ParserError


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
    sp_project_input_file = args.file
    sp_project_info = read_yaml_file(sp_project_input_file)
    project_dir = sp_project_info['project_folder_path']
    sp_dir = os.path.join(project_dir, 'sp')

    # 1. Analyze sp result
    # 1.1. Parse sp from log file
    # todo: parse turbomole log file to check if single point calc terminate normally
    valid_conformer_hash_ids = sp_project_info['valid_conformer_hash_ids']
    final_valid_conformer_hash_ids = list(valid_conformer_hash_ids)
    for fingerprint in valid_conformer_hash_ids:
        sp_file_path = sp_project_info['conformers'][fingerprint]['file_path']['output']['sp_after_opt']
        try:
            e_elect_dict = get_e_elect_from_log(sp_file_path)
        except (FileNotFoundError, ParserError):
            final_valid_conformer_hash_ids.remove(fingerprint)
            continue
        sp_project_info['conformers'][fingerprint]['energy']['sp_after_opt'] = e_elect_dict['hartree']

    sp_project_info['valid_conformer_hash_ids'] = tuple(final_valid_conformer_hash_ids)

    # 1.2. Save sp project info (assume all calcs are done and all energies are parsed)
    sp_result_info = deepcopy(sp_project_info)
    sp_result_info_path = os.path.join(project_dir, 'sp_result_info.yml')
    write_yaml_file(path=sp_result_info_path, content=sp_result_info)


if __name__ == '__main__':
    main()