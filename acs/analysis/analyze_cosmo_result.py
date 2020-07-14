"""
This module analyzes cosmo-rs solvation free energy calculation results.
"""

import os
import argparse

from acs.common import read_yaml_file, write_yaml_file, read_cosmo_gsolv, mkdir, gen_gaussian_irc_input_file
from copy import deepcopy
from acs.exceptions import ParserError
import shutil



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
    is_ts = cosmo_project_info['species']['is_ts']
    is_cal_sol = cosmo_project_info['calc_solvation_sp_correction']
    is_cal_fine = True if cosmo_project_info['level_of_theory']['fine_opt_freq'] else False

    # 1. Analyze cosmo result
    # 1.1. Parse Gsolv from cosmo tab file
    # todo: parse turbomole log file to check if single point calc terminate normally
    valid_conformer_hash_ids = cosmo_project_info['valid_conformer_hash_ids']
    final_valid_conformer_hash_ids = list(valid_conformer_hash_ids)
    for fingerprint in valid_conformer_hash_ids:
        tab_file_path = cosmo_project_info['conformers'][fingerprint]['file_path']['output']['solv_correction']
        try:
            g_solv = read_cosmo_gsolv(path=tab_file_path, use_hartree=True)
        except (FileNotFoundError, ParserError):
            final_valid_conformer_hash_ids.remove(fingerprint)
            continue
        cosmo_project_info['conformers'][fingerprint]['energy']['solv_correction'] = g_solv

        sp_after_opt = cosmo_project_info['conformers'][fingerprint]['energy'].get('sp_after_opt', None)
        # todo: deal with the case where sp_after_opt is None (maybe always parse sp result first)
        if sp_after_opt is not None:
            sp_include_solv_correction = sp_after_opt + g_solv
            cosmo_project_info['conformers'][fingerprint]['energy']['sp_include_solv_correction'] = \
                sp_include_solv_correction
        else:
            raise NotImplementedError
    cosmo_project_info['valid_conformer_hash_ids'] = tuple(final_valid_conformer_hash_ids)

    # 2. Final check
    # 2.1 Copy all Gaussian log files of valid conformers to a folder for user to conveniently check the TS manually
    # todo: flag low freq and move them to seperate folder
    check_dir = os.path.join(project_dir, 'check')
    mkdir(check_dir)
    if is_ts:
        check_neg_freq_dir = os.path.join(check_dir, 'check_neg_freq')
        mkdir(check_neg_freq_dir)
        neg_freq_warning_threshold = -800

    for fingerprint in final_valid_conformer_hash_ids:
        log_file_path = cosmo_project_info['conformers'][fingerprint]['file_path']['output']['fine_opt_freq']
        filename = os.path.split(log_file_path)[-1]
        check_file_path = os.path.join(check_dir, filename)
        shutil.copyfile(src=log_file_path, dst=check_file_path)

        if is_ts:
            neg_freq = cosmo_project_info['conformers'][fingerprint]['negative_frequencies'][0]
            if neg_freq > neg_freq_warning_threshold:
            # e.g., neg_freq = -300 might be problematic for H abstraction
                check_neg_freq_file_path = os.path.join(check_neg_freq_dir, filename)
                shutil.copyfile(src=log_file_path, dst=check_neg_freq_file_path)

    # 2.2 Prepare IRC input files for lowest energy conformer
    # 2.2.1 Find lowest energy conformer
    id_energy_list = list()
    for fingerprint in final_valid_conformer_hash_ids:
        if is_cal_sol:
            energy = cosmo_project_info['conformers'][fingerprint]['energy']['sp_include_solv_correction']
            id_energy_list.append((fingerprint, energy))
        else:
            energy = cosmo_project_info['conformers'][fingerprint]['energy']['sp_after_opt']
            id_energy_list.append((fingerprint, energy))

    low_id_energy = min(id_energy_list, key = lambda t: t[1])

    # 2.2.2 Write IRC input file
    irc_dir = os.path.join(project_dir, 'irc')
    mkdir(irc_dir)

    fingerprint = low_id_energy[0]

    forward_irc_file_name = '0_' + str(fingerprint) + '_forward_irc.gjf'
    reverse_irc_file_name = '1_' + str(fingerprint) + '_reverse_irc.gjf'

    forward_irc_file_path = os.path.join(irc_dir, forward_irc_file_name)
    reverse_irc_file_path = os.path.join(irc_dir, reverse_irc_file_name)

    charge = cosmo_project_info['species']['charge']
    multiplicity = cosmo_project_info['species']['multiplicity']
    if is_cal_fine:
        xyz_str = cosmo_project_info['conformers'][fingerprint]['xyz_str_after_fine_opt']
        level_of_theory = cosmo_project_info['level_of_theory']['fine_opt_freq']
    else:
        xyz_str = cosmo_project_info['conformers'][fingerprint]['xyz_str_after_opt']
        level_of_theory = cosmo_project_info['level_of_theory']['opt_freq']

    forward_irc_file_basename = forward_irc_file_name.split('.')[0]
    forward_irc_file = gen_gaussian_irc_input_file(name=forward_irc_file_basename,
                                                          xyz_str=xyz_str,
                                                          charge=charge,
                                                          multiplicity=multiplicity,
                                                          memory_mb=300000,
                                                          cpu_threads=40,
                                                          is_forward=True,
                                                          level_of_theory=level_of_theory,
                                                          comment=str(fingerprint),
                                                          )

    with open(forward_irc_file_path, 'w') as f:
        f.write(forward_irc_file)

    reverse_irc_file_basename = forward_irc_file_name.split('.')[0]
    reverse_irc_file = gen_gaussian_irc_input_file(name=reverse_irc_file_basename,
                                                          xyz_str=xyz_str,
                                                          charge=charge,
                                                          multiplicity=multiplicity,
                                                          memory_mb=300000,
                                                          cpu_threads=40,
                                                          is_forward=False,
                                                          level_of_theory=level_of_theory,
                                                          comment=str(fingerprint),
                                                          )

    with open(reverse_irc_file_path, 'w') as f:
        f.write(reverse_irc_file)

    # 3. Save final project info (assume all calcs are done and all energies are parsed)
    final_project_info = deepcopy(cosmo_project_info)
    all_ids = tuple(final_project_info['conformers'].keys())
    for fingerprint in all_ids:
        if fingerprint not in valid_conformer_hash_ids:
            del final_project_info['conformers'][fingerprint]

    del final_project_info['is_initial_sp_screening']
    del final_project_info['conformer_to_opt_hash_ids']
    del final_project_info['conformer_to_fine_opt_hash_ids']
    del final_project_info['colliding_conformer_after_opt_hash_ids']
    del final_project_info['crashing_conformer_in_opt_hash_ids']
    del final_project_info['colliding_conformer_after_fine_opt_hash_ids']
    del final_project_info['crashing_conformer_in_fine_opt_hash_ids']
    del final_project_info['conformer_to_calc_sp_after_opt_hash_ids']
    del final_project_info['conformer_to_calc_solv_hash_ids']

    full_project_info_path = os.path.join(project_dir, 'final_project_info.yml')
    write_yaml_file(path=full_project_info_path, content=final_project_info)


if __name__ == '__main__':
    main()