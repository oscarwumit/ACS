"""
This module analyzes geometry optimization and frequency calculation results.
"""

import os
import yaml
import argparse
from itertools import product
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Rectangle

from acs.common import read_yaml_file, write_yaml_file, mkdir, update_selected_keys, gen_molpro_sp_input_file, \
    process_opt_freq_output, cluster_confs_by_rmsd_with_id, gen_gaussian_cosmo_sp_input_file, \
    gen_cosmo_rs_input_file, gen_orca_dlpno_sp_input_file
from acs.script import default_job_info_dict_after_initial_sp_screening, molpro_slurm_array_script, \
    default_conformer_info_dict_after_initial_sp_screening, cosmo_slurm_array_script, \
    turbomole_cosmo_slurm_array_script, orca_slurm_array_script
from acs.exceptions import ParserError
from acs.converter.geom import xyz_dict_to_xyz_str
from copy import deepcopy


ESS_SOFTWARE = 'qchem'

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
    The analyze_opt_freq_result executable function
    """
    # todo: TS volume check
    # todo: check isomporphism, collision, bond breaking using distance matrix; for connected atoms,
    #  the bond distance should not become too large; for not-connected atoms, the bond disntance should not become too
    #  close; for TS bond, the threhold might be a bit larger; since we preserve atom order, these checks can be done via
    #  a pre-defined mask with entries deduced from the starting conformer (e.g., get nominal bond distance, then scale
    #  it with a percent threhold for change) and applies the mask on the distance matrix of new conformers
    # 0. Parse input
    # 0.1 Parse opt project info
    args = parse_command_line_arguments()
    opt_project_input_file = args.file
    opt_project_info = read_yaml_file(opt_project_input_file)
    project_dir = opt_project_info['project_folder_path']

    # 1. Analyze opt freq result
    # Pre-allocation and set ups
    is_ts = opt_project_info['species']['is_ts']
    charge = opt_project_info['species']['charge']
    multiplicity = opt_project_info['species']['multiplicity']

    # 1.1. Parse info from opt freq output. Assumes no fine opt freq job was run
    opted_conf_fingerprints = opt_project_info['conformer_to_opt_hash_ids']
    fingerprint_to_all_opt_log_info_dict = dict()
    crashing_conformer_hash_ids = list()
    normal_termination_conformer_hash_ids = list()

    # 1.1.1 load energy, freq, geom info from log
    for fingerprint in opted_conf_fingerprints:
        opt_input_file_path = opt_project_info['conformers'][fingerprint]['file_path']['input']['opt_freq']
        dir_name, file_name = os.path.split(opt_input_file_path)
        file_basename, file_extension = os.path.splitext(file_name)
        new_file_name = file_basename + '.log'
        opt_output_file_path = os.path.join(dir_name, new_file_name)
        if not os.path.exists(opt_output_file_path):
            raise
        opt_project_info['conformers'][fingerprint]['file_path']['output']['opt_freq'] = opt_output_file_path
        try:
            fingerprint_to_all_opt_log_info_dict[fingerprint] = \
                process_opt_freq_output(logfile=opt_output_file_path, ess_software=ESS_SOFTWARE, is_ts=is_ts)
            normal_termination_conformer_hash_ids.append(fingerprint)
        except ParserError as e:
            print(f'Parser Error! {e}')
            crashing_conformer_hash_ids.append(fingerprint)
        except TypeError:
            crashing_conformer_hash_ids.append(fingerprint)

    # 1.1.2 check if crash
    for fingerprint in opted_conf_fingerprints:
        opt_project_info['conformers'][fingerprint]['is_crashing'] = \
            False if fingerprint in normal_termination_conformer_hash_ids else True
    opt_project_info['crashing_conformer_in_opt_hash_ids'] = tuple(crashing_conformer_hash_ids)


    # 1.1.3 check if geometry makes sense
    # 1.1.3.1 todo: for non-TS, check isomorphism
    if not is_ts:
        pass
    # 1.1.3.2 todo: for TS, check volume, bond breaking, freq check based on reaction family
    else:
        pass

    # 1.1.4 cluster conformers
    labeled_xyzs = list()
    for fingerprint in normal_termination_conformer_hash_ids:
        labeled_xyzs.append((fingerprint, fingerprint_to_all_opt_log_info_dict[fingerprint]['xyz_dict']))
    distinct_id_xyz_pair = cluster_confs_by_rmsd_with_id(labeled_xyzs=tuple(labeled_xyzs),
                                                         rmsd_threshold=1e-2,
                                                         )
    distinct_conformer_hash_ids = [i[0] for i in distinct_id_xyz_pair]
    for fingerprint in normal_termination_conformer_hash_ids:
        opt_project_info['conformers'][fingerprint]['is_distinct'] = \
            True if fingerprint in distinct_conformer_hash_ids else False

    # 1.1.5 populate valid conformer information
    # checked conf should also pass other tests such as isomorphism, bond distance, volume etc. Here only use rmsd.
    # todo: save frequency for jobs that do not need fine opt
    checked_conformer_hash_ids = tuple(distinct_conformer_hash_ids)

    for fingerprint in checked_conformer_hash_ids:
        xyz_dict = fingerprint_to_all_opt_log_info_dict[fingerprint]['xyz_dict']
        opt_project_info['conformers'][fingerprint]['arc_xyz_after_opt'] = xyz_dict
        opt_project_info['conformers'][fingerprint]['xyz_str_after_opt'] = xyz_dict_to_xyz_str(xyz_dict)

        opt_project_info['conformers'][fingerprint]['energy']['end_of_opt'] = \
            fingerprint_to_all_opt_log_info_dict[fingerprint]['electronic_energy']['hartree']

        # store the frequencies from the opt job
        # currently, running a fine_opt job will override these values
        opt_project_info['conformers'][fingerprint]['frequencies'] = \
            fingerprint_to_all_opt_log_info_dict[fingerprint]['freq'][0]
        opt_project_info['conformers'][fingerprint]['negative_frequencies'] = \
            fingerprint_to_all_opt_log_info_dict[fingerprint]['freq'][1]

    # 2. Generate input files for high level sp
    sp_dir = os.path.join(project_dir, 'sp')
    mkdir(sp_dir)
    sp_level_of_theroy = opt_project_info['level_of_theory']['sp_after_opt']

    conformer_to_sp_hash_ids = deepcopy(checked_conformer_hash_ids)
    opt_project_info['conformer_to_sp_hash_ids'] = conformer_to_sp_hash_ids

    for i, fingerprint in enumerate(opt_project_info['conformer_to_sp_hash_ids']):
        sp_input_file_name = str(i) + '_' + str(fingerprint) + '_sp.in'
        sp_input_file_path = os.path.join(sp_dir, sp_input_file_name)

        xyz_str = opt_project_info['conformers'][fingerprint]['xyz_str_after_opt']

        sp_input_file_basename = sp_input_file_name.split('.')[0]
        sp_input_file = gen_molpro_sp_input_file(xyz_str=xyz_str,
                                                 method=sp_level_of_theroy.split('/')[0],
                                                 basis=sp_level_of_theroy.split('/')[1],
                                                 spin=multiplicity-1,
                                                 charge=charge,
                                                 )

        opt_project_info['conformers'][fingerprint]['file_path']['input']['sp_after_opt'] = sp_input_file_path

        with open(sp_input_file_path, 'w') as f:
            f.write(sp_input_file)

    # 2.1 Array job submission script
    name = opt_project_info['species']['name']
    last_job_num = len(opt_project_info['conformer_to_sp_hash_ids']) - 1
    sub_script = deepcopy(molpro_slurm_array_script)
    sub_script = sub_script.format(job_name=name, last_job_num=str(last_job_num))
    sub_script_file_path = os.path.join(sp_dir, 'submit_molpro_array.sh')
    with open(sub_script_file_path, 'w') as f:
        f.write(sub_script)

    # 2.2 save sp project info
    sp_project_info_path = os.path.join(sp_dir, 'sp_project_info.yml')
    write_yaml_file(path=sp_project_info_path, content=opt_project_info)

if __name__ == '__main__':
    main()