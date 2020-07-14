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

from acs.common import read_yaml_file, write_yaml_file, mkdir, update_selected_keys, gen_gaussian_optfreq_input_file, \
    process_gaussian_opt_freq_output, cluster_confs_by_rmsd_with_id, gen_gaussian_cosmo_sp_input_file, \
    gen_cosmo_rs_input_file
from acs.script import default_job_info_dict_after_initial_sp_screening, g16_slurm_array_script, \
    default_conformer_info_dict_after_initial_sp_screening, cosmo_slurm_array_script, turbomole_cosmo_slurm_array_script
from acs.exceptions import ParserError
from acs.converter.geom import xyz_dict_to_xyz_str
from copy import deepcopy



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
    # input_file_dir = os.path.abspath(os.path.dirname(args.file)) # todo: use this to deduce abs path if folder moved
    opt_project_input_file = args.file
    opt_project_info = read_yaml_file(opt_project_input_file)
    project_dir = opt_project_info['project_folder_path']
    # opt_dir = os.path.join(project_dir, 'opt')
    opt_fine_dir = os.path.join(project_dir, 'opt_fine')
    # 0.2 Determine fine opt status
    fine_opt_level_of_theory = opt_project_info['level_of_theory']['fine_opt_freq']
    if fine_opt_level_of_theory is None:
        has_fine_opted = None  # means fine opt is not specified/required
    else:
        if os.path.exists(opt_fine_dir):
            has_fine_opted = True  # assumes fine opt job has done, can use better method here e.g., check log
        else:
            has_fine_opted = False

    # 1. Analyze opt freq result
    # Pre-allocation and set ups

    is_ts = opt_project_info['species']['is_ts']
    charge = opt_project_info['species']['charge']
    multiplicity = opt_project_info['species']['multiplicity']

    # 1.1. Parse info from opt freq output
    # (a) Assume we are analyzing regular opt results
    if not has_fine_opted:
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
                    process_gaussian_opt_freq_output(logfile=opt_output_file_path, is_ts=is_ts)
                normal_termination_conformer_hash_ids.append(fingerprint)
            except ParserError:
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

    # (b) Assume we are analyzing fine opt results
    # todo: cosolidate code with regular opt
    else:
        opted_conf_fingerprints = opt_project_info['conformer_to_fine_opt_hash_ids']
        fingerprint_to_all_opt_log_info_dict = dict()
        crashing_conformer_hash_ids = list()
        normal_termination_conformer_hash_ids = list()

        # 1.1.1 load energy, freq, geom info from log
        for fingerprint in opted_conf_fingerprints:
            opt_input_file_path = opt_project_info['conformers'][fingerprint]['file_path']['input']['fine_opt_freq']
            dir_name, file_name = os.path.split(opt_input_file_path)
            file_basename, file_extension = os.path.splitext(file_name)
            new_file_name = file_basename + '.log'
            opt_output_file_path = os.path.join(dir_name, new_file_name)
            if not os.path.exists(opt_output_file_path):
                raise
            opt_project_info['conformers'][fingerprint]['file_path']['output']['fine_opt_freq'] = opt_output_file_path
            try:
                fingerprint_to_all_opt_log_info_dict[fingerprint] = \
                    process_gaussian_opt_freq_output(logfile=opt_output_file_path, is_ts=is_ts)
                normal_termination_conformer_hash_ids.append(fingerprint)
            except ParserError:
                crashing_conformer_hash_ids.append(fingerprint)

        # 1.1.2 check if crash
        for fingerprint in opted_conf_fingerprints:
            opt_project_info['conformers'][fingerprint]['is_crashing'] = \
                False if fingerprint in normal_termination_conformer_hash_ids else True
        opt_project_info['crashing_conformer_in_fine_opt_hash_ids'] = tuple(crashing_conformer_hash_ids)

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
        checked_conformer_hash_ids = tuple(distinct_conformer_hash_ids)

        for fingerprint in checked_conformer_hash_ids:
            xyz_dict = fingerprint_to_all_opt_log_info_dict[fingerprint]['xyz_dict']
            opt_project_info['conformers'][fingerprint]['arc_xyz_after_fine_opt'] = xyz_dict
            opt_project_info['conformers'][fingerprint]['xyz_str_after_fine_opt'] = xyz_dict_to_xyz_str(xyz_dict)

            opt_project_info['conformers'][fingerprint]['energy']['end_of_fine_opt'] = \
                fingerprint_to_all_opt_log_info_dict[fingerprint]['electronic_energy']['hartree']

            opt_project_info['conformers'][fingerprint]['frequencies'] = \
                fingerprint_to_all_opt_log_info_dict[fingerprint]['freq'][0]
            opt_project_info['conformers'][fingerprint]['negative_frequencies'] = \
                fingerprint_to_all_opt_log_info_dict[fingerprint]['freq'][1]

        valid_conformer_hash_ids = deepcopy(checked_conformer_hash_ids)
        opt_project_info['valid_conformer_hash_ids'] = valid_conformer_hash_ids

    # 2. Generate input files
    # (a) regular opt only, no need for fine opt
    if has_fine_opted is None:
    # 2.a.1 orca single point file
    # 2.a.2. cosmo input file
    # 2.a.3. save project info
        raise NotImplementedError
    # (b) generate fine opt input file
    elif has_fine_opted is False:
    # 2.b.1 Gaussian input file
        conformer_to_fine_opt_hash_ids = deepcopy(checked_conformer_hash_ids)
        opt_project_info['conformer_to_fine_opt_hash_ids'] = conformer_to_fine_opt_hash_ids
        mkdir(opt_fine_dir)

        level_of_theory = opt_project_info['level_of_theory']['fine_opt_freq']
        if level_of_theory.lower() not in ['cbs-qb3']:
            if multiplicity == 2:
                level_of_theory = 'u' + level_of_theory

        for i, fingerprint in enumerate(opt_project_info['conformer_to_fine_opt_hash_ids']):
            fine_opt_input_file_name = str(i) + '_' + str(fingerprint) + '_geom_fine_opt_freq.gjf'
            fine_opt_input_file_path = os.path.join(opt_fine_dir, fine_opt_input_file_name)

            xyz_str = opt_project_info['conformers'][fingerprint]['xyz_str_after_opt']

            fine_opt_input_file_basename = fine_opt_input_file_name.split('.')[0]
            fine_opt_input_file = gen_gaussian_optfreq_input_file(name=fine_opt_input_file_basename,
                                                                  xyz_str=xyz_str,
                                                                  charge=charge,
                                                                  multiplicity=multiplicity,
                                                                  memory_mb=300000,
                                                                  cpu_threads=40,
                                                                  is_ts=is_ts,
                                                                  level_of_theory=level_of_theory,
                                                                  comment=str(fingerprint),
                                                                  )

            opt_project_info['conformers'][fingerprint]['file_path']['input']['fine_opt_freq'] = fine_opt_input_file_path

            with open(fine_opt_input_file_path, 'w') as f:
                f.write(fine_opt_input_file)

    # 2.b.2 Array job submission script
        name = opt_project_info['species']['name']
        last_job_num = len(opt_project_info['conformer_to_fine_opt_hash_ids']) - 1
        sub_script = deepcopy(g16_slurm_array_script)
        sub_script = sub_script.format(last_job_num=str(last_job_num), name=name)
        sub_script_file_path = os.path.join(opt_fine_dir, 'submit_g16_array.sh')
        with open(sub_script_file_path, 'w') as f:
            f.write(sub_script)

    # 2.b.3 save fine opt project info
        fine_opt_project_info_path = os.path.join(opt_fine_dir, 'fine_opt_project_info.yml')
        write_yaml_file(path=fine_opt_project_info_path, content=opt_project_info)

    # (c) generate sp and solvation files using fine opted geometry
    elif has_fine_opted is True:
        # todo: deal with non cbs-qb3 methods
        # todo: deal with other solvation correction methods
        level_of_theory = opt_project_info['level_of_theory']['fine_opt_freq']
        if level_of_theory.lower() not in ['cbs-qb3']:
            raise NotImplementedError
        solvation_method = opt_project_info['level_of_theory']['solv_correction']
        if solvation_method not in ['cosmo']:
            raise NotImplementedError

        valid_conformer_hash_ids = opt_project_info['valid_conformer_hash_ids']

    # 2.c.1 generate sp input file
    # if sp_after_opt is None, use fine opt energy
        sp_level = opt_project_info['level_of_theory'].get('sp_after_opt', None)
        if not sp_level:
            for fingerprint in valid_conformer_hash_ids:
                opt_project_info['conformers'][fingerprint]['energy']['sp_after_opt'] = \
                    opt_project_info['conformers'][fingerprint]['energy']['end_of_fine_opt']
        else:
            raise NotImplementedError

    # 2.c.2 generate solvation input file
    # implemented for cosmo only
    # The following section uses turbomole to generate cosmo file
        cosmo_dir = os.path.join(project_dir, 'cosmo')
        mkdir(cosmo_dir)
        xyz_dir = os.path.join(cosmo_dir, 'xyz')
        mkdir(xyz_dir)

        for i, fingerprint in enumerate(valid_conformer_hash_ids):
            cosmo_input_file_name = str(i) + '_' + str(fingerprint) + '_cosmo.txt'
            cosmo_input_file_path = os.path.join(cosmo_dir, cosmo_input_file_name)
            cosmo_input_file_basename = cosmo_input_file_name.split('.')[0]

            cosmo_input_file = f"{cosmo_input_file_basename} {charge} {multiplicity}"

            with open(cosmo_input_file_path, 'w') as f:
                f.write(cosmo_input_file)

            xyz_str = opt_project_info['conformers'][fingerprint]['xyz_str_after_fine_opt']
            xyz_file = f"{len(xyz_str.splitlines())}\n\n{xyz_str}"
            xyz_file_name = cosmo_input_file_basename + '.xyz'
            xyz_file_name_path = os.path.join(xyz_dir, xyz_file_name)

            with open(xyz_file_name_path, 'w') as f:
                f.write(xyz_file)

    # 2.c.2.2 Generate cosmo-rs input file
            cosmo_rs_input_file_name = str(i) + '_' + str(fingerprint) + '_cosmo.inp'
            cosmo_rs_input_file_path = os.path.join(cosmo_dir, cosmo_rs_input_file_name)
            opt_project_info['conformers'][fingerprint]['file_path']['input']['solv_correction'] = cosmo_rs_input_file_path

            cosmo_rs_input_file_basename = cosmo_rs_input_file_name.split('.')[0]
            cosmo_rs_input_file = gen_cosmo_rs_input_file(name=cosmo_rs_input_file_basename)

            cosmo_rs_output_file_name = str(i) + '_' + str(fingerprint) + '_cosmo.tab'
            cosmo_rs_output_file_path = os.path.join(cosmo_dir, cosmo_rs_output_file_name)
            opt_project_info['conformers'][fingerprint]['file_path']['output']['solv_correction'] = cosmo_rs_output_file_path

            with open(cosmo_rs_input_file_path, 'w') as f:
                f.write(cosmo_rs_input_file)

    # 2.c.2.3 Array job submission script
    # 2.c.2.3.1 for Gaussian cosmo input files
        last_job_num = len(valid_conformer_hash_ids) - 1
        sub_script = deepcopy(turbomole_cosmo_slurm_array_script)
        sub_script = sub_script.format(last_job_num=str(last_job_num))
        sub_script_file_path = os.path.join(cosmo_dir, 'submit_turbomole_array.sh')
        with open(sub_script_file_path, 'w') as f:
            f.write(sub_script)

    # 2.c.2.3.2 for cosmo-rs input files
        cosmo_sub_script = deepcopy(cosmo_slurm_array_script)
        cosmo_sub_script = cosmo_sub_script.format(last_job_num=str(last_job_num))
        cosmo_sub_script_file_path = os.path.join(cosmo_dir, 'submit_cosmo_array.sh')
        with open(cosmo_sub_script_file_path, 'w') as f:
            f.write(cosmo_sub_script)

    # 2.c.2.4 save cosmo project info
        cosmo_project_info_path = os.path.join(cosmo_dir, 'cosmo_project_info.yml')
        write_yaml_file(path=cosmo_project_info_path, content=opt_project_info)



    # # Commented section uses gaussian to generate cosmo file
    # # 2.c.2.1 Generate cosmo gaussian input file
    #     cosmo_dir = os.path.join(project_dir, 'cosmo')
    #     mkdir(cosmo_dir)
    #
    #     for i, fingerprint in enumerate(valid_conformer_hash_ids):
    #         cosmo_input_file_name = str(i) + '_' + str(fingerprint) + '_cosmo.gjf'
    #         cosmo_input_file_path = os.path.join(cosmo_dir, cosmo_input_file_name)
    #
    #         xyz_str = opt_project_info['conformers'][fingerprint]['xyz_str_after_fine_opt']
    #
    #         cosmo_input_file_basename = cosmo_input_file_name.split('.')[0]
    #         cosmo_input_file = gen_gaussian_cosmo_sp_input_file(name=cosmo_input_file_basename,
    #                                                  xyz_str=xyz_str,
    #                                                  charge=charge,
    #                                                  multiplicity=multiplicity,
    #                                                  memory_mb=300000,
    #                                                  cpu_threads=40,
    #                                                  comment=str(fingerprint),
    #                                                  )
    #
    #         opt_project_info['conformers'][fingerprint]['file_path']['input']['solv_correction'] = cosmo_input_file_path
    #
    #         with open(cosmo_input_file_path, 'w') as f:
    #             f.write(cosmo_input_file)
    #
    # # 2.c.2.2 Generate cosmo-rs input file
    #         cosmo_rs_input_file_name = str(i) + '_' + str(fingerprint) + '_cosmo.inp'
    #         cosmo_rs_input_file_path = os.path.join(cosmo_dir, cosmo_rs_input_file_name)
    #
    #         cosmo_rs_input_file_basename = cosmo_rs_input_file_name.split('.')[0]
    #         cosmo_rs_input_file = gen_cosmo_rs_input_file(name=cosmo_rs_input_file_basename)
    #
    #         cosmo_rs_output_file_name = str(i) + '_' + str(fingerprint) + '_cosmo.tab'
    #         cosmo_rs_output_file_path = os.path.join(cosmo_dir, cosmo_rs_output_file_name)
    #         opt_project_info['conformers'][fingerprint]['file_path']['output']['solv_correction'] = cosmo_rs_output_file_path
    #
    #         with open(cosmo_rs_input_file_path, 'w') as f:
    #             f.write(cosmo_rs_input_file)
    #
    # # 2.c.2.3 Array job submission script
    # # 2.c.2.3.1 for Gaussian cosmo input files
    #     name = opt_project_info['species']['name']
    #     last_job_num = len(valid_conformer_hash_ids) - 1
    #     sub_script = deepcopy(g16_slurm_array_script)
    #     sub_script = sub_script.format(last_job_num=str(last_job_num), name=name)
    #     sub_script_file_path = os.path.join(cosmo_dir, 'submit_g16_array.sh')
    #     with open(sub_script_file_path, 'w') as f:
    #         f.write(sub_script)
    #
    # # 2.c.2.3.2 for cosmo-rs input files
    #     cosmo_sub_script = deepcopy(cosmo_slurm_array_script)
    #     cosmo_sub_script = cosmo_sub_script.format(last_job_num=str(last_job_num))
    #     cosmo_sub_script_file_path = os.path.join(cosmo_dir, 'submit_cosmo_array.sh')
    #     with open(cosmo_sub_script_file_path, 'w') as f:
    #         f.write(cosmo_sub_script)
    #
    # # 2.c.2.4 save cosmo project info
    #     cosmo_project_info_path = os.path.join(cosmo_dir, 'cosmo_project_info.yml')
    #     write_yaml_file(path=cosmo_project_info_path, content=opt_project_info)

if __name__ == '__main__':
    main()