"""
This module analyzes initial sp screening results.
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
from copy import deepcopy

from acs.common import read_yaml_file, write_yaml_file, mkdir, update_selected_keys, gen_gaussian_input_file
from acs.script import default_job_info_dict_after_initial_sp_screening, g16_slurm_array_script, \
    default_conformer_info_dict_after_initial_sp_screening


def get_step_to_adjacent_points(fsize, dim=2, cutoff=np.inf):
    one_d_points = list(range(- fsize, fsize + 1))
    var_combinations = product(*[one_d_points] * dim)
    for points in var_combinations:
        dist = np.linalg.norm(np.array(points))
        if dist <= cutoff:
            yield points


def get_energy(coord, energies):
    try:
        return energies[coord]
    except IndexError:
        new_coord = tuple(x if x < energies.shape[i] else
                          x - energies.shape[i] for i, x in enumerate(coord))
        return energies[new_coord]


def compare_to_adjacent_point(coord, energies, unchecked_points, filters):
    # each element is a coordinate
    new_coords = [tuple(x + var_x for x, var_x in zip(coord, var))
                  for var in filters]

    # Get the energies of adjacent points
    energies = [get_energy(new_coord, energies) for new_coord in new_coords]

    # Sorted
    energies, new_coords = zip(*sorted(zip(energies, new_coords)))

    # Find the current point index and points that has higher energy than this point
    # Will be removed from unchecked points list
    cur_point_ind = new_coords.index(coord)
    for new_coord in new_coords[cur_point_ind:]:
        try:
            unchecked_points.remove(new_coord)
        except ValueError:
            # ValueError if coord_min is not in unchecked_points
            pass
    return new_coords[0]


def search_for_a_minimum(coord, energies, unchecked_points, filters):
    while True:
        next_point = compare_to_adjacent_point(coord, energies,
                                               unchecked_points, filters)
        next_point = tuple(x if x >= 0 else energies.shape[i] + x
                           for i, x in enumerate(next_point))
        if next_point == coord:
            return coord
        elif next_point not in unchecked_points:
            return
        else:
            coord = next_point


def search_minimum(energies, fsize, cutoff=np.inf):
    minimum = []

    dim = len(energies.shape)
    filters = list(get_step_to_adjacent_points(fsize, dim, cutoff))

    oned_points = [list(range(energies.shape[i])) for i in range(dim)]
    unchecked_points = list(product(*oned_points))

    while True:
        if not unchecked_points:
            break
        coord = unchecked_points[np.random.randint(len(unchecked_points))]
        new_min = search_for_a_minimum(coord, energies,
                                       unchecked_points, filters)
        if new_min:
            minimum.append(new_min)
    return minimum

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
    # todo: parse opt info from input file
    # todo: ML approach to explore higher dimensions
    """
    The analyze_screen executable function
    """
    # 0. Parse input
    # 0.1 Parse screening project info
    args = parse_command_line_arguments()
    screening_input_file = args.file
    # input_file_dir = os.path.abspath(os.path.dirname(args.file)) # todo: use this to decude abs path if folder moved
    screen_project_info = read_yaml_file(screening_input_file)
    n_rotors_to_couple = screen_project_info.get('n_rotors_to_couple', 2)
    project_dir = screen_project_info['project_folder_path']
    screen_dir = os.path.join(project_dir, 'initial_sp_screening')
    # 0.2 Initialize optimization project info
    opt_dir = os.path.join(project_dir, 'opt')
    mkdir(opt_dir)
    opt_project_info = deepcopy(default_job_info_dict_after_initial_sp_screening)
    opt_project_info = update_selected_keys(opt_project_info,
                                            screen_project_info,
                                            keys_to_update=('project',
                                                            'project_folder_path',
                                                            'species',
                                                            'level_of_theory',
                                                            ))

    # 1. Analyze screening result
    # pre-allocation
    conformer_to_opt_hash_ids = list()

    # Process results one by one
    sub_folder_info = read_yaml_file(os.path.join(screen_dir, 'sub_folder_info.yml'))
    for sub in sub_folder_info.keys():
        subfolder = str(sub)
        screen_result_file_path = os.path.join(screen_dir, subfolder, 'initial_conf_screening_result.yml')
        screen_result = read_yaml_file(screen_result_file_path)

    # 1.1. Populate energy array for further analyses
    # 1.1.1 Pre-allocation based on num of torsions and sampling points
        num_torsions = len(screen_result['dihedrals_considered_in_this_file'])
        num_sampling = screen_result['n_point_each_torsion']
        energies = np.empty([num_sampling for i in range(num_torsions)])


    # 1.1.2 Record energies
        index_to_id = dict()
        collide = screen_result.get('colliding_conformer_hash_ids')
        collide = collide if collide is not None else tuple()
        crash = screen_result.get('crashing_conformer_hash_ids')
        crash = crash if crash is not None else tuple()

        for fingerprint in screen_result['conformers'].keys():
            dihedral = screen_result['conformers'][fingerprint]['dihedral']
            index = tuple([d[-1] for d in dihedral])
            # todo: better treat different dimensions
            if n_rotors_to_couple == 1:
                index_1 = tuple([[d[-1] for d in dihedral][0], 0])
                index_to_id[index_1] = fingerprint
            else:
                index_to_id[index] = fingerprint

            if fingerprint not in (collide + crash):
                energies[index] = screen_result['conformers'][fingerprint]['initial_screening_sp_energy']
            elif fingerprint in collide:
                energies[index] = None
            elif fingerprint in crash:
                energies[index] = np.nan

        # todo: make plots using mask
        # mask = np.zeros_like(energies)
        # mask[energies == None] = 1

    # 1.1.3 For colliding conformers, set the energy (None) to be very high
        max_energy = np.nanmax(energies)
        replacer_threshold = 0.01
        replacer = (1 - replacer_threshold) * max_energy if max_energy < 0 else (1 + replacer_threshold) * max_energy
        # noinspection PyComparisonWithNone
        energies[energies == None] = replacer  # energies is None does not work for numpy in this case, do not replace!

    # 1.1.4 For crashing conformers, set the energy (np.nan) via linear interpolation
        # todo: save dataframe, also deal with non-even sampling
        # todo: better treat different dimensions
        if not n_rotors_to_couple == 0:
            label = [str(int(n * 360 / num_sampling)) for n in range(num_sampling)]
            if n_rotors_to_couple == 2:
                df = pd.DataFrame(energies, index=label, columns=label)
            elif n_rotors_to_couple == 1:
                df = pd.DataFrame(energies, index=label, columns=['0'])
            else:
                raise NotImplementedError
            df = df.interpolate(method='linear', axis=0, limit_direction='both').interpolate(method='linear', axis=1,
                                                                                             limit_direction='both')
            energies = df.to_numpy()

        # 1.1.5 Rescale energies by minimum
            # This will not change the result of search but make detailed view more clear
            energies = energies - np.min(energies)

        # 2. Search for minima in energy array
            # todo: alternative method via function fitting and optimization
        # 2.1 Greedy search for minima in array
            # todo: better treat different dimensions, the output when fsize=1 should not be a tuple of two indices
            minimum_points = search_minimum(energies, fsize=n_rotors_to_couple)
        # 2.2 Set an energy ceiling
            ceiling = 1.0  # 1 hartree = 627.509 kcal/mol, it seems high but the geom is not optimized
            minimum_points = [i for i in minimum_points if energies[i] < ceiling]
        else:
            minimum_points = list(index_to_id.keys())

    # 3. Save minima to optimization project dictionary
        for index in minimum_points:
            fingerprint = index_to_id[index]
            conformer_to_opt_hash_ids.append(fingerprint)

            conformer_from_screen = screen_result['conformers'][fingerprint]
            conformer_to_opt = deepcopy(default_conformer_info_dict_after_initial_sp_screening)

            conformer_to_opt['rotor_dimension'] = conformer_from_screen['rotor_dimension']
            conformer_to_opt['dihedral_before_opt'] = conformer_from_screen['dihedral']
            conformer_to_opt['xyz_str_before_opt'] = conformer_from_screen['xyz_str']
            conformer_to_opt['arc_xyz_before_opt'] = conformer_from_screen['arc_xyz']
            conformer_to_opt['energy']['initial_screening_sp'] = conformer_from_screen['initial_screening_sp_energy']
            conformer_to_opt['file_path']['input']['initial_screening_sp'] = \
                os.path.join(screen_dir, subfolder, 'initial_conf_coords.yml')
            conformer_to_opt['file_path']['output']['initial_screening_sp'] = screen_result_file_path

            opt_project_info['conformers'][fingerprint] = conformer_to_opt

    opt_project_info['conformer_to_opt_hash_ids'] = tuple(conformer_to_opt_hash_ids)

    # 4. Generate opt input file
    # 4.1 Gaussian input file
    charge = opt_project_info['species']['charge']
    multiplicity = opt_project_info['species']['multiplicity']
    is_ts = opt_project_info['species']['is_ts']
    # todo: deal with other levels and multiplicity = 3
    # assume multiplicity = 1 or 2 here
    level_of_theory = opt_project_info['level_of_theory']['opt_freq']
    if level_of_theory.lower() not in ['cbs-qb3']:
        if multiplicity == 2:
            level_of_theory = 'u' + level_of_theory

    for i, fingerprint in enumerate(opt_project_info['conformer_to_opt_hash_ids']):
        opt_input_file_name = str(i) + '_' + str(fingerprint) + '_geom_opt_freq.gjf'
        opt_input_file_path = os.path.join(opt_dir, opt_input_file_name)

        xyz_str = opt_project_info['conformers'][fingerprint]['xyz_str_before_opt']

        opt_input_file_basename = opt_input_file_name.split('.')[0]
        opt_input_file = gen_gaussian_input_file(name=opt_input_file_basename ,
                                                    xyz_str=xyz_str,
                                                    charge=charge,
                                                    multiplicity=multiplicity,
                                                    memory_mb=300000,
                                                    cpu_threads=40,
                                                    is_ts=is_ts,
                                                    level_of_theory=level_of_theory,
                                                    comment=str(fingerprint),
                                                    )

        opt_project_info['conformers'][fingerprint]['file_path']['input']['opt_freq'] = opt_input_file_path

        with open(opt_input_file_path, 'w') as f:
            f.write(opt_input_file)

    # 4.2 Array job submission script
    # todo: make submission script more general
    name = opt_project_info['species']['name']
    last_job_num = len(opt_project_info['conformer_to_opt_hash_ids']) - 1
    sub_script = deepcopy(g16_slurm_array_script)
    sub_script = sub_script.format(last_job_num=str(last_job_num), name=name)
    sub_script_file_path = os.path.join(opt_dir, 'submit_g16_array.sh')
    with open(sub_script_file_path, 'w') as f:
        f.write(sub_script)


    # 5. save opt project info
    opt_project_info_path = os.path.join(opt_dir, 'opt_project_info.yml')
    write_yaml_file(path=opt_project_info_path, content=opt_project_info)



if __name__ == '__main__':
    main()


