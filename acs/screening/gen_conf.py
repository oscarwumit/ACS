"""
This module generates conformers for initial sp screening.
"""

import argparse
import os
import random
import yaml
from itertools import combinations, product
from typing import List, Optional, Union

import numpy as np
from rdkit import Chem
from copy import deepcopy

from acs.converter.geom import (xyz_str_to_xyz_dict,
                                xyz_dict_to_xyz_str,
                                xyz_dict_to_xyz_file)
from acs.script import default_conformer_info_dict_for_initial_sp_screening, \
    default_job_info_dict_for_initial_sp_screening
from acs.common import read_yaml_file, write_yaml_file
from openbabel import pybel

from rdmc.mol import RDKitMol, RDKitTS
from rdmc.ts import get_formed_and_broken_bonds

def get_separable_angle_list(conf,
                             torsions,
                             samplings: Union[list, tuple],
                             from_angles: Union[list, tuple] = None):
    """
    Get a angle list for each input torsional dimension. For each dimension
    The input can be a int, indicating the angles will be evenly sampled;
    Or a list, indicate the angles to be sampled. You can also input a
    ``from_angles`` to customize the starting angles.
    Examples for the ``sampling``:
    [[120, 240,], 4, 0] => [np.array([120, 240,]),
                            np.array([0, 90, 180, 270,]),
                            np.array([0])]

    Args:
        conf (RDKitConf): The conformer
        samplings (Union[list, tuple]): An array of sampling information.
                  For each element, it can be either list or int.
        from_angles (Union[list, tuple]): An array of initial angles.
                    If not set, angles will begin at orignal dihedral angles.

    Returns:
        list: A list of sampled angles sets.
    """
    from_angles = list() if from_angles is None else from_angles

    if not from_angles:  # Assuming angles are rotating from the original value
        for tor in torsions:
            from_angles.append(conf.GetTorsionDeg(tor))

    angle_list = []
    for ind, angles in enumerate(samplings):
        if isinstance(angles, (int, float)):
            # Only provide a number
            # This is the step number of the angles
            try:
                step = 360 // angles
            except ZeroDivisionError:
                # Does not change
                angles = from_angles[ind] + np.array([0])
            else:
                angles = from_angles[ind] + \
                         np.array([step * i for i in range(angles)])
        elif isinstance(angles, list):
            angles = from_angles[ind] + np.array(angles)

        # Set to 0 - 360 range
        for i in range(angles.shape[0]):
            while angles[i] < 0.:
                angles[i] += 360
            while angles[i] > 360.:
                angles[i] -= 360

        angle_list.append(angles.tolist())
    return angle_list


def conformers_by_change_torsions(conf: 'RDKitConf',
                                  angle_mesh,
                                  bookkeep: dict,
                                  torsions = None,
                                  on_the_fly_check = False):
    """
    Generate conformers by rotating the angles of the torsions. The result will be saved into
    ``bookkeep``. A on-the-fly check can be applied, which identifies the conformers with colliding
    atoms.

    Args:
        conf (RDkitConf): A RDKit Conformer to be used.
        angle_mesh (iterable): An iterable contains the angle_list for conformers to be generated from.
        bookkeep (dict): A dictionary to save the coords.
        torsions (list): A list of four-atom-index lists indicating the torsional modes.
        on_the_fly_filter (bool): Whether to check colliding atoms on the fly.
    """
    # todo: extend this function to accomodate different dimensionality
    if not torsions:
        # Torsions are not set, assuming changing all of the torsions
        torsions = conf.GetTorsionalModes()

    ref_xyz_dict = bookkeep['species']['coord']['arc_xyz']
    n_torsions = len(torsions)
    if n_torsions not in (0, 1, 2):
        raise NotImplementedError

    if not n_torsions == 0:
        index_1 = -1
        index_2 = 0
        lookup = set()
        for angles in angle_mesh:

            if angles[0] not in lookup:
                index_1 += 1
                index_2 = 0
                lookup.add(angles[0])

            hash_id = random.getrandbits(128)
            hash_key = hex(hash_id)

            bookkeep['conformers'][hash_key] = {'rotor_dimension': n_torsions,
                                                'dihedral': list(), }

            for angle, tor in zip(angles, torsions):
                conf.SetTorsionDeg(tor, angle)
                bookkeep['conformers'][hash_key]['dihedral'].append([tuple(tor), angle])

            bookkeep['conformers'][hash_key]['dihedral'][0].append(index_1)
            if n_torsions == 2:
                bookkeep['conformers'][hash_key]['dihedral'][1].append(index_2)

            bookkeep['conformers'][hash_key]['dihedral'][0] = tuple(bookkeep['conformers'][hash_key]['dihedral'][0])
            if n_torsions == 2:
                bookkeep['conformers'][hash_key]['dihedral'][1] = tuple(bookkeep['conformers'][hash_key]['dihedral'][1])

            bookkeep['conformers'][hash_key]['dihedral'] = tuple(bookkeep['conformers'][hash_key]['dihedral'])
            index_2 += 1

            xyz_dict = deepcopy(ref_xyz_dict)
            xyz_list = conf.GetPositions().tolist()
            xyz_dict['coords'] = xyz_list
            xyz_str = xyz_dict_to_xyz_str(xyz_dict)
            bookkeep['conformers'][hash_key].update({
                'xyz_str': xyz_str,
                'arc_xyz': xyz_dict,
                'is_colliding': conf.HasCollidingAtoms() if on_the_fly_check else None,
                'is_crashing': None,
            })
    else:
        hash_id = random.getrandbits(128)
        hash_key = hex(hash_id)
        bookkeep['conformers'][hash_key] = {'rotor_dimension': n_torsions,
                                            'dihedral': list(), }
        xyz_dict = deepcopy(ref_xyz_dict)
        xyz_list = conf.GetPositions().tolist()
        xyz_dict['coords'] = xyz_list
        xyz_str = xyz_dict_to_xyz_str(xyz_dict)
        bookkeep['conformers'][hash_key].update({
            'xyz_str': xyz_str,
            'arc_xyz': xyz_dict,
            'is_colliding': conf.HasCollidingAtoms() if on_the_fly_check else None,
            'is_crashing': None,
        })


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
    The gen_conf executable function
    """
    # todo: write log
    # todo: exceptions
    # todo: ts volume check
    # 0. Parse input
    # 0.1 Basic project info
    args = parse_command_line_arguments()
    input_file = args.file
    project_directory = os.path.abspath(os.path.dirname(args.file))
    input_dict = read_yaml_file(path=input_file)

    project_info = deepcopy(default_job_info_dict_for_initial_sp_screening)
    project_info.update(input_dict)
    project_info['project_folder_path'] = project_directory

    is_ts = project_info['species']['is_ts']

    # 0.2 Load geometry from input file
    # todo: deal with other coord inputs
    # here assume xyz_str is always there, obviously we need to improve this
    xyz_str =  project_info['species']['coord']['xyz_str']
    xyz_dict = xyz_str_to_xyz_dict(xyz_str)
    project_info['species']['coord']['arc_xyz'] = xyz_dict
    xyz_file = f"{len(xyz_str.splitlines())}\n\n{xyz_str}"

    # 1. Generate conformers
    # 1.1.a Process TS conformer
    if is_ts:
        # 1.1.a.1  Perceive TS
        rdkitts = RDKitTS.FromOBMol(pybel_mol.OBMol)
        rdkitts.EmbedConformer()
        conf = rdkitts.GetConformer()
        conf.SetPositions(xyz_dict['coords'])

        # 1.1.a.2 Set the missing bonds existing in the TS
        bonds = list()
        threshold = 1.6

        if not bonds:
            dist_mat = np.triu(rdkitts.GetDistanceMatrix())
            covl_mat = rdkitts.GetCovalentMatrix()

            for multiplier in np.arange(1.1, threshold, 0.1):
                atom1s, atom2s = np.where((dist_mat - multiplier * covl_mat) < 0)
                bonds = [(int(atom1s[i]), int(atom2s[i])) for i in range(len(atom1s))]
            if not bonds:
                raise ValueError('Cannot id TS bond!')

        # 1.1.a.3 Overwrite the RDKitTS with new bonding info
        rw_mol = rdkitts.ToRDMol()

        for bond in bonds:
            # Use BondType.OTHER if you want to avoid to be counted as a torsional mode
            # If you want to include it, please use BondType.SINGLE
            rw_mol.AddBond(*bond, Chem.BondType.OTHER)

        rdkitts = rdkitts.FromRDMol(rw_mol)
    # 1.1.b Process non-TS conformer
    # no need to guess connectivity from OpenBabel. Instead, use the SMILES for stable species
    else:
        smi = project_info['species']['smiles']
        rdkitmol = RDKitMol.FromSmiles(smi)

    # 1.2 Use RDKit to generate conformers
    # 1.2.1.a Initialize a TS conformer instance
    if is_ts:
        rdkitts.EmbedConformer()
        conf = rdkitts.GetConformer()
    # 1.2.1.b Initialize a non-TS conformer instance
    else:
        rdkitmol.EmbedConformer()
        rdkitmol.SetPositions(xyz_dict['coords'])
        conf = rdkitmol.GetConformer()

    # 1.2.2 Get the torsional mode and the original angles
    # You can set the correct (all) torsions, otherwise RDKit will perceive.
    ######################################
    # INPUT
    torsions = None
    ######################################
    if not torsions:
        if is_ts:
            torsions = rdkitts.GetTorsionalModes()
        else:
            torsions = rdkitmol.GetTorsionalModes(excludeMethyl=False)
        # print(f'RDKit perceived torsions: {torsions}')

    conf.SetTorsionalModes(torsions)
    num_torsions = len(torsions)
    original_angles = conf.GetAllTorsionsDeg()
    # # print('Torsions highlighted in the molecule:')
    # display(rdkitts)
    # print(f'The original dihedral angles is: {original_angles}')

    # Save to dict
    project_info['species']['1d_torsions'] = torsions
    if is_ts:
        project_info['species']['coord']['connectivity'] = tuple([sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
                                                                  for bond in rdkitts.GetBonds()])
    else:
        project_info['species']['coord']['connectivity'] = tuple([sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
                                                                  for bond in rdkitmol.GetBonds()])

    # 1.2.3 Generate conformers according to the mangle mesh
    # todo: make this a function
    n_point_each_torsion = project_info.get('n_point_each_torsion', 20)
    project_info['n_point_each_torsion'] = n_point_each_torsion
    n_dimension = project_info.get('n_rotors_to_couple', 2)
    project_info['n_rotors_to_couple'] = n_dimension

    bookkeeps = []

    torsion_pairs = list(combinations(torsions, n_dimension))

    for torsion_pair in torsion_pairs:
        # Reset the geometry
        conf.SetPositions(xyz_dict['coords'])
        # Get angles
        sampling = [n_point_each_torsion for i in range(n_dimension)]
        angles_list = get_separable_angle_list(conf,
                                               torsion_pair,
                                               sampling, )
        angle_mesh = product(*angles_list)
        # Generate conformers
        project_info_copy = deepcopy(project_info)
        conformers_by_change_torsions(conf,
                                      angle_mesh,
                                      bookkeep=project_info_copy,
                                      torsions=torsion_pair,
                                      on_the_fly_check=True)
        project_info_copy['dihedrals_considered_in_this_file'] = torsion_pair
        bookkeeps.append(project_info_copy)

    # 1.4 Save to files
    os.mkdir(os.path.join(project_directory, 'initial_sp_screening'))
    sub_folder_info = dict()
    for i, bookkeep in enumerate(bookkeeps):
        sub_folder_info[i] = bookkeep['dihedrals_considered_in_this_file']
        conformer_to_screen_hash_ids = list()
        colliding_conformer_hash_ids = list()
        for k, v in bookkeep['conformers'].items():
            if not bookkeep['conformers'][k]['is_colliding']:
                conformer_to_screen_hash_ids.append(k)
            else:
                colliding_conformer_hash_ids.append(k)

        bookkeep['conformer_to_screen_hash_ids'] = tuple(conformer_to_screen_hash_ids)
        bookkeep['colliding_conformer_hash_ids'] = tuple(colliding_conformer_hash_ids)

        sub_dir_path = os.path.join(project_directory, 'initial_sp_screening', str(i))
        os.mkdir(sub_dir_path)
        outfile_path = os.path.join(sub_dir_path, 'initial_conf_coords.yml')
        write_yaml_file(path=outfile_path, content=bookkeep)
    sub_folder_info_path = os.path.join(project_directory, 'initial_sp_screening', 'sub_folder_info.yml')
    write_yaml_file(path=sub_folder_info_path, content=sub_folder_info)



if __name__ == '__main__':
    main()
