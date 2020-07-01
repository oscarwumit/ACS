"""
This module contains functions which are shared across multiple ACS modules.
"""


import os

import yaml
from typing import Any, List, Optional, Tuple, Union, Dict
import shutil

from copy import deepcopy
import numpy as np
import pandas as pd

import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
from scipy import interpolate
from scipy.optimize import minimize

# from arc.species.species import ARCSpecies
# from arc.species.converter import xyz_to_xyz_file_format, str_to_zmat, zmat_to_xyz, modify_coords


def read_yaml_file(path: str,
                   ) -> dict:
    """
    Read a YAML file (usually an input / restart file, but also conformers file)
    and return the parameters as python variables.

    Args:
        path (str): The YAML file path to read.

    Returns:
        dict or list: The content read from the file.
    """
    with open(path, 'r') as f:
        content = yaml.load(stream=f, Loader=yaml.Loader)
    return content


def write_yaml_file(path: str,
                    content: dict,
                   ) -> None:
    """
    Write a YAML file (usually an input / restart file, but also conformers file)
    and return the parameters as python variables.

    Args:
        path (str): The YAML file path to write.
        content (dict): Information to write.

    Returns:
        dict or list: The content read from the file.
    """
    with open(path, 'w') as f:
        yaml.dump(content, f, default_flow_style=False)

def mkdir(path: str,
          force: bool = True,
          ) -> None:
    if force and os.path.exists(path):
        shutil.rmtree(path)

    os.mkdir(path)

def update_selected_keys(dict1: dict,
                         dict2: dict,
                         keys_to_ignore: tuple = tuple(),
                         keys_to_update: tuple = tuple(),
                         ) -> dict:

    dict1_cp = deepcopy(dict1)
    s = set(dict2.keys())
    t1 = set(keys_to_ignore)
    t2 = set(keys_to_update)
    s.difference_update(t1)
    s.intersection_update(t2)

    for k in s:
        if k in dict1_cp:
            dict1_cp[k] = dict2[k]
    return dict1_cp

def gen_gaussian_input_file(name: str,
                            xyz_str: str,
                            charge: int,
                            multiplicity: int,
                            memory_mb: int,
                            cpu_threads: int,
                            is_ts: bool,
                            level_of_theory: str,
                            comment: str = '',
                            ) -> str:
    if is_ts:
        title_card = "#p opt=(calcfc,noeigentest,maxcycles=120) freq guess=mix scf=xqc iop(2/9=2000)"
    else:
        title_card = "#p opt=(calcfc,ts,noeigentest,maxcycles=120) freq guess=mix scf=xqc iop(2/9=2000)"

    script = f"""%chk={name}.chk
%mem={memory_mb}mb
%NProcShared={cpu_threads}

{title_card} {level_of_theory}

{comment}

{charge} {multiplicity}
{xyz_str}



"""
    return script

# def convert_gaussian_zmat_to_arc_zmat(zmat_file_path: str,
#                                       ) -> Dict:
#     """
#     Convert Gaussian Z-matrix to ARC Z-matrix.
#
#     Args:
#         zmat_file_path: Gaussian input file with Z-matrix.
#
#     Returns:
#         ARC Z-matrix dictionary.
#     """
#
#     # read gaussian input
#     atom_dict = dict()
#     hydrogen_list = list()
#     connectivity_dict = dict()
#
#     with open(zmat_file_path, 'r') as f:
#         index = 1
#         line = f.readline()
#         flag0 = True
#         flag1 = False
#         flag2 = False
#         while line:
#             if flag0:
#                 try:
#                     if not line.split()[0] == "0":
#                         line = f.readline()
#                     else:
#                         flag0 = False
#                         flag1 = True
#                         line = f.readline()
#                 except IndexError:
#                     line = f.readline()
#
#             if flag1:
#                 line_content = line.split()
#                 atom_dict[index] = dict()
#                 atom_dict[index].update(atom=line_content[0])
#                 if line_content[0] == 'H':
#                     hydrogen_list.append(index)
#                 if index > 1:
#                     atom_dict[index].update(bond=(line_content[1], line_content[2]))
#                 if index > 2:
#                     atom_dict[index].update(angle=(line_content[3], line_content[4]))
#                 if index > 3:
#                     atom_dict[index].update(dihedral=(line_content[5], line_content[6]))
#                 line = f.readline()
#                 index += 1
#                 if not line.split():
#                     flag1 = False
#                     flag2 = True
#                     line = f.readline()
#             elif flag2:
#                 line_content = line.split()
#                 if not line_content:
#                     break
#                 key = line_content[0]
#                 val = line_content[1]
#                 connectivity_dict[key] = val
#                 line = f.readline()
#
#     # convert to ARC zmat
#     arc_zmat = dict()
#     symbols_list = list()
#     coords_list = list()
#     vars_dict = dict()
#     map_dict = dict()
#
#     for atom_id in atom_dict.keys():
#         atom_num = atom_id - 1
#         map_dict[atom_num] = atom_num
#
#         atom_id_dict = atom_dict[atom_id]
#
#         symbols_list.append(atom_id_dict['atom'])
#
#         bond_info_tuple = atom_id_dict.get('bond', None)
#         angle_info_tuple = atom_id_dict.get('angle', None)
#         dihedral_info_tuple = atom_id_dict.get('dihedral', None)
#
#         R = None
#         A = None
#         D = None
#         r = None
#         a = None
#
#         if bond_info_tuple is not None:
#             r = str(int(bond_info_tuple[0]) - 1)
#             R = '_'.join(['R', str(atom_num), r])
#             vars_dict[R] = float(connectivity_dict[bond_info_tuple[1]])
#
#         if angle_info_tuple is not None:
#             a = str(int(angle_info_tuple[0]) - 1)
#             A = '_'.join(['A', str(atom_num), r, a])
#             vars_dict[A] = float(connectivity_dict[angle_info_tuple[1]])
#
#         if dihedral_info_tuple is not None:
#             d = str(int(dihedral_info_tuple[0]) - 1)
#             D = '_'.join(['D', str(atom_num), r, a, d])
#             vars_dict[D] = float(connectivity_dict[dihedral_info_tuple[1]])
#
#         coords_list.append((R, A, D))
#
#     arc_zmat['symbols'] = tuple(symbols_list)
#     arc_zmat['coords'] = tuple(coords_list)
#     arc_zmat['vars'] = vars_dict
#     arc_zmat['map'] = map_dict
#
#     return arc_zmat
#
#
# def composite_modify_coords(scan_res: int,
#                             arc_zmat: Dict,
#                             pivots: Tuple[Tuple[int], Tuple[int]],
#                             arc_spc: ARCSpecies,
#                             ) ->  dict:
#     """
#     Modify the coordinates two times.
#
#     Args:
#         scan_res: Number of sampling points for each 1d rotor.
#         arc_zmat: ARC's Z-matrix dictionary.
#         pivots: Atom number representing pivots in the molecule.
#         arc_spc: ARC species.
#
#     Returns:
#         A dictionary with keys represent dihedral combination and values represent xyz.
#     """
#
#     scan_deg = 360 / scan_res
#     xyz_new_dict = dict()
#
#     for i in range(scan_res):
#         zmat_1 = modify_coords(coords=arc_zmat,
#                                indices=pivots[0],
#                                new_value=i * scan_deg,
#                                modification_type='groups',
#                                mol=arc_spc.mol,
#                                index=1,
#                                output_zmat=True,
#                                )
#
#         for j in range(scan_res):
#             zmat_2 = modify_coords(coords=zmat_1,
#                                    indices=pivots[1],
#                                    new_value=j * scan_deg,
#                                    modification_type='groups',
#                                    mol=arc_spc.mol,
#                                    index=1,
#                                    output_zmat=True,
#                                    )
#
#             xyz_new_dict[(i, j)] = xyz_to_xyz_file_format(zmat_to_xyz(zmat_2))
#     return xyz_new_dict
#
#
# def process_gaussian_opt_freq_output(logfile):
#     if not check_gaussian_normal_termination(logfile):
#         raise ParserError('Gaussian error termination.')
#     info = dict()
#     info['freq'] = get_gaussian_freq(logfile, checkneg=True, ts=False)
#     info['xyz'] = get_gaussian_geometry(logfile, plot=False)
#     info['energy'] = get_gaussian_energy(logfile)
#     return info
#
#
# def get_gaussian_freq(logfile, checkneg=True, ts=True):
#     freq = parse_frequencies(logfile, software='gaussian')
#     neg_freq = tuple([float(x) for x in freq if x < 0])
#     if checkneg:
#         if ts:
#             if len(neg_freq) == 0:
#                 raise ParserError('Did not find any negative frequencies.')
#             elif len(neg_freq) > 1:
#                 raise ParserError(f'Find more than one negative frequencies: {neg_freq}')
#         else:
#             if len(neg_freq):
#                 raise ParserError(f'Find negative frequencies for non-TS species: {neg_freq}')
#     return (freq, neg_freq)
#
#
# def check_gaussian_normal_termination(logfile):
#     with open(logfile, 'r') as f:
#         lines = f.readlines()
#         forward_lines = tuple(lines)
#     for line in forward_lines[-1:-20:-1]:
#         if 'Normal termination' in line:
#             return True
#     else:
#         return False
#
#
# def get_gaussian_energy(logfile):
#     energy_dict = dict()
#     e_j_mol = parse_e_elect(logfile)
#     energy_dict['J/mol'] = e_j_mol
#     e_kj_mol = e_j_mol / 1000
#     energy_dict['kJ/mol'] = e_kj_mol
#     e_kcal_mol = e_j_mol / 4184
#     energy_dict['kcal/mol'] = e_kcal_mol
#     e_scf = round(e_j_mol/(constants.E_h * constants.Na / 1000), 9)
#     energy_dict['scf'] = e_scf
#     return energy_dict
#
#
# def get_gaussian_geometry(logfile, plot=False):
#     xyz = parse_geometry(logfile)
#     if plot:
#         show_sticks(xyz)
#     return xyz
#
#
# def write_gaussian_input_file():
#     script = """%chk={name}.chk
#     %mem=300000mb
#     %NProcShared=40
#
#     #p opt=(calcfc,noeigentest,maxcycles=120) freq guess=mix uwb97xd def2svp iop(2/9=2000) scf=xqc
#
#     {name}
#
#     0 2
#     {xyz}
#
#
#
#
#     """
#
#     all_xyz_dict = dict()
#
#     fingerprint = 0
#     scan_pts = 45
#     scan_deg = 360 / scan_pts
#
#     for grid_search_file_name in os.listdir(psi4_scan_dir):
#
#         indices_1 = tuple(
#             [int(x) for x in re.search('_oo_(.*)_n_', grid_search_file_name).group(1).split('_') if x.isnumeric()])
#         indices_2 = tuple(
#             [int(x) for x in re.search('_n_(.*)_coord', grid_search_file_name).group(1).split('_') if x.isnumeric()])
#
#         print('----------------------------')
#         print(f'Considering dihedral combinations: {indices_1} and {indices_2}')
#
#         with open(psi4_scan_dir + '/' + grid_search_file_name, 'r') as outfile:
#             energy = yaml.load(outfile, Loader=yaml.FullLoader)
#
#         a = np.zeros((scan_pts, scan_pts))
#         for k in energy.keys():
#             a[k[0], k[1]] = energy[k]
#
#         label = [str(n * scan_deg) for n in range(scan_pts)]
#         df = pd.DataFrame(a, index=label, columns=label)
#
#         df = df.interpolate(method='linear', axis=0, limit_direction='both').interpolate(method='linear', axis=1,
#                                                                                          limit_direction='both')
#
#         g = interpolate.RectBivariateSpline(range(scan_pts), range(scan_pts), df)
#
#         local_minima_locations = detect_local_minima(df)
#
#         x0s = (np.array([x, y]) for x, y in zip(local_minima_locations[0], local_minima_locations[1]))
#
#         res_list = list()
#         for x0 in x0s:
#             res = minimize(run_2d_params, x0=x0, args=g, method='Nelder-Mead', tol=1e-12)
#             res_list.append(res)
#         res_tuple = tuple(res_list)
#
#         res_result = tuple([(r.x[0], r.x[1], r.fun) for r in res_tuple])
#         print('Fitted local minima')
#         print(res_result)
#
#         for r in res_result:
#             new_val_1 = r[0]
#             new_val_2 = r[1]
#
#             xyz_1_new_dihedral = new_val_1 * scan_deg
#             xyz_2_new_dihedral = new_val_2 * scan_deg
#
#             xyz_1 = modify_coords(coords=arc_zmat,
#                                   indices=indices_1,
#                                   new_value=xyz_1_new_dihedral,
#                                   modification_type='groups',
#                                   mol=spc.mol,
#                                   index=1,
#                                   output_zmat=True,
#                                   )
#
#             xyz_2 = modify_coords(coords=xyz_1,
#                                   indices=indices_2,
#                                   new_value=xyz_2_new_dihedral,
#                                   modification_type='groups',
#                                   mol=spc.mol,
#                                   index=1,
#                                   output_zmat=True,
#                                   )
#
#             fingerprint += 1
#
#             all_xyz_dict[(indices_1, int(round(xyz_1_new_dihedral)), indices_2, int(round(xyz_2_new_dihedral)),
#                           fingerprint)] = zmat_to_xyz(xyz_2)
#
#             non_colliding_xyz = [xyz for xyz in tuple(all_xyz_dict.values()) if
#                                  not colliding_atoms(xyz=xyz, threshold=0.65)]
#
#             isomorphic_xyz = list()
#             for xyz in tuple(non_colliding_xyz):
#                 try:
#                     spc_to_check = ARCSpecies(label='check', xyz=xyz, is_ts=False, multiplicity=2,
#                                               smiles='NCCC(O[O])N(C)C')
#                     spc_to_check.mol_from_xyz(xyz)
#
#                     if check_isomorphism(spc.mol, spc_to_check.mol):
#                         isomorphic_xyz.append(xyz)
#                 except:
#                     continue
#
#             all_xyz_distinct_tuple = cluster_confs_by_rmsd(tuple(isomorphic_xyz))
#
#             file_counter = 0
#             save_batch_size = 200000
#             batch_folder_counter = 1
#
#             for distinct_xyz in all_xyz_distinct_tuple:
#
#                 k = key_by_val(all_xyz_dict, distinct_xyz)
#
#                 indices_1 = k[0]
#                 deg_1 = k[1]
#
#                 indices_2 = k[2]
#                 deg_2 = k[3]
#
#                 if not file_counter % save_batch_size:
#                     batch_foler = 'batch_' + str(batch_folder_counter)
#                     if not os.path.exists(save_folder + '/' + batch_foler):
#                         os.mkdir(save_folder + '/' + batch_foler)
#                         batch_folder_counter += 1
#
#                 file_counter += 1
#
#                 xyz_str = xyz_to_xyz_file_format(distinct_xyz)
#                 xyz_str = '\n'.join(xyz_str.split('\n')[2:-1])
#
#                 d1str = "{0:.4g}".format(deg_1)
#                 d2str = "{0:.4g}".format(deg_2)
#
#                 d1name = '_'.join([str(elem) for elem in indices_1])
#                 d2name = '_'.join([str(elem) for elem in indices_2])
#                 comb_name_list = ['d1', d1name, 'deg1', d1str, 'n', 'd2', d2name, 'deg2', d2str]
#                 comb_name = '_'.join(comb_name_list)
#                 g16_file_base_name = str(file_counter) + '_' + ts_name + '_' + comb_name + '_g16'
#
#                 g16_file_path = save_folder + '/' + batch_foler + '/' + g16_file_base_name + '.gjf'
#                 with open(g16_file_path, 'wt') as f:
#                     f.write(script.format(name=g16_file_base_name, xyz=xyz_str))
#
#
# def detect_local_minima(arr):
#     # https://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
#     """
#     Takes an array and detects the troughs using the local maximum filter.
#     Returns a boolean mask of the troughs (i.e. 1 when
#     the pixel's value is the neighborhood maximum, 0 otherwise)
#     """
#     # define an connected neighborhood
#     # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#generate_binary_structure
#     neighborhood = morphology.generate_binary_structure(len(arr.shape),2)
#     # apply the local minimum filter; all locations of minimum value
#     # in their neighborhood are set to 1
#     # http://www.scipy.org/doc/api_docs/SciPy.ndimage.filters.html#minimum_filter
#     local_min = (filters.minimum_filter(arr, footprint=neighborhood)==arr)
#     # local_min is a mask that contains the peaks we are
#     # looking for, but also the background.
#     # In order to isolate the peaks we must remove the background from the mask.
#     #
#     # we create the mask of the background
#     background = (arr==0)
#     #
#     # a little technicality: we must erode the background in order to
#     # successfully subtract it from local_min, otherwise a line will
#     # appear along the background border (artifact of the local minimum filter)
#     # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#binary_erosion
#     eroded_background = morphology.binary_erosion(
#         background, structure=neighborhood, border_value=1)
#     #
#     # we obtain the final mask, containing only peaks,
#     # by removing the background from the local_min mask
#     detected_minima = local_min ^ eroded_background
#     return np.where(detected_minima)
#
#
# def run_2d_params(params, func):
#     a, b = params
#     return func(a, b)
#
#
# def optimizer(x0s):
#     res_list = list()
#     for x0 in x0s:
#         res = minimize(run_2d_params, x0=x0, args=g, method='Nelder-Mead', tol=1e-12)
#         res_list.append(res)
#     return res_list
#
#
# def highlight_max(data, color='yellow'):
#     """
#     highlight the maximum in a Series or DataFrame
#     """
#     attr = 'background-color: {}'.format(color)
#     if data.ndim == 1:  # Series from .apply(axis=0) or axis=1
#         is_max = data == data.max()
#         return [attr if v else '' for v in is_max]
#     else:  # from .apply(axis=None)
#         is_max = data == data.max().max()
#         return pd.DataFrame(np.where(is_max, attr, ''),
#                             index=data.index, columns=data.columns)
#
#
# def highlight_min(data, color='lightgreen'):
#     """
#     highlight the minimum in a Series or DataFrame
#     """
#     attr = 'background-color: {}'.format(color)
#     if data.ndim == 1:  # Series from .apply(axis=0) or axis=1
#         is_min = data == data.min()
#         return [attr if v else '' for v in is_min]
#     else:  # from .apply(axis=None)
#         is_min = data == data.min().min()
#         return pd.DataFrame(np.where(is_min, attr, ''),
#                             index=data.index, columns=data.columns)