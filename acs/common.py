"""
This module contains functions which are shared across multiple ACS modules.
"""


import os

import yaml
from typing import Any, List, Optional, Tuple, Union, Dict, Iterable
import shutil

from copy import deepcopy
from itertools import product
import numpy as np
import pandas as pd

import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
from scipy import interpolate
from scipy.optimize import minimize
from acs.backend import calculate_rmsd as calrmsd

# from arc.species.species import ARCSpecies
# from arc.species.converter import xyz_to_xyz_file_format, str_to_zmat, zmat_to_xyz, modify_coords
from arkane.ess import ess_factory, GaussianLog, MolproLog, OrcaLog, QChemLog, TeraChemLog
from acs.exceptions import ParserError, InputError, ConverterError
from arkane.exceptions import LogError
import rmgpy.constants as constants
from arkane.common import get_element_mass, mass_by_symbol, symbol_by_number
import qcelemental as qcel

from acs.converter.geom import xyz_dict_to_xyz_file


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

def gen_gaussian_optfreq_input_file(name: str,
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
        title_card = "#p opt=(calcall,ts,noeigentest,maxcycles=120) freq guess=mix scf=xqc iop(2/9=2000)"
    else:
        title_card = "#p opt=(calcall,noeigentest,maxcycles=120) freq guess=mix scf=xqc iop(2/9=2000)"

    script = f"""%chk={name}.chk
%mem={memory_mb}mb
%NProcShared={cpu_threads}

{title_card} {level_of_theory}

{comment}

{charge} {multiplicity}
{xyz_str}



"""
    return script

def gen_gaussian_irc_input_file(name: str,
                                    xyz_str: str,
                                    charge: int,
                                    multiplicity: int,
                                    memory_mb: int,
                                    cpu_threads: int,
                                    is_forward: bool,
                                    level_of_theory: str,
                                    comment: str = '',
                                    ) -> str:
    if is_forward:
        title_card = "#p irc=(calcfc,maxpoints=50,recalc=3,maxcycles=50, forward) scf=xqc iop(2/9=2000)"
    else:
        title_card = "#p irc=(calcfc,maxpoints=50,recalc=3,maxcycles=50, reverse) scf=xqc iop(2/9=2000)"

    script = f"""%chk={name}.chk
%mem={memory_mb}mb
%NProcShared={cpu_threads}

{title_card} {level_of_theory}

{comment}

{charge} {multiplicity}
{xyz_str}



"""
    return script


def gen_qchem_optfreq_input_file(name: str,
                                 xyz_str: str,
                                 charge: int,
                                 multiplicity: int,
                                 memory_mb: int,
                                 cpu_threads: int,
                                 is_ts: bool,
                                 level_of_theory: str,
                                 comment: str = '',
                                 ) -> str:
    # todo: add TS syntax for qchem
    # if is_ts:
    #     title_card = "#p opt=(calcall,ts,noeigentest,maxcycles=120) freq guess=mix scf=xqc iop(2/9=2000)"
    # else:
    #     title_card = "#p opt=(calcall,noeigentest,maxcycles=120) freq guess=mix scf=xqc iop(2/9=2000)"

    script = f"""
$molecule
{charge} {multiplicity}
{xyz_str}
$end

$rem
   JOBTYPE                   Opt
   METHOD                    {level_of_theory.split('/')[0]}
   BASIS                     {level_of_theory.split('/')[1]}
   SCF_ALGORITHM             DIIS_GDM
   MAX_DIIS_CYCLES           50
   THRESH_DIIS_SWITCH        11
   MAX_SCF_CYCLES            150
   SCF_CONVERGENCE           8
   THRESH                    11
   MEM_TOTAL                 {memory_mb}
$end

@@@

$molecule
   read
$end

$rem
   JOBTYPE                   FREQ
   METHOD                    {level_of_theory.split('/')[0]}
   BASIS                     {level_of_theory.split('/')[1]}
   SCF_ALGORITHM             DIIS_GDM
   MAX_DIIS_CYCLES           50
   THRESH_DIIS_SWITCH        11
   MAX_SCF_CYCLES            150
   SCF_CONVERGENCE           8
   THRESH                    11
   MEM_TOTAL                 {memory_mb}
$end


"""
    return script


def read_cosmo_gsolv(path: str,
                     use_hartree: bool = True,
                     ) -> float:
    """
    Read Gsolv from cosmo tab output. If the tab file contains results from multiple temperatures, only the first one
    is read.

    Args:
        path: cosmo tab file path.
        use_hartree: whether to return Gsolv in hartree (True, Default) or kcal/mol.

    Returns:
        Float of Gsolv in hartree.
    """

    with open(path, 'r') as f:
        lines = f.readlines()

    spc_name = os.path.splitext(os.path.basename(path))[0]
    hartree_to_kcal_per_mol = 627.5094740631

    for l in tuple(lines):
        if spc_name in l:
            g_solv = float(l.split()[5])
            break
    else:
        raise ParserError(f'Did not found Gsolv for species {spc_name} in file {path}')

    if use_hartree:
        return g_solv/hartree_to_kcal_per_mol
    else:
        return g_solv


def gen_gaussian_cosmo_sp_input_file(name: str,
                                    xyz_str: str,
                                    charge: int,
                                    multiplicity: int,
                                    memory_mb: int,
                                    cpu_threads: int,
                                    comment: str = '',
                                    ) -> str:
    """
    Generate gaussian BP-TZVP cosmo single point input file.
    Notice: Gaussian has strict rules on the position of blank lines in the input file. Be aware and do test if you
            are changing the script variable below.
    Args:
        name:
        xyz_str:
        charge:
        multiplicity:
        memory_mb:
        cpu_threads:
        comment:

    Returns:

    """

    title_card = "#P BVP86/TZVP/DGA1 scf=(tight,novaracc) SCRF=COSMORS NoSymm"

    script = f"""%chk={name}.chk
%mem={memory_mb}mb
%NProcShared={cpu_threads}

{title_card}

{comment}

{charge} {multiplicity}
{xyz_str}
{name}.cosmo"""
    return script


def gen_cosmo_rs_input_file(name: str,
                            ) -> str:
    """
    Generate cosmo rs input file for delta G solvation calculation. Currently implemented specifically for BP-TZVP
    calculation in ternary mixture with 70 mol % water, 30 mol % methanol as solvents, and infinite dilution solute.

    Args:
        name: gaussian cosmo output.

    Returns:

    """

    script = f"""ctd = BP_TZVPD_FINE_20.ctd cdir = "/home/gridsan/groups/RMG/Software/COSMOtherm2020/COSMOtherm/CTDATA-FILES" ldir = "/home/gridsan/groups/RMG/Software/COSMOtherm2020/licensefiles"
unit notempty wtln ehfile
!! generated by ACS !!       
f = "h2o_c0.cosmo" fdir="/home/gridsan/oscarwu/RMG_group_folder/COSMO_database/COSMObase2020/BP-TZVPD-FINE/h" VPfile 
f = "methanol_c0.cosmo" fdir="/home/gridsan/oscarwu/RMG_group_folder/COSMO_database/COSMObase2020/BP-TZVPD-FINE/m" VPfile 
f = "{name}.cosmo" fdir="." VPfile 
henry  xh={{ 70 30 0 }}  tc=25.0 GSOLV # Automatic Henry Law coefficient Calculation 
    """
    return script

def gen_orca_dlpno_sp_input_file(xyz_str: str,
                                charge: int,
                                multiplicity: int,
                                memory_mb: int,
                                cpu_threads: int,
                                comment: str = '',
                                ) -> str:
    """
    # todo: consider other methods.
    Generate ORCA DLPNO single point input file.
    Args:
        xyz_str:
        charge:
        multiplicity:
        memory_mb:
        cpu_threads:
        comment:

    Returns:

    """

    if multiplicity == 1:
        reference = 'rHF'
    elif multiplicity == 2:
        reference = 'uHF'
    else:
        raise NotImplementedError

    script = f"""!{reference} dlpno-ccsd(t) def2-tzvp def2-tzvp/c NormalPNO
!NRSCF # using Newton Raphson SCF algorithm 
!sp 

%maxcore {memory_mb}
%pal
nprocs {cpu_threads}
end
%scf # recommended SCF settings 
NRMaxIt 400
NRStart 0.00005
MaxIter 500
end

{comment}

* xyz {charge} {multiplicity}
{xyz_str}*"""
    return script

def process_opt_freq_output(logfile, ess_software, is_ts=True, check_neg_freq=True):
    if ess_software == 'gaussian':
        if not check_gaussian_normal_termination(logfile):
            raise ParserError('Gaussian error termination.')
    elif ess_software == 'qchem':
        if not check_qchem_normal_termination(logfile):
            raise ParserError('QChem error termination.')
    info = dict()
    info['freq'] = get_freq(logfile, ess_software, check_neg_freq=check_neg_freq, ts=is_ts)
    info['xyz_dict'] = get_geometry(logfile)
    info['electronic_energy'] = get_e_elect_from_log(logfile)
    info['unscaled_zpe'] = get_unscaled_zpe(logfile)

    return info


def get_freq(logfile, ess_software, check_neg_freq=True, min_neg_freq=-2400, max_neg_freq=-200, ts=True):
    freq = parse_frequencies(logfile, software=ess_software)
    neg_freq = tuple([float(x) for x in freq if x < 0])
    if check_neg_freq:
        if ts:
            if len(neg_freq) == 0:
                raise ParserError('Did not find any negative frequencies.')
            elif len(neg_freq) > 1:
                raise ParserError(f'Find more than one negative frequencies: {neg_freq}')
            elif neg_freq[0] < min_neg_freq:
                raise ParserError(f'Value of negative frequency {neg_freq} lower than minimum {min_neg_freq}')
            elif neg_freq[0] > max_neg_freq:
                raise ParserError(f'Value of negative frequency {neg_freq} higher than maximum {max_neg_freq}')
        else:
            if len(neg_freq):
                raise ParserError(f'Find negative frequencies for non-TS species: {neg_freq}')
    return freq, neg_freq


def parse_frequencies(path: str,
                      software: str,
                      ) -> np.ndarray:
    """
    Parse the frequencies from a freq job output file.

    Args:
        path (str): The log file path.
        software (str): The ESS.

    Returns:
        np.ndarray: The parsed frequencies (in cm^-1).
    """
    lines = _get_lines_from_file(path)
    freqs = np.array([], np.float64)
    if software.lower() == 'qchem':
        for line in lines:
            if ' Frequency:' in line:
                items = line.split()
                for i, item in enumerate(items):
                    if i:
                        freqs = np.append(freqs, [(float(item))])
    elif software.lower() == 'gaussian':
        with open(path, 'r') as f:
            line = f.readline()
            while line != '':
                if 'and normal coordinates' in line:
                    freqs = np.array([], np.float64)
                if 'Frequencies --' in line:
                    freqs = np.append(freqs, [float(frq) for frq in line.split()[2:]])
                line = f.readline()
    elif software.lower() == 'molpro':
        read = False
        for line in lines:
            if 'Nr' in line and '[1/cm]' in line:
                continue
            if read:
                if line == os.linesep:
                    read = False
                    continue
                freqs = np.append(freqs, [float(line.split()[-1])])
            if 'Low' not in line and 'Vibration' in line and 'Wavenumber' in line:
                read = True
    elif software.lower() == 'orca':
        with open(path, 'r') as f:
            line = f.readline()
            read = True
            while line:
                if 'VIBRATIONAL FREQUENCIES' in line:
                    while read:
                        if not line.strip():
                            line = f.readline()
                        elif not line.split()[0] == '0:':
                            line = f.readline()
                        else:
                            read = False
                    while line.strip():
                        if float(line.split()[1]) != 0.0:
                            freqs = np.append(freqs, [float(line.split()[1])])
                        line = f.readline()
                    break
                else:
                    line = f.readline()
    elif software.lower() == 'terachem':
        read_output = False
        for line in lines:
            if '=== Mode' in line:
                # example: '=== Mode 1: 1198.526 cm^-1 ==='
                freqs = np.append(freqs, [float(line.split()[3])])
            elif 'Vibrational Frequencies/Thermochemical Analysis After Removing Rotation and Translation' in line:
                read_output = True
                continue
            elif read_output:
                if 'Temperature (Kelvin):' in line or 'Frequency(cm-1)' in line:
                    continue
                if not line.strip():
                    break
                # example:
                # 'Mode  Eigenvalue(AU)  Frequency(cm-1)  Intensity(km/mol)   Vib.Temp(K)      ZPE(AU) ...'
                # '  1     0.0331810528   170.5666870932      52.2294230772  245.3982965841   0.0003885795 ...'
                freqs = np.append(freqs, [float(line.split()[2])])

    else:
        raise ParserError(f'parse_frequencies() can currently only parse Gaussian, Molpro, Orca, QChem and TeraChem '
                          f'files, got {software}')
    return freqs


def check_gaussian_normal_termination(logfile):
    with open(logfile, 'r') as f:
        lines = f.readlines()
        forward_lines = tuple(lines)
    for line in forward_lines:
        if 'Error termination' in line:
            return False
    else:
        return True


def check_qchem_normal_termination(logfile):
    with open(logfile, 'r') as f:
        lines = f.readlines()
        forward_lines = tuple(lines)
    for line in forward_lines:
        # todo: how does qchem error?
        if 'Error termination' in line:
            return False
    else:
        return True


def parse_zpe(path: str) -> Optional[float]:
    """
    Determine the calculated ZPE from a frequency output file

    Args:
        path (str): The path to a frequency calculation output file.

    Returns:
        Optional[float]: The calculated zero point energy in kJ/mol.
    """
    if not os.path.isfile(path):
        raise InputError('Could not find file {0}'.format(path))
    log = ess_factory(fullpath=path)
    try:
        zpe = log.load_zero_point_energy() * 0.001  # convert to kJ/mol
    except (LogError, NotImplementedError):
        zpe = None
    return zpe


def get_unscaled_zpe(logfile):
    energy_dict = dict()
    try:
        zpe_kj_mol = parse_zpe(logfile)
    except:
        raise ParserError('Cannot parse energy from log file.')

    zpe_j_mol = zpe_kj_mol * 1000
    zpe_kcal_mol = zpe_j_mol / 4184
    zpe_unscaled = round(zpe_j_mol/(constants.E_h * constants.Na), 9)

    energy_dict['J/mol'] = zpe_j_mol
    energy_dict['kJ/mol'] = zpe_kj_mol
    energy_dict['kcal/mol'] = zpe_kcal_mol
    energy_dict['hartree'] = zpe_unscaled
    return energy_dict


def get_e_elect_from_log(logfile):
    energy_dict = dict()
    try:
        e_kj_mol = parse_e_elect(logfile)
    except:
        raise ParserError('Cannot parse energy from the log file.')

    e_j_mol = e_kj_mol * 1000
    e_kcal_mol = e_j_mol / 4184
    e_elect = round(e_j_mol/(constants.E_h * constants.Na), 9)

    energy_dict['J/mol'] = e_j_mol
    energy_dict['kJ/mol'] = e_kj_mol
    energy_dict['kcal/mol'] = e_kcal_mol
    energy_dict['hartree'] = e_elect
    return energy_dict


def parse_e_elect(path: str,
                  zpe_scale_factor: float = 1.,
                  ) -> Optional[float]:
    """
    Parse the electronic energy from an sp job output file.

    Args:
        path (str): The ESS log file to parse from.
        zpe_scale_factor (float): The ZPE scaling factor, used only for composite methods in Gaussian via Arkane.

    Returns:
        Optional[float]: The electronic energy in kJ/mol.
    """
    if not os.path.isfile(path):
        raise InputError(f'Could not find file {path}')
    log = ess_factory(fullpath=path)
    try:
        e_elect = log.load_energy(zpe_scale_factor) * 0.001  # convert to kJ/mol
    except (LogError, NotImplementedError):
        e_elect = None
    return e_elect


def get_geometry(logfile):
    xyz = parse_geometry(logfile)
    return xyz


def parse_geometry(path: str) -> Optional[Dict[str, tuple]]:
    """
    Parse the xyz geometry from an ESS log file.

    Args:
        path (str): The ESS log file to parse from.

    Returns:
        Optional[Dict[str, tuple]]: The cartesian geometry.
    """
    log = ess_factory(fullpath=path)
    try:
        coords, number, _ = log.load_geometry()
    except LogError:
        raise ParserError(f'Could not parse xyz from {path}')
    return xyz_from_data(coords=coords, numbers=number)

def xyz_from_data(coords, numbers=None, symbols=None, isotopes=None):
    """
    Get the ARC xyz dictionary format from raw data.
    Either ``numbers`` or ``symbols`` must be specified.
    If ``isotopes`` isn't specified, the most common isotopes will be assumed for all elements.

    Args:
        coords (tuple, list): The xyz coordinates.
        numbers (tuple, list, optional): Element nuclear charge numbers.
        symbols (tuple, list, optional): Element symbols.
        isotopes (tuple, list, optional): Element isotope numbers.

    Returns:
        dict: The ARC dictionary xyz format.

    Raises:
        ConverterError: If neither ``numbers`` nor ``symbols`` are specified, if both are specified,
                        or if the input lengths aren't consistent.
    """
    if isinstance(coords, np.ndarray):
        coords = tuple(tuple(coord.tolist()) for coord in coords)
    elif isinstance(coords, list):
        coords = tuple(tuple(coord) for coord in coords)
    if numbers is not None and isinstance(numbers, np.ndarray):
        numbers = tuple(numbers.tolist())
    elif numbers is not None and isinstance(numbers, list):
        numbers = tuple(numbers)
    if symbols is not None and isinstance(symbols, list):
        symbols = tuple(symbols)
    if isotopes is not None and isinstance(isotopes, list):
        isotopes = tuple(isotopes)
    if not isinstance(coords, tuple):
        raise ConverterError('Expected coords to be a tuple, got {0} which is a {1}'.format(coords, type(coords)))
    if numbers is not None and not isinstance(numbers, tuple):
        raise ConverterError('Expected numbers to be a tuple, got {0} which is a {1}'.format(numbers, type(numbers)))
    if symbols is not None and not isinstance(symbols, tuple):
        raise ConverterError('Expected symbols to be a tuple, got {0} which is a {1}'.format(symbols, type(symbols)))
    if isotopes is not None and not isinstance(isotopes, tuple):
        raise ConverterError('Expected isotopes to be a tuple, got {0} which is a {1}'.format(isotopes, type(isotopes)))
    if numbers is None and symbols is None:
        raise ConverterError('Must set either "numbers" or "symbols". Got neither.')
    if numbers is not None and symbols is not None:
        raise ConverterError('Must set either "numbers" or "symbols". Got both.')
    if numbers is not None:
        symbols = tuple(symbol_by_number[number] for number in numbers)
    if len(coords) != len(symbols):
        raise ConverterError(f'The length of the coordinates ({len(coords)}) is different than the length of the '
                             f'numbers/symbols ({len(symbols)}).')
    if isotopes is not None and len(coords) != len(isotopes):
        raise ConverterError(f'The length of the coordinates ({len(coords)}) is different than the length of isotopes '
                             f'({len(isotopes)}).')
    if isotopes is None:
        isotopes = tuple(get_most_common_isotope_for_element(symbol) for symbol in symbols)
    xyz_dict = {'symbols': symbols, 'isotopes': isotopes, 'coords': coords}
    return xyz_dict

def get_most_common_isotope_for_element(element_symbol):
    """
    Get the most common isotope for a given element symbol.

    Args:
        element_symbol (str): The element symbol.

    Returns:
        int: The most common isotope number for the element.
             Returns ``None`` for dummy atoms ('X').
    """
    if element_symbol == 'X':
        # this is a dummy atom (such as in a zmat)
        return None
    mass_list = mass_by_symbol[element_symbol]
    if len(mass_list[0]) == 2:
        # isotope contribution is unavailable, just get the first entry
        isotope = mass_list[0][0]
    else:
        # isotope contribution is unavailable, get the most common isotope
        isotope, isotope_contribution = mass_list[0][0], mass_list[0][2]
        for iso in mass_list:
            if iso[2] > isotope_contribution:
                isotope_contribution = iso[2]
                isotope = iso[0]
    return isotope

def _get_lines_from_file(path: str) -> List[str]:
    """
    A helper function for getting a list of lines from a file.

    Args:
        path (str): The file path.

    Raises:
        InputError: If the file could not be read.

    Returns:
        List[str]: Entries are lines from the file.
    """
    if os.path.isfile(path):
        with open(path, 'r') as f:
            lines = f.readlines()
    else:
        raise InputError(f'Could not find file {path}')
    return lines


def cluster_confs_by_rmsd_with_id(labeled_xyzs: Iterable[Tuple[str, Dict[str, tuple]]],
                                  rmsd_threshold: float = 1e-2,
                                  ) -> Tuple[Tuple[str, Dict[str, tuple]]]:
    """
    Cluster conformers with the same atom orders using RMSD of distance matrices. Work for both TS and non-TS conformers.

    Intended for finding structurally distinct conformers from a pool of conformers.
    Suitable scenarios:
        1. filter a pool of conformers with their geometry optimized at some level.
    Not suitable for:
        1. cluster conformers (not optimized) that are sampling of a well or a saddle point (these conformers may have
           large difference in RMSE, but they really should be representing the same well or saddle point).

    Args:
        labeled_xyzs (Iterable): Conformers with the same atom orders with ids.
        rmsd_threshold (float): The minimum RMSD to consider two conformers as distinct
                                (i.e., if rmsd > rmsd_threshold, then two conformers are considered distinctive).

    Returns:
        Tuple[Dict[str, tuple]]: Conformers with distinctive geometries.
    """
    labeled_xyzs = tuple(labeled_xyzs)
    fingerprint = labeled_xyzs[0][0]
    xyz = labeled_xyzs[0][1]
    distinct_id_xyz_pair = [(fingerprint, xyz)]
    distinct_xyzs = [xyz]

    for id_xyz_pair in labeled_xyzs:
        fingerprint = id_xyz_pair[0]
        xyz = id_xyz_pair[1]
        rmsd_list = [compare_confs(xyz, distinct_xyz, rmsd_score=True) for distinct_xyz in tuple(distinct_xyzs)]
        if all([rmsd > rmsd_threshold for rmsd in tuple(rmsd_list)]):
            distinct_xyzs.append(xyz)
            distinct_id_xyz_pair.append((fingerprint, xyz))
    return tuple(distinct_id_xyz_pair)

def compare_confs(xyz1: dict,
                  xyz2: dict,
                  rtol: float = 1e-5,
                  atol: float = 1e-5,
                  rmsd_score: bool = False,
                  ) -> Union[float, bool]:
    """
    Compare two Cartesian coordinates representing conformers using distance matrices.

    The relative difference (``rtol`` * abs(value in xyz2)) and the absolute difference ``atol``
    are added together to compare against the absolute difference between (value in xyz1) and (value in xyz2).

    Args:
        xyz1 (dict): Conformer 1.
        xyz2 (dict): Conformer 2.
        rtol (float): The relative tolerance parameter (see Notes).
        atol (float): The absolute tolerance parameter (see Notes).
        rmsd_score (bool): Whether to output a root-mean-square deviation score of the two distance matrices.

    Returns:
        Union[float, bool]:
            - If ``rmsd_score`` is ``False`` (default): Whether the two conformers have almost equal atom distances.
              ``True`` if they do.
            - If ``rmsd_score`` is ``True``: The RMSD score of two distance matrices.
    """
    xyz1, xyz2 = check_xyz_dict(xyz1), check_xyz_dict(xyz2)
    dmat1, dmat2 = xyz_to_dmat(xyz1), xyz_to_dmat(xyz2)
    if rmsd_score:
        # method 1: use distance matrix
        # distance matrix is symmetric, only need the upper triangular part to compute rmsd
        rmsd_1 = calc_rmsd(np.triu(dmat1), np.triu(dmat2))

        # method 2: use Kabsch algorithm (works better for indistinguishable H atoms)
        # https://github.com/charnley/rmsd
        rmsd_2 = calc_rmsd_wrapper(xyz1, xyz2)

        rmsd = min((rmsd_1, rmsd_2))

        return rmsd
    else:
        return almost_equal_lists(dmat1, dmat2, rtol=rtol, atol=atol)

def almost_equal_lists(iter1: list or tuple or np.ndarray,
                       iter2: list or tuple or np.ndarray,
                       rtol: float = 1e-05,
                       atol: float = 1e-08,
                       ) -> bool:
    """
    A helper function for checking whether two iterables are almost equal.

    Args:
        iter1 (list, tuple, np.array): An iterable.
        iter2 (list, tuple, np.array): An iterable.
        rtol (float, optional): The relative tolerance parameter.
        atol (float, optional): The absolute tolerance parameter.

    Returns:
        bool: ``True`` if they are almost equal, ``False`` otherwise.
    """
    if len(iter1) != len(iter2):
        return False
    for entry1, entry2 in zip(iter1, iter2):
        if isinstance(entry1, (list, tuple, np.ndarray)) and isinstance(entry2, (list, tuple, np.ndarray)):
            return almost_equal_lists(iter1=entry1, iter2=entry2, rtol=rtol, atol=atol)
        else:
            if isinstance(entry1, (int, float)) and isinstance(entry2, (int, float)):
                if not np.isclose([entry1], [entry2], rtol=rtol, atol=atol):
                    return False
            else:
                if entry1 != entry2:
                    return False
    return True

def calc_rmsd(x: np.array,
              y: np.array,
              ) -> float:
    """
    Compute the root-mean-square deviation between two matrices.

    Args:
        x (np.array): Matrix 1.
        y (np.array): Matrix 2.

    Returns:
        float: The RMSD score of two matrices.
    """
    d = x - y
    n = x.shape[0]
    sqr_sum = (d**2).sum()
    rmsd = np.sqrt(sqr_sum/n)
    return float(rmsd)

def calc_rmsd_wrapper(xyz_1: dict,
                      xyz_2: dict,
                      ) -> float:
    """
    A wrapper for https://github.com/charnley/rmsd.
    Calculate RMSD from xyz dict.

    Args:
        xyz_1 (dict): XYZ coordinate of species 1.
        xyz_2 (dict): XYZ coordinate of species 2.

    Returns:
        float: The RMSD between two species.
    """

    p_all_atoms, p_all = np.array(xyz_1['symbols']), np.array(xyz_1['coords'])
    q_all_atoms, q_all = np.array(xyz_2['symbols']), np.array(xyz_2['coords'])

    p_size = p_all.shape[0]
    q_size = q_all.shape[0]

    if not p_size == q_size:
        raise ParserError('Does not make sense to compute RMSD of conformers with different sizes.')

    p_coord = deepcopy(p_all)
    q_coord = deepcopy(q_all)
    p_atoms = deepcopy(p_all_atoms)
    q_atoms = deepcopy(q_all_atoms)

    p_cent = calrmsd.centroid(p_coord)
    q_cent = calrmsd.centroid(q_coord)
    p_coord -= p_cent
    q_coord -= q_cent

    rmsds = list()
    rotation_methods = (calrmsd.kabsch_rmsd, calrmsd.quaternion_rmsd)
    reorder_methods = (calrmsd.reorder_hungarian, calrmsd.reorder_inertia_hungarian)

    for rotation_method, reorder_method in product(rotation_methods, reorder_methods):
        try:
            result_rmsd, q_swap, q_reflection, q_review = calrmsd.check_reflections(
                p_atoms,
                q_atoms,
                p_coord,
                q_coord,
                reorder_method=reorder_method,
                rotation_method=rotation_method,
                keep_stereo=True)
        except ValueError:
            continue
        rmsds.append(float(result_rmsd))

    rmsd = min(rmsds)
    return rmsd


def check_xyz_dict(xyz):
    """
    Check that the xyz dictionary entered is valid.
    If it is a string, convert it.
    If it is a Z matrix, convert it to cartesian coordinates,
    If isotopes are not in xyz_dict, common values will be added.

    Args:
        xyz (dict, str): The xyz dictionary.

    Raises:
        ConverterError: If ``xyz`` is of wrong type or is missing symbols or coords.
    """
    xyz_dict = xyz
    if not isinstance(xyz_dict, dict):
        raise ConverterError(f'Expected a dictionary, got {type(xyz_dict)}')
    if 'vars' in list(xyz_dict.keys()):
        raise ConverterError('Expected xyz got zmat')
    if 'symbols' not in list(xyz_dict.keys()):
        raise ConverterError(f'XYZ dictionary is missing symbols. Got:\n{xyz_dict}')
    if 'coords' not in list(xyz_dict.keys()):
        raise ConverterError(f'XYZ dictionary is missing coords. Got:\n{xyz_dict}')
    if len(xyz_dict['symbols']) != len(xyz_dict['coords']):
        raise ConverterError(f'Got {len(xyz_dict["symbols"])} symbols and {len(xyz_dict["coords"])} '
                             f'coordinates:\n{xyz_dict}')
    if 'isotopes' not in list(xyz_dict.keys()):
        xyz_dict = xyz_from_data(coords=xyz_dict['coords'], symbols=xyz_dict['symbols'])
    if len(xyz_dict['symbols']) != len(xyz_dict['isotopes']):
        raise ConverterError(f'Got {len(xyz_dict["symbols"])} symbols and {len(xyz_dict["isotopes"])} '
                             f'isotopes:\n{xyz_dict}')
    return xyz_dict


def xyz_to_dmat(xyz_dict: dict) -> np.array:
    """
    Convert Cartesian coordinates to a distance matrix.

    Args:
        xyz_dict (dict): The Cartesian coordinates,

    Returns:
        list: the distance matrix.
    """
    xyz_dict = check_xyz_dict(xyz_dict)
    dmat = qcel.util.misc.distance_matrix(a=np.array(xyz_to_coords_list(xyz_dict)),
                                          b=np.array(xyz_to_coords_list(xyz_dict)))
    return dmat


def xyz_to_coords_list(xyz_dict):
    """
    Get the coords part of an xyz dict as a (mutable) list of lists (rather than a tuple of tuples)

    Args:
        xyz_dict (dict): The ARC xyz format.

    Returns:
        list: The coordinates.
    """
    xyz_dict = check_xyz_dict(xyz_dict)
    coords_tuple = xyz_dict['coords']
    coords_list = list()
    for coords_tup in coords_tuple:
        coords_list.append([coords_tup[0], coords_tup[1], coords_tup[2]])
    return coords_list

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