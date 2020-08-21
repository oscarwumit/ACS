"""
This module computes MSTST rate coefficient.
"""

import logging
import shutil
import os
import yaml
import argparse
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from arkane.input import process_model_chemistry
from arkane.common import symbol_by_number, get_principal_moments_of_inertia
from arkane.encorr.corr import assign_frequency_scale_factor, get_atom_correction
from arkane.ess import ess_factory
from arkane.modelchem import LOT, LevelOfTheory, CompositeLevelOfTheory, model_chem_to_lot

from rmgpy import constants
from rmgpy.kinetics.tunneling import Eckart, Wigner
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule.element import get_element
from rmgpy.molecule.molecule import Molecule
from rmgpy.qm.qmdata import QMData
from rmgpy.qm.symmetry import PointGroupCalculator
from rmgpy.species import Species, TransitionState
from rmgpy.statmech import (IdealGasTranslation,
                            NonlinearRotor,
                            LinearRotor,
                            HarmonicOscillator,
                            Conformer)

from acs.common import read_yaml_file


def get_symmetry(coords, atom_numbers, scr_dir=None):

    scr_dir = scr_dir or os.path.join('.', 'scratch')
    os.makedirs(scr_dir, exist_ok=True)
    
    symmetry = optical_isomers = 1
    try:
        qmdata = QMData(
            groundStateDegeneracy=1,  # Only needed to check if valid QMData
            numberOfAtoms=len(atom_numbers),
            atomicNumbers=atom_numbers,
            atomCoords=(coords, 'angstrom'),
            energy=(0.0, 'kcal/mol')  # Only needed to avoid error
        )
        settings = type('', (),
                        dict(symmetryPath='symmetry',
                             scratchDirectory=scr_dir))()
        pgc = PointGroupCalculator(settings, '0', qmdata)  # '0' is an unique id used for calculator
        pg = pgc.calculate()
        if pg is not None:
            symmetry = pg.symmetry_number
            optical_isomers = 2 if pg.chiral else 1
            logging.debug(f"Symmetry algorithm found {optical_isomers} optical isomers "
                          f"and a symmetry number of {symmetry}")
        else:
            logging.warning("Symmetry algorithm errored when computing point group. "
                            "Using symmetry number=1 and optical isomers = 1 for "
                            "further calculations, which may not be true.")
        return symmetry, optical_isomers
    finally:
        shutil.rmtree(scr_dir)

def get_lot_and_freq_scale(energy_level: str,
                           freq_level: str,
                           energy_software: str = 'orca',
                           freq_scale: Optional[float] = None):
    # Get energy level and assign software
    energy_level = LevelOfTheory(energy_level)
    energy_software = energy_software or energy_log.get_software()
    energy_level = energy_level.update(software=energy_software)

    # Get freq level
    freq_level = LevelOfTheory(freq_level)

    # Assign level of theory and frequency scale factor
    if energy_level.to_model_chem() in ['cbsqb3'] \
           or energy_level == freq_level:
        level_of_theory = energy_level
    else:
        level_of_theory = CompositeLevelOfTheory(freq=freq_level, energy=energy_level)
    if freq_scale is None:
        try:
            freq_scale = assign_frequency_scale_factor(level_of_theory)
        except:
            logging.warning('Setting freq scale to 1.0')
            freq_scale = 1.0
    logging.warning(f'freq scale is {freq_scale}')
    return level_of_theory, freq_scale

def get_rotational_mode(coords, number, external_symmetry=None):
    if not external_symmetry:
        external_symmetry, _ = get_symmetry(coords, number)
    
    # Rotational
    moments_of_inertia = get_principal_moments_of_inertia(coords=coords,
                                                          numbers=number,)[0]
    if any([moment_of_inertia == 0.0 for moment_of_inertia in moments_of_inertia]):
        # this is a linear rotor
        moments_of_inertia = [moment_of_inertia for moment_of_inertia in moments_of_inertia
                              if moment_of_inertia != 0.0]
        if abs(moments_of_inertia[0] - moments_of_inertia[1]) > 0.01:
            raise Exceptions(f'Expected two identical moments of inertia for a linear rigis rotor, '
                             f'but got {moments_of_inertia}')
        return LinearRotor(inertia=(moments_of_inertia[0], "amu*angstrom^2"),
                           symmetry=external_symmetry)
    else:
        # this is a non-linear rotor
        return NonlinearRotor(inertia=(moments_of_inertia, "amu*angstrom^2"),
                              symmetry=external_symmetry)
    
def get_element_counts(number):
    # Get atoms count
    atoms = {}
    for atom_num in number:
        try:
            symbol = symbol_by_number[atom_num]
        except KeyError:
            raise ElementError('Could not recognize element number {0}.'.format(atom_num))
        atoms[symbol] = atoms.get(symbol, 0) + 1
    return atoms

def get_rmg_conformer_from_logs(label,
                                 energy_log_path,
                                 freq_log_path,
                                 level_of_theory,
                                 freq_scale=1,
                                 multiplicity=1,
                                 molecule=None,
                                 use_atom_corrections=True,
                                 use_bond_corrections=False,):
    
    energy_log = ess_factory(energy_log_path)
    freq_log = ess_factory(freq_log_path)
    
    # Get the coords, atom_numbers, and mass
    coords, number, mass = freq_log.load_geometry()
    # Get the symmetry info of the molecule
    external_symmetry, optical_isomers = get_symmetry(coords, number)
    
    # This script will automatically assign modes but with wrong E0 and unscaled freqs
    conformer, unscaled_frequencies = freq_log.load_conformer(symmetry=external_symmetry,
                                                              spin_multiplicity=multiplicity,
                                                              optical_isomers=optical_isomers,
                                                              label=label)
    
    conformer.coordinates = (coords, "angstroms")
    conformer.number = number
    conformer.mass = (mass, "amu")
    
    zpe_scale_factor = freq_scale / 1.014
    e_electronic = energy_log.load_energy(zpe_scale_factor)

    # Atom energy corrections
    if use_atom_corrections:
        # Get atoms count
        atoms = get_element_counts(number)
        atom_corrections = get_atom_correction(level_of_theory, atoms)
    else:
        atom_corrections = 0

    # Bond energy corrections
    if use_bond_corrections:
        # Get bonds count
        try:
            bonds = molecule.enumerate_bonds()
            bond_corrections = get_bac(level_of_theory, bonds, coords, number,
                                       bac_type='p', multiplicity=multiplicity)
        except AttributeError:
            raise ValueError('Cannot get BAC, since argument ``molecule`` is not provided.')
    else:
        bond_corrections = 0

    e_electronic_with_corrections = e_electronic + atom_corrections + bond_corrections
    zpe = freq_log.load_zero_point_energy() * zpe_scale_factor if len(number) > 1 else 0
    conformer.E0 = ((e_electronic_with_corrections + zpe) * 0.001, 'kJ/mol')

    # Correct the frequencies
    for mode in conformer.modes:
        if isinstance(mode, HarmonicOscillator):
            mode.frequencies = (np.array(unscaled_frequencies) * freq_scale, "cm^-1")
    
    return conformer

def get_rmg_conformer(label,
                      level_of_theory,
                      e_electronic,
                      frequencies,
                      coords,
                      numbers,
                      mass,
                      multiplicity=1,
                      freq_scale=1,
                      molecule=None,
                      use_atom_corrections=True,
                      use_bond_corrections=False,):
    
    external_symmetry, optical_isomers = get_symmetry(coords, numbers,)
    
    modes = []
    # Translational
    translation = IdealGasTranslation(mass=mass)
    modes.append(translation)
    
    # Rotational
    rotation = get_rotational_mode(coords, numbers,
                                   external_symmetry=external_symmetry)
    modes.append(rotation)
    
    # Vibrational
    frequencies = np.array(frequencies)
    frequencies = frequencies[frequencies >= 0]
    vibration = HarmonicOscillator(frequencies=(frequencies * freq_scale, "cm^-1"))
    modes.append(vibration)
    
    # Atom energy corrections
    if use_atom_corrections:
        atoms = get_element_counts(numbers)
        atom_corrections = get_atom_correction(level_of_theory, atoms)
    else:
        atom_corrections = 0

    # Bond energy corrections
    if use_bond_corrections:
        # Get bonds count
        try:
            bonds = molecule.enumerate_bonds()
            bond_corrections = get_bac(level_of_theory, bonds, coords, number,
                                       bac_type='p', multiplicity=multiplicity)
        except AttribureError:
            raise ValueError('Cannot get BAC, since argument ``molecule`` is not provided.')
    else:
        bond_corrections = 0
        
    e_electronic_with_corrections = e_electronic + atom_corrections + bond_corrections
    
    if len(numbers) > 1:
        zpe_scale_factor = freq_scale / 1.014
        scaled_zpe = 0.5 * constants.h * constants.c * constants.Na \
                     * np.sum(frequencies) * 100 * zpe_scale_factor
    else:
        scaled_zpe = 0

    e0 = ((e_electronic_with_corrections + scaled_zpe) * 0.001, 'kJ/mol')
    
    return Conformer(E0=e0, modes=modes,
                     spin_multiplicity=multiplicity,
                     optical_isomers=optical_isomers)


def get_rmg_conformer_from_yml(path, freq_scale=None):
    # No bond correction at this moment
    
    results = read_yaml_file(path)
    energy_level = results['level_of_theory']['sp_after_opt']
    freq_level = results['level_of_theory']['fine_opt_freq']

    level_of_theory, freq_scale = get_lot_and_freq_scale(energy_level=energy_level,
                                                         freq_level=freq_level,
                                                         freq_scale=freq_scale)

    multiplicity = results['species']['multiplicity']

    rmg_conformers = []
    for hash_id in results['valid_conformer_hash_ids']:
        cur_conformer = results['conformers'][hash_id]
        xyz = cur_conformer['arc_xyz_after_fine_opt']
        coords = np.array(xyz['coords'])
        atom_numbers = [get_element(symbol).number for symbol in xyz['symbols']]
        mass = (sum([get_element(symbol).mass for symbol in xyz['symbols']])/constants.Na, 'kg')

        e_electronic = cur_conformer['energy']['sp_include_solv_correction'] * 2625500  # hartree to J/mol 

        rmg_conformer = get_rmg_conformer(label=hash_id,
                                          level_of_theory=level_of_theory,
                                          e_electronic=e_electronic,
                                          frequencies=cur_conformer['frequencies'],
                                          coords=coords,
                                          numbers=atom_numbers,
                                          mass=mass,
                                          multiplicity=multiplicity,
                                          freq_scale=freq_scale,
                                          molecule=None,
                                          use_atom_corrections=True,
                                          use_bond_corrections=False)
        rmg_conformers.append(rmg_conformer)
    return rmg_conformers

def get_E0min_and_Q(Ts, conformers, option='ms'):
    if option not in ['min', 'max', 'ms']:
        raise
        
    if option == 'min':
        E0_min = min([conf.E0.value_si for conf in conformers])
        min_idx = np.argmin([conf.E0.value_si for conf in conformers])
        E0_ref = E0_min
        Q = np.zeros_like(Ts)
        for idx in range(Ts.shape[0]):
            T = Ts[idx]
            Q[idx] = np.array([conformers[min_idx].get_partition_function(T)])
        
    elif option == 'max': 
        E0_max = max([conf.E0.value_si for conf in conformers])
        max_idx = np.argmax([conf.E0.value_si for conf in conformers])
        E0_ref = E0_max
        Q = np.zeros_like(Ts)
        for idx in range(Ts.shape[0]):
            T = Ts[idx]
            Q[idx] = np.array([conformers[max_idx].get_partition_function(T)])
        
    elif option == 'ms':
        E0_min = min([conf.E0.value_si for conf in conformers])
        E0_ref = E0_min
        Uis = np.array([conf.E0.value_si - E0_ref for conf in conformers])
        Q = np.zeros_like(Ts)
        for idx in range(Ts.shape[0]):
            T = Ts[idx]
            qis = np.array([conf.get_partition_function(T) for conf in conformers])
            Q[idx] = np.sum(qis * np.exp(- Uis / constants.R / T))
            
    return E0_ref, Q


def calc_ms_tst_rate_coefficient(reactants_path,
                                products_path,
                                ts_path,
                                Temps,
                                freq_scale,
                                ):

    reactants = [get_rmg_conformer_from_yml(path, freq_scale) for path in reactants_path]
    products = [get_rmg_conformer_from_yml(path, freq_scale) for path in products_path]
    ts = get_rmg_conformer_from_yml(ts_path, freq_scale)

    E0_reac = 0
    Q_reac = 1
    for reactant in reactants:
        E0_min, Q = get_E0min_and_Q(Temps, reactant)
        E0_reac += E0_min
        Q_reac *= Q
        
    E0_prod = 0
    for product in products:
        E0_min, _ = get_E0min_and_Q(Temps, product)
        E0_prod += E0_min

    E0_TS, Q_TS = get_E0min_and_Q(Temps, ts)
    low_idx = np.argmin([conf.E0.value_si for conf in ts])
    low_ts_conf = list(read_yaml_file(ts_path)['conformers'].values())[low_idx]
    neg_frequency = (low_ts_conf['negative_frequencies'][0], 'cm^-1')

    eckart = Eckart(frequency = neg_frequency,
                    E0_reac=(E0_reac, 'J/mol'),
                    E0_TS=(E0_TS, 'J/mol'),
                    E0_prod=(E0_prod, 'J/mol'))


    dE0 = E0_TS - E0_reac

    ks_w_tunnel = np.zeros_like(Temps)
    ks_wo_tunnel = np.zeros_like(Temps)

    for idx in range(Temps.shape[0]):
        T = Temps[idx]
        ks_wo_tunnel[idx] = (constants.kB * T / constants.h * Q_TS[idx] / Q_reac[idx]) * np.exp(-dE0 / constants.R / T)
        ks_w_tunnel[idx] = ks_wo_tunnel[idx] * eckart.calculate_tunneling_factor(T)
    return (ks_w_tunnel, ks_wo_tunnel)


def parse_command_line_arguments(command_line_args=None):
    """
    Parse command-line arguments.

    Args:
        command_line_args: The command line arguments.

    Returns:
        The parsed command-line arguments by key words.
    """

    parser = argparse.ArgumentParser(description='Automatic Conformer Search (ACS)')
    parser.add_argument('--r1', 
                        type=str, 
                        default=None,
                        help='Final project info file name for reactant 1')
    parser.add_argument('--r2', 
                        type=str, 
                        default=None,
                        help='Final project info file name for reactant 2')
    parser.add_argument('--p1', 
                        type=str, 
                        default=None,
                        help='Final project info file name for product 1')
    parser.add_argument('--p2', 
                        type=str, 
                        default=None,
                        help='Final project info file name for product 2')
    parser.add_argument('--ts', 
                        type=str, 
                        default=None,
                        help='Final project info file name for TS')
    parser.add_argument('--freq_scale', 
                        type=float, 
                        default=None,
                        help='Optional frequency scale factor to use')   
    parser.add_argument('--index', 
                        type=int, 
                        default=0,
                        help='Optional reaciton index to use')     
    args = parser.parse_args(command_line_args)

    return args

def main():
    """
    The analyze_cosmo_result executable function
    """
    # 0. Parse input
    args = parse_command_line_arguments()
    cwd = os.getcwd()
    reactant_list = list()
    product_list = list()
    if args.r1 is not None:
        reactant_list.append(args.r1)
    if args.r2 is not None:
        reactant_list.append(args.r2)
    if args.p1 is not None:
        product_list.append(args.p1)
    if args.p2 is not None:
        product_list.append(args.p2)
    freq_scale = args.freq_scale

    Temps = np.arange(275., 351., 1.)

    # 1. Compute MS-TST rate coefficient

    fwd_reactants_path = [os.path.join(cwd, f) for f in reactant_list]
    fwd_products_path = [os.path.join(cwd, f) for f in product_list]
    fwd_rxn_order = len(fwd_reactants_path)
    rev_rxn_order = len(fwd_products_path)
    ts_path = os.path.join(cwd, args.ts)

    fwd_reactant_names = ' + '.join(reactant_list)
    fwd_product_names = ' + '.join(product_list)

    fwd_ks_w_tunnel, fwd_ks_wo_tunnel = calc_ms_tst_rate_coefficient(reactants_path=fwd_reactants_path,
                                                                    products_path=fwd_products_path,
                                                                    ts_path=ts_path,
                                                                    Temps=Temps,
                                                                    freq_scale=freq_scale,)

    rev_ks_w_tunnel, rev_ks_wo_tunnel = calc_ms_tst_rate_coefficient(reactants_path=fwd_products_path,
                                                                products_path=fwd_reactants_path,
                                                                ts_path=ts_path,
                                                                Temps=Temps,
                                                                freq_scale=freq_scale,)

    # 1.2. Save Rate coefficient
    rmg_lib_str = '''entry(
index = {index},
label = "{r} <=> {p}",
kinetics = {kinetics}
longDesc = 
"""
Calculated by ACS using multiple-structure local-harmonic 
conventional transition state theory with Eckart tunneling 
(MS-LH-CTST/Eckart).

Optfreq: wb97xd/def2tzvp in vacuum (freq scale factor: {freq_scale})
SP: dlpno-ccsd(t)/def2-tzvp normalPNO + Cosmo-RS TZVPD-Fine 
Solvent: H2O:MeOH = 0.7:0.3 (mol%)
"""
    '''

    k_units_dict = {1: 's^-1', 2: 'm^3/(mol*s)', 3: 'm^6/(mol^2*s)'}

    # 1.2.1 forward
    fwd_ks_w_tunnel_fit = Arrhenius().fit_to_data(Temps, fwd_ks_w_tunnel, kunits=k_units_dict[fwd_rxn_order], three_params=False) 
    with open('fwd_ks_w_tunnel_fit.txt', 'w') as f:
        f.write(rmg_lib_str.format(index=args.index, r=fwd_reactant_names, p=fwd_product_names, kinetics=fwd_ks_w_tunnel_fit, freq_scale=freq_scale))

    fwd_ks_wo_tunnel_fit = Arrhenius().fit_to_data(Temps, fwd_ks_wo_tunnel, kunits=k_units_dict[fwd_rxn_order], three_params=False) 
    with open('fwd_ks_wo_tunnel_fit.txt', 'w') as f:
        f.write(rmg_lib_str.format(index=args.index, r=fwd_reactant_names, p=fwd_product_names, kinetics=fwd_ks_wo_tunnel_fit, freq_scale=freq_scale))

     # 1.2.2 reverse
    rev_ks_w_tunnel_fit = Arrhenius().fit_to_data(Temps, rev_ks_w_tunnel, kunits=k_units_dict[rev_rxn_order], three_params=False) 
    with open('rev_ks_w_tunnel_fit.txt', 'w') as f:
        f.write(rmg_lib_str.format(index=args.index, r=fwd_product_names, p=fwd_reactant_names, kinetics=rev_ks_w_tunnel_fit, freq_scale=freq_scale))

    rev_ks_wo_tunnel_fit = Arrhenius().fit_to_data(Temps, rev_ks_wo_tunnel, kunits=k_units_dict[rev_rxn_order], three_params=False) 
    with open('rev_ks_wo_tunnel_fit.txt', 'w') as f:
        f.write(rmg_lib_str.format(index=args.index, r=fwd_product_names, p=fwd_reactant_names, kinetics=rev_ks_wo_tunnel_fit, freq_scale=freq_scale))


if __name__ == '__main__':
    main()