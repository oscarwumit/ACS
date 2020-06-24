#!/usr/bin/env python3

"""
A module for performing various geometry-related format conversions.

This module is majorly based on the converter.py from ARC/arc/species
"""

from arkane.common import mass_by_symbol
import os
from typing import Optional


def xyz_str_to_xyz_dict(xyz_str: str) -> dict:
    """
    Convert a string xyz format to the ARC dict xyz style.
    Note: The ``xyz_str`` argument could also direct to a file path to parse the data from.
    The xyz string format may have optional Gaussian-style isotope specification, e.g.::

        C(Iso=13)    0.6616514836    0.4027481525   -0.4847382281
        N           -0.6039793084    0.6637270105    0.0671637135
        H           -1.4226865648   -0.4973210697   -0.2238712255
        H           -0.4993010635    0.6531020442    1.0853092315
        H           -2.2115796924   -0.4529256762    0.4144516252
        H           -1.8113671395   -0.3268900681   -1.1468957003

    which will also be parsed into the ARC xyz dictionary format, e.g.::

        {'symbols': ('C', 'N', 'H', 'H', 'H', 'H'),
         'isotopes': (13, 14, 1, 1, 1, 1),
         'coords': ((0.6616514836, 0.4027481525, -0.4847382281),
                    (-0.6039793084, 0.6637270105, 0.0671637135),
                    (-1.4226865648, -0.4973210697, -0.2238712255),
                    (-0.4993010635, 0.6531020442, 1.0853092315),
                    (-2.2115796924, -0.4529256762, 0.4144516252),
                    (-1.8113671395, -0.3268900681, -1.1468957003))}

    Args:
        xyz_str (str): The string xyz format to be converted.

    Returns:
        dict: The ARC xyz format.

    Raises:
        ConverterError: If xyz_str is not a string or does not have four space-separated entries per non empty line.
    """
    xyz_dict = {key: [] for key in ['symbols', 'isotopes', 'coords']}
    for line in xyz_str.strip().splitlines():
        splits = line.split()
        if len(splits) != 4:
            raise ValueError(f'Input has an incorrect format, expected 4 elements in each line, '
                             f'got "{line}" in:\n{xyz_str}')
        symbol = splits[0]
        if '(iso=' in symbol.lower():
            isotope = int(symbol.split('=')[1].strip(')'))
            symbol = symbol.split('(')[0]
        else:
            # no specific isotope is specified in str_xyz
            isotope = get_most_common_isotope_for_element(symbol)
        coord = [float(num) for num in splits[1:4]]
        xyz_dict['symbols'].append(symbol)
        xyz_dict['isotopes'].append(isotope)
        xyz_dict['coords'].append(coord)
    return xyz_dict


def xyz_str_to_xyz_file(xyz_str: str,
                        save_path: Optional[str] = None,
                        comment: str = '',
                        ) -> str:
    xyz_atom_list = xyz_str.splitlines()
    xyz_file = f'{len(xyz_atom_list)}\n{comment}\n'
    xyz_file += xyz_str

    if save_path:
        try:
            with open(save_path, 'w') as f:
                f.write(xyz_file)
        except TypeError:
            raise ValueError('Bad save path.')
    return xyz_file


def xyz_dict_to_xyz_str(xyz_dict: dict,
                        isotope_format: Optional[str] = None,
                        ) -> str:
    """
    Convert an ARC xyz dictionary format, e.g.::

        {'symbols': ('C', 'N', 'H', 'H', 'H', 'H'),
         'isotopes': (13, 14, 1, 1, 1, 1),
         'coords': ((0.6616514836, 0.4027481525, -0.4847382281),
                    (-0.6039793084, 0.6637270105, 0.0671637135),
                    (-1.4226865648, -0.4973210697, -0.2238712255),
                    (-0.4993010635, 0.6531020442, 1.0853092315),
                    (-2.2115796924, -0.4529256762, 0.4144516252),
                    (-1.8113671395, -0.3268900681, -1.1468957003))}

    to a string xyz format with optional Gaussian-style isotope specification, e.g.::

        C(Iso=13)    0.6616514836    0.4027481525   -0.4847382281
        N           -0.6039793084    0.6637270105    0.0671637135
        H           -1.4226865648   -0.4973210697   -0.2238712255
        H           -0.4993010635    0.6531020442    1.0853092315
        H           -2.2115796924   -0.4529256762    0.4144516252
        H           -1.8113671395   -0.3268900681   -1.1468957003

    Args:
        xyz_dict (dict): The ARC xyz format to be converted.
        isotope_format (str, optional): The format for specifying the isotope if it is not the most abundant one.
                                        By default, isotopes will not be specified. Currently the only supported
                                        option is 'gaussian'.

    Returns:
        str: The string xyz format.

    Raises:
        ConverterError: If input is not a dict or does not have all attributes.
    """

    recognized_isotope_formats = {
        'gaussian': lambda symbol, isotope: f'{symbol}Iso={isotope}'}
    if isotope_format and isotope_format not in recognized_isotope_formats:
        raise NotImplementedError

    xyz_str = str()
    for symbol, isotope, coord in zip(xyz_dict['symbols'], xyz_dict['isotopes'], xyz_dict['coords']):
        if isotope_format:
            common_isotope = get_most_common_isotope_for_element(symbol)
            symbol = symbol if common_isotope == isotope else \
                recognized_isotope_formats[isotope_format](symbol, isotope)
        xyz_str += f'{symbol:14}'
        xyz_str += '{0:14.8f}{1:14.8f}{2:14.8f}\n'.format(*coord)
    return xyz_str


def xyz_dict_to_xyz_file(xyz_dict: dict,
                         save_path: Optional[str] = None,
                         comment: str = '',
                         ) -> str:
    """
    Convert xyz dict to xyz_file.

    Args:
        xyz_dict (dict): The ARC xyz format to be converted.
        save_path (Optional[str]): The path to save the xyz file.
        comment (str): The comment append to the xyz file

    Returns:
        str: The xyz file format str.
    """
    xyz_file = f'{len(xyz_dict["symbols"])}\n{comment}\n'
    xyz_file += xyz_dict_to_xyz_str(xyz_dict)
    return xyz_file


def xyz_file_to_xyz_str(xyz_file: str) -> str:
    """
    Convert xyz file (or xyz file format string) to xyz str.

    Args:
        xyz_file (str): The path of ARC xyz file or xyz file format str.

    Returns:
        str: The xyz str.
    """
    if os.path.isfile(xyz_file):
        with open(xyz_file) as f:
            atom_list = f.readlines()
    elif isinstance(xyz_file, str):
        atom_list = xyz_file.splitlines()
    else:
        raise ValueError("The input does not satisfy a valid xyz file format.")

    try:
        xyz_str = '\n'.join(atom_list[2:])
    except IndexError:
        raise ValueError("The input does not satisfy a valid xyz file format.")

    return xyz_str


def xyz_file_to_xyz_dict(xyz_file: str) -> dict:
    """
    Convert xyz file (or xyz file format string) to xyz_dict.

    Args:
        xyz_file (str): The path of ARC xyz file or xyz file format str.

    Returns:
        dict: The xyz dict.
    """
    xyz_str = xyz_file_to_xyz_str(xyz_file)
    return xyz_str_to_xyz_dict(xyz_str)


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
