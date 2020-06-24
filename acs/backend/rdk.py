#!/usr/bin/env python3

"""
This module provides class and methods for dealing with RDKit molecule.
"""

import logging
from typing import List, Optional, Union

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms as rdMT
from rdkit.Chem.rdchem import BondType, Conformer, Mol

# openbabel import is currently put behind rdkit.
# This relates to https://github.com/rdkit/rdkit/issues/2292
# For Mac, with conda build:
# rdkit -  2020.03.2.0 from rdkit channel
# openbabel -  2.4.1 from rmg channel
# works fine without sequence problem
# For Linux,
# rdkit - 2020.03.3 from rdkit channel
# openbabel - 2.4.1 from rmg channel does not work
import openbabel as ob

# Bond order dictionary for RDKit, numbers are the bond order.
# Note: There is a bond type 'OTHER' which may be helpful to treat TSs
ORDERS = {1: BondType.SINGLE, 2: BondType.DOUBLE, 3: BondType.TRIPLE, 1.5: BondType.AROMATIC,
          4: BondType.QUADRUPLE,
          'S': BondType.SINGLE, 'D': BondType.DOUBLE, 'T': BondType.TRIPLE, 'B': BondType.AROMATIC,
          'Q': BondType.QUADRUPLE}

# Van de Waals radii dictionary.
# The key refers to the symbol of the element or the atomic index
# obtained from Wikipedia
VDW_RADII = {'H': 1.09, 'C': 1.7, 'N': 1.55, 'O': 1.52,
             1: 1.2, 6: 1.7, 7: 1.55, 8: 1.52}

# The rotational bond definition in RDkit
# It is the same as rdkit.Chem.Lipinski import RotatableBondSmarts
ROTATABLE_BOND_SMARTS = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')

# Keep the representation method from rdchem.Mol
KEEP_RDMOL_ATTRIBUTES = ['_repr_html_',
                         '_repr_png_',
                         '_repr_svg_']


class RDKitMol(object):
    """
    A helpful wrapper for rdchem.Mol.
    The method nomenclature follows the Camel style to be consistent with RDKit.
    It keeps almost all of the orignal method of Chem.rdchem.Mol, but add few useful
    shortcuts, so that one doesn't need to import other functions to do the operation.
    """
    def __init__(self, rd_mol):
        """
        Generate an RDKitMol Molecule instance from a RDKit Chem.rdchem.Mol molecule.

        Args:
            rd_mol (Chem.rdchem.Mol): The RDKit Chem.rdchem.Mol molecule to be converted.
        """
        # keep the link to original rdchem.Mol so we can easily recover it if needed.
        self._rd_mol = rd_mol
        # Link methods of rdchem.Mol to the new instance
        for attr in dir(self._rd_mol):
            # Not reset private properties and repeated properties
            if not attr.startswith('_') and not hasattr(self, attr):
                setattr(self, attr, getattr(self._rd_mol, attr,))
            elif attr in KEEP_RDMOL_ATTRIBUTES:
                setattr(self, attr, getattr(self._rd_mol, attr,))

    @ classmethod
    def FromSmiles(cls,
                   smiles: str,
                   remove_h: bool = False,
                   sanitize: bool = True,
                  ) -> 'RDKitMol':
        """
        Convert a smiles to an RDkit Mol object.

        Args:
            smiles (str): A SMILES representation of the molecule.
            remove_h (bool, optional): Whether to remove hydrogen atoms from the molecule, ``True`` to remove.
            sanitize (bool, optional): Whether to sanitize the RDKit molecule, ``True`` to sanitize.

        Returns:
            RDKitMol: An RDKit molecule object corresponding to the SMILES.
        """
        rd_mol = Chem.MolFromSmiles(smiles)
        if not remove_h:
            rd_mol = Chem.AddHs(rd_mol)
        if sanitize:
            Chem.SanitizeMol(rd_mol)
        return cls(rd_mol)

    @ classmethod
    def FromRMGMol(cls,
                   rmg_mol: 'rmgpy.molecule.Molecule',
                   remove_h: bool = False,
                   sanitize: bool = True,
                   ) -> 'RDKitMol':
        """
        Convert a RMG molecular structure to an RDKit Mol object. Uses
        `RDKit <http://rdkit.org/>`_ to perform the conversion.
        Perceives aromaticity.
        Adopted from rmgpy/molecule/converter.py

        Args:
            rmg_mol (Molecule): An RMG Molecule object for the conversion.
            remove_h (bool, optional): Whether to remove hydrogen atoms from the molecule, ``True`` to remove.
            sanitize (bool, optional): Whether to sanitize the RDKit molecule, ``True`` to sanitize.

        Returns:
            RDKitMol: An RDKit molecule object corresponding to the input RMG Molecule object.
        """
        atom_id_map = dict()

        # only manipulate a copy of ``mol``
        mol_copy = rmg_mol.copy(deep=True)
        if not mol_copy.atom_ids_valid():
            mol_copy.assign_atom_ids()
        for i, atom in enumerate(mol_copy.atoms):
            # keeps the original atom order before sorting
            atom_id_map[atom.id] = i
        # mol_copy.sort_atoms()  # Sort the atoms before converting to ensure output is consistent between different runs
        atoms_copy = mol_copy.vertices

        rd_mol = Chem.rdchem.EditableMol(Chem.rdchem.Mol())
        for rmg_atom in atoms_copy:
            rd_atom = Chem.rdchem.Atom(rmg_atom.element.symbol)
            if rmg_atom.element.isotope != -1:
                rd_atom.SetIsotope(rmg_atom.element.isotope)
            rd_atom.SetNumRadicalElectrons(rmg_atom.radical_electrons)
            rd_atom.SetFormalCharge(rmg_atom.charge)
            if rmg_atom.element.symbol == 'C' and rmg_atom.lone_pairs == 1 and mol_copy.multiplicity == 1:
                # hard coding for carbenes
                rd_atom.SetNumRadicalElectrons(2)
            if not (remove_h and rmg_atom.symbol == 'H'):
                rd_mol.AddAtom(rd_atom)

        # Add the bonds
        for atom1 in atoms_copy:
            for atom2, bond12 in atom1.edges.items():
                if bond12.is_hydrogen_bond():
                    continue
                if atoms_copy.index(atom1) < atoms_copy.index(atom2):
                    rd_mol.AddBond(
                        atom_id_map[atom1.id], atom_id_map[atom2.id], orders[bond12.get_order_str()])

        # Make editable mol and rectify the molecule
        rd_mol = rd_mol.GetMol()
        if sanitize:
            Chem.SanitizeMol(rd_mol)
        if remove_h:
            rd_mol = Chem.RemoveHs(rd_mol, sanitize=sanitize)
        return cls(rd_mol)

    @ classmethod
    def FromOBMol(cls,
                  ob_mol: 'openbabel.OBMol',
                  remove_h: bool = False,
                  sanitize: bool = True,
                  ) -> 'RDKitMol':
        """
        Convert a OpenBabel molecular structure to an RDKit RDMol object.

        Args:
            ob_mol (Molecule): An OpenBabel Molecule object for the conversion.
            remove_h (bool, optional): Whether to remove hydrogen atoms from the molecule, ``True`` to remove.
            sanitize (bool, optional): Whether to sanitize the RDKit molecule, ``True`` to sanitize.

        Returns:
            RDKitMol: An RDKit molecule object corresponding to the input OpenBabel Molecule object.
        """
        rd_mol = Chem.rdchem.EditableMol(Chem.rdchem.Mol())
        for obatom in ob.OBMolAtomIter(ob_mol):
            rd_atom = Chem.rdchem.Atom(obatom.GetAtomicNum())
            isotope = obatom.GetIsotope()
            if isotope != 0:
                rd_atom.SetIsotope(isotope)
            spin = obatom.GetSpinMultiplicity()
            if spin == 2:  # radical
                rd_atom.SetNumRadicalElectrons(1)
            elif spin in [1, 3]:  # carbene
                # TODO: Not sure if singlet and triplet are distinguished
                rd_atom.SetNumRadicalElectrons(2)
            rd_atom.SetFormalCharge(obatom.GetFormalCharge())
            if not (remove_h and obatom.GetAtomicNum == 1):
                rd_mol.AddAtom(rd_atom)

        for bond in ob.OBMolBondIter(ob_mol):
            # Atom indexes in Openbael is 1-indexed, so we need to convert them to 0-indexed
            atom1_ind = bond.GetBeginAtomIdx() - 1
            atom2_ind = bond.GetEndAtomIdx() - 1
            # Get the bond order. For aromatic molecules, the bond order is not
            # 1.5 but 1 or 2. Manually set them to 1.5
            bond_order = bond.GetBondOrder() if not bond.IsAromatic() else 1.5

            rd_mol.AddBond(atom1_ind, atom2_ind, ORDERS[bond_order])

        rd_mol = rd_mol.GetMol()
        if sanitize:
            Chem.SanitizeMol(rd_mol)
        if remove_h:
            rd_mol = Chem.RemoveHs(rd_mol, sanitize=sanitize)
        return cls(rd_mol)

    @classmethod
    def FromRDMol(cls,
                  rd_mol: 'Chem.rdchem.Mol',
                  ) -> 'RDKitMol':
        """
        Convert a RDKit Chem.rdchem.Mol molecule to RDKitMol Molecule.

        Args:
            rd_mol (Chem.rdchem.Mol): The RDKit Chem.rdchem.Mol molecule to be converted.

        Returns:
            RDKitMol: An RDKitMol molecule.
        """
        return cls(rd_mol)

    def ToRDMol(self):
        """
        Convert the RDKitMol Molecule back to a RDKit Chem.rdchem.Mol.

        returns:
            Chem.rdchem.Mol: A RDKit Chem.rdchem.Mol molecule.
        """
        return self._rd_mol

    def GetTorsionalModes(self) -> list:
        """
        Get all of the torsional modes (rotors) from the molecule.

        Returns:
            list: A list of four-atom-indice to indicating the torsional modes.
        """
        return find_internal_torsions(self._rd_mol)

    def EmbedConformer(self):
        """
        Embed a conformer to the RDKitMol. This will overwrite current conformers.
        """
        AllChem.EmbedMolecule(self._rd_mol)

    def EmbedMultipleConfs(self, n: int = 1):
        """
        Embed conformers to the RDKitMol. This will overwrite current conformers.

        Args:
            n (int): The number of conformers to be embedded. The default is 1.
        """
        AllChem.EmbedMultipleConfs(self._rd_mol, numConfs=n)

    def GetConformer(self,
                     id: int = 0) -> 'RDKitConf':
        """
        Get the embedded conformer according to ID.

        Args:
            id (int): The ID of the conformer to be obtained. The default is 0.

        Raises:
            ValueError: Bad id assigned.

        Returns:
            RDKitConf: A conformer corresponding to the ID.
        """
        try:
            rd_org_conformer = self._rd_mol.GetConformer(id)
        except ValueError as e:
            raise ValueError(f"{e}: {id}")
        rd_conformer = RDKitConf(rd_org_conformer)
        rd_conformer.SetOwningMol(self)
        return rd_conformer

    def GetConformers(self,
                      ids: Union[list, tuple] = [0],
                      ) -> List['RDKitConf']:
        """
        Get the embedded conformers according to IDs.

        Args:
            ids (Union[list, tuple]): The ids of the conformer to be obtained.
                                      The default is [0].

        Raises:
            ValueError: Bad id assigned.

        Returns:
            List[RDKitConf]: A list of conformers corresponding to the IDs.
        """
        conformers = list()
        for id in ids:
            try:
                rd_org_conformer = self.GetConformer(id)
                conformers.append(rd_org_conformer)
            except ValueError as e:
                raise
        return conformers

    def GetAllConformers(self) -> List['RDKitConf']:
        """
        Get all of the embedded conformers.

        Returns:
            List['RDKitConf']: A list all of conformers.
        """
        return self.GetConformers(list(range(self.GetNumConformers())))

    def GetVdwMatrix(self) -> Optional[np.ndarray]:
        """
        Get the derived Van de Waals matrix. The dervied Van de Waals matrix can be used to analyze
        the collision of atoms. More information can be found from ``generate_vdw_mat``.

        Returns:
            Optional[np.ndarray]: A 2D array of the derived Van de Waals Matrix, if the
                                  the matrix exists, otherwise ``None``.
        """
        try:
            return self._vdw_mat
        except AttributeError:
            self.SetVdwMatrix()
            return self._vdw_mat

    def SetVdwMatrix(self,
                     threshold: float = 0.4,
                     vdw_radii: dict = VDW_RADII):
        """
        Set the derived Van de Waals matrix. The derived Van de Waals matrix is a upper
        triangle matrixcalculated from a threshold usually around 0.4 and Van de Waals
        Radii. Its diagonal elements are all zeros. The element (i, j) is calculated by
        threshold * sum( R(atom i) + R(atom j) ). If two atoms are bonded, the value is
        set to be zero. When threshold = 0.4, the value is close to the covalent bond
        length.

        Args:
            threshold (float): The threshold used to calculate the derived Van de Waals
                               matrix. A larger value results in a matrix with larger values;
                               When compared with distance matrix, it may overestiate the
                               overlapping between atoms. The default value is 0.4.
            vdw_radii (dict): A dict stores the Van de Waals radii of different elements.

        Raises:
            ValueError: Invalid threshold is supplied.
        """
        self._vdw_mat = generate_vdw_mat(self, threshold, vdw_radii)

    def GetDistanceMatrix(self, id: int = 0) -> np.ndarray:
        return Chem.rdmolops.Get3DDistanceMatrix(self._rd_mol, confId=id)

def find_internal_torsions(rd_mol) -> list:
    """
    Find the internal torsions from RDkit molecule.

    Args:
        rd_mol (rdkit.Chem.rdchem.Mol): RDKit molecule.

    Returns:
        list: A list of internal torsions.
    """
    torsions = list()
    rot_atom_pairs = rd_mol.GetSubstructMatches(ROTATABLE_BOND_SMARTS)

    for atoms_ind in rot_atom_pairs:
        pivots = [rd_mol.GetAtomWithIdx(i) for i in atoms_ind]
        first_atom_ind = determine_smallest_atom_index_in_scan(*pivots)
        pivots.reverse()
        last_atom_ind = determine_smallest_atom_index_in_scan(*pivots)
        torsions.append([first_atom_ind, *atoms_ind, last_atom_ind])
    return torsions


def determine_smallest_atom_index_in_scan(atom1: 'rdkit.Chem.rdchem.Atom',
                                          atom2: 'rdkit.Chem.rdchem.Atom',
                                          ) -> int:
    """
    Determine the smallest atom index in mol connected to ``atom1`` which is not ``atom2``.
    Returns a heavy atom if available, otherwise a hydrogen atom.
    Useful for deterministically determining the indices of four atom in a scan.
    This function assumes there ARE additional atoms connected to ``atom1``, and that ``atom2`` is not a hydrogen atom.

    Args:
        atom1 (Atom): The atom who's neighbors will be searched.
        atom2 (Atom): An atom connected to ``atom1`` to exclude (a pivotal atom).

    Returns:
        int: The smallest atom index (1-indexed) connected to ``atom1`` which is not ``atom2``.
    """
    neighbor = [a for a in atom1.GetNeighbors() if a.GetIdx()
                != atom2.GetIdx()]
    atomic_num_list = sorted([nb.GetAtomicNum() for nb in neighbor])
    min_atomic, max_atomic = atomic_num_list[0], atomic_num_list[-1]
    if min_atomic == max_atomic or min_atomic > 1:
        return min([nb.GetIdx() for nb in neighbor])
    else:
        return min([nb.GetIdx() for nb in neighbor if nb.GetAtomicNum() != 1])


def generate_vdw_mat(rd_mol,
                     threshold: float = 0.4,
                     vdw_radii: dict = VDW_RADII):
    """
    Generate a derived Van de Waals matrix. The derived Van de Waals matrix is a upper
    triangle matrixcalculated from a threshold usually around 0.4 and Van de Waals
    Radii. Its diagonal elements are all zeros. The element (i, j) is calculated by
    threshold * sum( R(atom i) + R(atom j) ). If two atoms are bonded, the value is
    set to be zero. When threshold = 0.4, the value is close to the covalent bond
    length.

    Args:
        threshold (float): The threshold used to calculate the derived Van de Waals
                            matrix. A larger value results in a matrix with larger values;
                            When compared with distance matrix, it may overestiate the
                            overlapping between atoms. The default value is 0.4.
        vdw_radii (dict): A dict stores the Van de Waals radii of different elements.

    Raises:
        ValueError: Invalid threshold is supplied.
    """
    if threshold <= 0:
        raise ValueError("The provided threshold is invalid.")

    # Initialize a vdw matrix
    num_atom = rd_mol.GetNumAtoms()
    vdw_mat = np.zeros((num_atom, num_atom))
    # Get all of the atom index
    atom_idx_list = range(num_atom)

    for atom1_ind in atom_idx_list:

        atom1 = rd_mol.GetAtomWithIdx(atom1_ind)
        bonded_atom_number = [nb.GetIdx() for nb in atom1.GetNeighbors()]

        for atom2_ind in atom_idx_list[atom1_ind + 1:]:
            if atom2_ind in bonded_atom_number:
                    # Set a ridiculously small number to bonded atoms,
                    # So that they won't arise the collision detector
                    vdw_mat[atom1_ind, atom2_ind] = 0.
            else:
                atom2 = rd_mol.GetAtomWithIdx(atom2_ind)
                vdw_mat[atom1_ind, atom2_ind] = threshold * \
                    (vdw_radii[atom1.GetAtomicNum()] +
                        vdw_radii[atom2.GetAtomicNum()])
    return vdw_mat
