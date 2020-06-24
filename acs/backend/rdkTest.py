#!/usr/bin/env python3

"""
This module contains unit tests of the :mod:`ccg.backend` module.
"""

import unittest

import numpy as np
import pybel
from rdkit import Chem

from ccg.backend import rdk


class RDKitMolTest(unittest.TestCase):
    """
    Contains unit tests of RDKitMol Class
    """
    def setUp(self):
        self.smiles = 'C'
        rd_mol = Chem.MolFromSmiles('C')
        self.rd_mol = Chem.AddHs(rd_mol)
        self.vdw_mat = np.array([[0., 0., 0., 0., 0.],
                                [0., 0., 0.96, 0.96, 0.96],
                                [0., 0., 0., 0.96, 0.96],
                                [0., 0., 0., 0., 0.96],
                                [0., 0., 0., 0., 0.]])
        self.xyz = '''14

O       1.85073500   -1.14075800   -0.80808500
O       1.64413700   -0.98958800    0.61920900
C       0.89993200    0.19167500    0.83332300
C       1.57997200    1.46553500    0.32458600
C       0.87059800    2.73251700    0.76632100
O       2.94131500    1.54377500    0.86706900
O       3.78534700    0.75777800    0.22699800
H       2.74479300   -0.76185300   -0.87545700
H       0.78501500    0.22753800    1.92086600
H      -0.09098800    0.11998200    0.36665200
H       1.69085200    1.41865800   -0.75897600
H       0.84097800    2.80358600    1.85617500
H      -0.15579100    2.73441900    0.39232200
H       1.37972400    3.61444100    0.37512000
'''
        self.ob_mol = pybel.readstring('xyz', self.xyz).OBMol
        self.rdkitmol = rdk.RDKitMol.FromOBMol(self.ob_mol)
        self.torsions = [[7, 0, 1, 2], [0, 1, 2, 3], [1, 2, 3, 4],
                         [2, 3, 4, 11], [2, 3, 5, 6]]

    def test_init(self):
        methane = rdk.RDKitMol(self.rd_mol)
        self.assertEqual(methane.GetNumAtoms(), 5)
        self.assertEqual(methane.GetNumBonds(), 4)
        atomic_nums = [atom.GetAtomicNum() for atom in methane.GetAtoms()]
        self.assertEqual(atomic_nums.count(6), 1)
        self.assertEqual(atomic_nums.count(1), 4)

    def test_FromSmiles(self):
        methane = rdk.RDKitMol.FromSmiles(self.smiles)
        self.assertEqual(methane.GetNumAtoms(), 5)
        self.assertEqual(methane.GetNumBonds(), 4)
        atomic_nums = [atom.GetAtomicNum() for atom in methane.GetAtoms()]
        self.assertEqual(atomic_nums.count(6), 1)
        self.assertEqual(atomic_nums.count(1), 4)

    def test_FromOBMol(self):
        self.assertEqual(self.rdkitmol.GetNumAtoms(), 14)
        self.assertEqual(self.rdkitmol.GetNumBonds(), 13)
        atomic_nums = [atom.GetAtomicNum() for atom in self.rdkitmol.GetAtoms()]
        self.assertEqual(atomic_nums.count(6), 3)
        self.assertEqual(atomic_nums.count(8), 4)
        self.assertEqual(atomic_nums.count(1), 7)

    def test_GetTorsionalModes(self):
        # TODO: Needs more examples
        self.assertEqual(self.rdkitmol.GetTorsionalModes(), self.torsions)

    def test_EmbedConformer(self):
        self.rdkitmol.EmbedConformer()
        self.assertEqual(self.rdkitmol.GetNumConformers(), 1)

    def test_EmbedMultipleConfs(self):
        self.rdkitmol.EmbedMultipleConfs(10)
        self.assertEqual(self.rdkitmol.GetNumConformers(), 10)

    def test_GetConformer(self):
        self.rdkitmol.RemoveAllConformers()
        with self.assertRaises(ValueError):
            self.rdkitmol.GetConformer()
        self.rdkitmol.EmbedConformer()
        conf = self.rdkitmol.GetConformer()
        self.assertTrue(isinstance(conf, rdk.RDKitConf))

    def test_GetConformers(self):
        self.rdkitmol.RemoveAllConformers()
        with self.assertRaises(ValueError):
            self.rdkitmol.GetConformers([1, 5])
        self.rdkitmol.EmbedMultipleConfs(8)
        confs = self.rdkitmol.GetConformers([1, 5])
        self.assertEqual(len(confs), 2)
        for conf in confs:
            self.assertTrue(isinstance(conf, rdk.RDKitConf))
            self.assertIn(conf.GetId(), [1, 5])

    def test_GetAllConformers(self):
        self.rdkitmol.RemoveAllConformers()
        confs = self.rdkitmol.GetAllConformers()
        self.assertEqual(len(confs), 0)
        self.rdkitmol.EmbedMultipleConfs(8)
        confs = self.rdkitmol.GetAllConformers()
        self.assertEqual(len(confs), 8)

    def test_ToRDMol(self):
        rd_mol = self.rdkitmol.ToRDMol()
        self.assertTrue(isinstance(rd_mol, Chem.rdchem.Mol))
        self.assertEqual(rd_mol.GetNumAtoms(), 14)
        self.assertEqual(rd_mol.GetNumBonds(), 13)
        atomic_nums = [atom.GetAtomicNum()
                       for atom in rd_mol.GetAtoms()]
        self.assertEqual(atomic_nums.count(6), 3)
        self.assertEqual(atomic_nums.count(8), 4)
        self.assertEqual(atomic_nums.count(1), 7)

    def test_GetVdwMatrix(self):
        methane = rdk.RDKitMol(self.rd_mol)
        methane.SetVdwMatrix()
        m = methane.GetVdwMatrix()
        self.assertTrue(np.allclose(m,
                                    self.vdw_mat))

    def test_SetVdwMatrix(self):
        methane = rdk.RDKitMol(self.rd_mol)
        methane.SetVdwMatrix(threshold=0.4)
        m1 = methane.GetVdwMatrix()
        self.assertTrue(np.allclose(m1, self.vdw_mat))
        methane.SetVdwMatrix(threshold=0.5)
        m2 = methane.GetVdwMatrix()
        self.assertTrue(np.all(m1 <= m2))
        self.assertTrue(np.allclose(m1 / 4 * 5, m2))


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
