""" Tests for parsing TB2J inputs """

import os
from pathlib import Path

import pytest

from pyspinw.path import Path as SQWPath
from pyspinw.tb2j import TB2J_Input
import numpy as np

tb2j_directory = Path(__file__).resolve().parent / "tb2j_files"
tb2j_subdirs = [dirname for dirname in os.listdir(tb2j_directory)]

@pytest.mark.parametrize("tb2j_subdir", tb2j_subdirs)
def test_tb2j_reading(tb2j_subdir):
    loader_txt = TB2J_Input(tb2j_directory / tb2j_subdir / "exchange.out")
    loader_pickle = TB2J_Input(tb2j_directory / tb2j_subdir / "TB2J.pickle")
    path = SQWPath([[0,0,1], [0.5,0.5,1]])
    en_txt, _ = loader_txt.to_hamiltonian().energies_and_intensities(path)
    en_pickle, _ = loader_pickle.to_hamiltonian().energies_and_intensities(path)
    np.testing.assert_allclose(en_txt, en_pickle, atol=1e-4, rtol=1e-2)
    np.testing.assert_allclose(en_txt[1],
        [-0.32918, -0.32888, 0, 0, 0, 0, 0, 0, 0, 0, 0.32888, 0.32918], atol=1e-5)
