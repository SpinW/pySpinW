""" Tests for cif reading """

import os
from pathlib import Path

import pytest

from pyspinw.cif import load_cif

cif_directory = Path(__file__).resolve().parent / "cif_files"
cif_filenames = [filename for filename in os.listdir(cif_directory) if filename.lower().endswith(".cif")]

@pytest.mark.parametrize("cif_filename", cif_filenames)
def test_cif_reading(cif_filename):

    load_cif(cif_directory / cif_filename)


