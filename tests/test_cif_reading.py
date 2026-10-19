""" Tests for cif reading """

import os
from pathlib import Path

import pytest

from pyspinw import load, load_cif

cif_directory = Path(__file__).resolve().parent / "cif_files"
cif_filenames = [filename for filename in os.listdir(cif_directory) if filename.lower().endswith(".cif")]

@pytest.mark.parametrize("loader", [load, load_cif], ids=["load", "load_cif"])
@pytest.mark.parametrize("cif_filename", cif_filenames)
def test_cif_reading(cif_filename, loader):

    loader(cif_directory / cif_filename)


