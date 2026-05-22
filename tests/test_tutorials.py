"""Runs the scripts in the tutorials folder to ensure they run

We just want to check the syntax does not throw errors and do not
check for output values or graphs
"""

from unittest.mock import MagicMock, patch
import pytest
import sys
import os
import importlib
import builtins

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'examples', 'tutorials'))

@patch('matplotlib.pyplot.show', MagicMock())
@patch('pyspinw.view', MagicMock())
@pytest.mark.parametrize("tutorial",
                         ['antiferromagnetic_chain', 
                          'antiferromagnetic_chain_with_field',
                          'ferromagnetic_chain',
                          'ferromagnetic_chain_series',
                          'high_symmetry',
                          'kagome_antiferromagnet',
                          'kagome_ferromagnet',
                          'kagome_supercell',
                          'powder_fitting',
                          'powder_spectrum',
                          'reciprocal_space_map',
                          'square_antiferromagnet',
                          'swapping_chain_series',
                          'triangular_antiferro',
                          'triangular_antiferro_rotframe',
                          'twin_example',
                         ])
def test_tutorials(tutorial):
    """Test tutorials by importing them, but with graphics mocked out so it can run headless"""
    importlib.import_module(tutorial)
