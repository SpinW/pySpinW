"""Runs the scripts in the tutorial_inputs folder to ensure they run

We just want to check the syntax does not throw errors and do not
check for output values or graphs
"""

import os
import runpy
import pytest


def test_tutorial_docs_build(monkeypatch):
    """ Run the tutorials docs build"""
    path = os.path.join(os.path.dirname(__file__), '..', 'examples', 'tutorials')
    monkeypatch.chdir(path)
    monkeypatch.syspath_prepend(path)
    
    runpy.run_path(
        os.path.join(
            os.path.dirname(__file__),
            '..',
            'examples',
            'tutorials',
            'generate_tutorial_files.py'
        ), run_name="__main__")