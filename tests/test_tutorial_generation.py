"""Runs the scripts in the tutorial_inputs folder to ensure they run

We just want to check the syntax does not throw errors and do not
check for output values or graphs
"""

import os
import runpy
import pytest
import shutil

import pyspinw

path = os.path.join(os.path.dirname(__file__), '..', 'examples', 'tutorials')
output_directory = os.path.join(path, "tutorial_outputs")

@pytest.fixture
def preserve_directory(tmp_path):
    """ We need to back up the tutorials directory"""
    directory = os.path.abspath(output_directory)
    backup = tmp_path / "tutorials_backup"

    if os.path.exists(directory):
        shutil.copytree(directory, backup)

    yield

    if os.path.exists(directory):
        shutil.rmtree(directory)

    shutil.copytree(backup, directory)

def test_tutorial_docs_build(preserve_directory, monkeypatch):
    """ Run the tutorials docs build"""

    monkeypatch.chdir(path)
    monkeypatch.syspath_prepend(path)

    # Patch out snapshot as running it on a runner causes all sorts of issues!!!!!
    monkeypatch.setattr(pyspinw, "snapshot", lambda *args, **kwargs: None)

    runpy.run_path(
        os.path.join(
            os.path.dirname(__file__),
            '..',
            'examples',
            'tutorials',
            'generate_tutorial_files.py'
        ), run_name="__main__")