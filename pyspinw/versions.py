from importlib.metadata import PackageNotFoundError, version
import tomllib
from pathlib import Path

def get_version():
    """ Get the current version of pyspinw"""
    try:
        v = version("pyspinw")
    except PackageNotFoundError:
        pyproject = Path(__file__).parents[1] / "Cargo.toml"
        data = tomllib.loads(pyproject.read_text())
        v = data["package"]["version"]

    return v