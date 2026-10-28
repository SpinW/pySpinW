"""Loading files in general."""
import pathlib

from pyspinw.symmetry.supercell import TiledSupercell
from pyspinw.serialisation import load as load_spw
from pyspinw.cif import load_cif
from pyspinw.interface import load_tb2j
from pyspinw.symmetry.supercell import Supercell


def load(filename: str | pathlib.Path,
         cif_supercell: Supercell | None = None,
         cif_entry_index: int | None = None,
         verbose: bool = False):
    """Loading of data from a file, infers from extension.

    Supported formats:

    - ``.spw``: loads as pySpinW internal format.
    - ``.cif``: loads as Crystallographic Interchange Format.
    - For TB2J, provide a directory containing ``tb2j.pickle`` or ``exchange.out``.

    Parameters
    ----------
    filename : str or pathlib.Path
        Path to the file or TB2J directory to load.
    cif_supercell : Supercell or None, optional
        Supercell specification when loading a ``.cif`` file. Defaults to
        ``TiledSupercell()`` if ``None``.
    cif_entry_index : int or None, optional
        Entry index when loading a ``.cif`` file. Defaults to ``0`` if ``None``.
    verbose : bool, optional
        If ``True``, print diagnostic information when loading a ``.cif`` file.
        Default is ``False``.

    Returns
    -------
    Any
        Deserialised pySpinW object, `Structure`, or `Hamiltonian` depending on
        the input format.

    Raises
    ------
    ValueError
        If the file has no extension or an unknown extension.
    """
    # Check for TB2J
    filename_path = pathlib.Path(filename)

    if filename_path.is_dir():

        pickle = filename_path / "tb2j.pickle"
        if pickle.exists():
            return load_tb2j(pickle)

        out = filename_path / "exchange.out"
        if out.exists():
            return load_tb2j(out)

    # Extension based
    filename_string = str(filename_path.name)

    if filename_string.lower().endswith(".spw"):
        return load_spw(filename)

    elif filename_string.lower().endswith(".cif"):
        cif_supercell = TiledSupercell() if cif_supercell is None else cif_supercell
        cif_entry_index = 0 if cif_entry_index is None else cif_entry_index

        return load_cif(filename, supercell=cif_supercell, entry_index=cif_entry_index, verbose=verbose)

    else:
        extension = filename_string.rsplit(".", maxsplit=1)[-1]

        if len(filename_string) == 1:
            raise ValueError("File has no extension")

        raise ValueError(f"Unknown extension '{extension}'")
