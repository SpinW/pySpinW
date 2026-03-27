""" System for dealing with windows parallelisation issues """

from multiprocessing import freeze_support
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import platform

windows_python_parallelisation_enabled = False

def get_Executor():
    """ Get an Executor appropriate for the current system """
    if platform.system() == "Windows":
        if not windows_python_parallelisation_enabled:
            return ThreadPoolExecutor

    return ProcessPoolExecutor



def set_up_windows_python_parallelisation():
    """ This needs to be run to enable python parallelisation on windows

    Call it like this at the start of your script:

    ```
    if __name__ == "__main__":
        set_up_windows_python_parallelisation()

        [your code]
    ```


    Applies freeze_support() and sets a flag to be used elsewhere

    """
    # ruff: noqa: PLW0603
    global windows_python_parallelisation_enabled

    windows_python_parallelisation_enabled = True

    freeze_support()
