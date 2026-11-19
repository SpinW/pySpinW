"""Mathematical constants used by pySpinW.

Attributes
----------
MU_B : float
    Bohr magneton in units of meV/T.
ELECTRON_G : float
    Electron g-factor (defined as 2.0 to agree with MATLAB SpinW).

LEVI_CIVITA : numpy array
     Levi-Civita symbol, such that
     $ (A\times B) = \epsilon_{ijk} A_j B_k $
"""



# Bohr magneton in units meV/T
MU_B = 0.057883818066

# g-factor for electrons (anomalous moment + 2)
ELECTRON_G = 2.#00231930436092 # Defined as 2 to agree with Matlab

LEVI_CIVITA = np.array([[[0, 0, 0], [0, 0, 1], [0, -1, 0]],
                        [[0, 0, -1], [0, 0, 0], [1, 0, 0]],
                        [[0, 1, 0], [-1, 0, 0], [0, 0, 0]]])