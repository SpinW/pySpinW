import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.io import loadmat
from scipy.optimize import least_squares, Bounds
from pyspinw import *

""" Set up model for BiMn2O5 for powder fitting """
cif_struc = load_cif(os.path.join(os.path.dirname(__file__), 'data_files', 'ICSD_CollCode117638.cif'))
# Magnetic structure from PRB77 134434 Table IV
S0 = np.array([[2.10, -0.33,  0.25j], [ 2.10, -0.33, -0.25j], [ 2.07,  0.56,  0.08j], [ 2.07,  0.56, -0.08j],
               [-2.83, -0.23,  0.00], [-2.83,  0.33,  0.00], [ 2.80, -0.34,  0.00], [-2.74, -0.64,  0.00]])
# Note that CIF has no magnetic moment information and PySpinW objects are *immutable*
# So we have to set the magnetic moments in a copy of the sites list and create a new structure
sites = cif_struc.sites_by_name('Mn1') + cif_struc.sites_by_name('Mn2')
perm, scales = [3, 2, 0, 1, 7, 4, 5, 6], [1.5]*4 + [2]*4
for mom, idx, scale in zip(S0, perm, scales):
    sites[idx].spin = mom / np.linalg.norm(np.real(mom)) * scale
structure = Structure(sites, unit_cell=cif_struc.unit_cell, spacegroup=cif_struc.spacegroup,
                      supercell=SummationSupercell([CommensuratePropagationVector(0.5, 0, 0.5)]))

# Set the exchange interactions Jx, Kx using notation of PRB84 054444
# J1 and J2 along c, J3 (inter-chain), J4 (intra Mn3-Mn4), J5 (intra Mn3-Mn3)
Jvec = [-0.31, 0.45, 0.2, 3.544, 1.07, 0.59, 2.85]  # J1-J5, K1, K2
exchanges = generate_exchanges(structure, bond=1, j=Jvec[0], naming_pattern='J1') \
          + generate_exchanges(structure, bond=3, j=Jvec[1], naming_pattern='J2') \
          + generate_exchanges(structure, bond=4, j=Jvec[2], color=[1,0,0], naming_pattern='J3') \
          + generate_exchanges(structure, bond=5, j=Jvec[3], color=[0,1,0], naming_pattern='J4') \
          + generate_exchanges(structure, bond=2, j=Jvec[4], color=[0,0,1], naming_pattern='J5')

anisotropies = axis_anisotropies(sites[:4], a=Jvec[5], name="K1") \
             + axis_anisotropies(sites[4:], a=Jvec[6], name="K2")

hamiltonian = Hamiltonian(structure, exchanges, anisotropies)
hamiltonian.print_summary()
view(hamiltonian)

# Try to optimise the structure (currently buggy)
#optstruc = hamiltonian.ground_state(initial_randomisation='randomised', planar=hamiltonian.structure.sites)
#view(optstruc)
#optstruc.print_summary()

# Compute a "spaghetti plot"
path = Path([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1]])
hamiltonian.spaghetti_plot(path)

# Compute the powder spectrum
sample = Powder(hamiltonian)
path = Path1D(0.01, 4, n_points=101)
fig = sample.show_spectrum(path, n_energy_bins=100, n_samples=500, energy_stddev=0.4,
                           ignore_imaginary=True, show=False)
fig.axes[0].collections[0].set_clim(0, 5)
cb = plt.colorbar(fig.axes[0].collections[0])
plt.show()

# Loads experimental data
mat = loadmat(os.path.join(os.path.dirname(__file__), 'data_files', 'bimn2o5_ei30.mat'))
en, qm, dat, err = (np.squeeze(v) for v in [mat['x'][0][0], mat['x'][0][1], mat['y'], mat['e']])
path = Path1D(np.min(qm), np.max(qm), n_points=len(qm))

# Set up a fitting function - the parameters take the form of "<Name.ParameterName>"
# Where "Name" is the name used in their creation (e.g. J1, K1 etc above)
# and "ParameterName" is the name of the variable used in the constructor (e.g. "j" or "a")
# Note the use of n_samples=50 here to reduce running time (to ~10min)
specfunc = sample.parameterized_spectrum(parameters=["J1.j", "J2.j", "J4.j", "J5.j", "K1.a", "K2.a"],
                                         path=path, n_energy_bins=len(en), n_samples=50,
                                         min_energy=np.min(en), max_energy=np.max(en),
                                         energy_stddev=0.4, ignore_imaginary=True)
# Test that evaluating the fit function works
spectrum = specfunc(-0.31, 0.45, 3.544, 1.07, 0.59, 2.85)
spectrum[np.where(np.isnan(dat))] = np.nan
# Plot the spectrum over the data
fig, ax = plt.subplots()
mesh = ax.pcolormesh(qm, en, dat.T, shading='nearest', vmin=0, vmax=5)
Xi, Yi = np.meshgrid(qm, en)
ax.contour(Xi, Yi, spectrum.T, levels=np.linspace(0, 5, 10), linewidths=1.0, colors='k')
cb = fig.colorbar(mesh, ax=ax)
plt.show()

# Define mimization function which automatically fits a scale factor
nnan = np.where(~np.isnan(dat))
def minfun(x):
    calc = specfunc(*tuple(x))
    def scalefun(y):
        return (dat[nnan] - y*calc[nnan]).reshape(-1)
    y_opt = least_squares(scalefun, x0=100., bounds=Bounds([0.0], [1e9]), diff_step=1e-2)
    return (dat[nnan] - y_opt.x*calc[nnan]).reshape(-1)

def callback(intermediate_result):
    """Callback to print out what is going on"""
    print(intermediate_result)

fitpars = least_squares(minfun, x0=[-0.31, 0.45, 3.544, 1.07, 0.59, 2.85], callback=callback, diff_step=1e-2)
print(fitpars)

# Plot fitted spectrum
spectrum = specfunc(*tuple(fitpars.x))
fig, ax = plt.subplots()
mesh = ax.pcolormesh(qm, en, dat.T, shading='nearest', vmin=0, vmax=5)
Xi, Yi = np.meshgrid(qm, en)
ax.contour(Xi, Yi, spectrum.T, levels=np.linspace(0, 5, 10), linewidths=1.0, colors='k')
cb = fig.colorbar(mesh, ax=ax)
plt.show()
