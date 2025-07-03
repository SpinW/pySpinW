# Example syntax to get energies and scattering for a simple ferromagnetic chain powder sample

the_only_site = Site([0.5, 0.5, 0.5], m=[1, 0, 0])

sites = [the_only_site]

symmetry = Symmetry(cell=UnitCell(1, 1, 1), magnetic_group="P1.1")

couplings = batch_coupling(the_only_site, symmetry, max_distance=1.5, direction=(0, 0, 1), j=-1)

hamiltonian = Hamiltonian(sites, couplings, anisotropy=None)

# Spaghetti diagram data

energies = hamiltonian.energies(path=symmetry.path("GX", n=100))

# Powder spectrum

sample = Powder(hamiltonian)

instrument = Instrument.from_ResINS("Let")

experiment = Experiment(sample, instrument)

spectrum = experiment.simulate1d(npts=201)


