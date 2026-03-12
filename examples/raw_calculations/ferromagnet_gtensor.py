"""Basic Heisenberg ferromagnetic chain with magnetic field and g-tensor specified.

Equivalent MATLAB code:
```
% spin chain
FMchain = spinw;
FMchain.genlattice('lat_const',[3 8 8],'angled',[90 90 90])
FMchain.addatom('r', [0 0 0],'S', 1,'label','MCu1','color','blue')

% magnetic structure
FMchain.gencoupling('maxDistance',7)
FMchain.addmatrix('value',-eye(3),'label','Ja','color','green')
FMchain.addcoupling('mat','Ja','bond',1);
FMchain.genmagstr('mode','direct', 'k',[0 0 0],'n',[1 0 0],'S',[0; 1; 0]);

% magnetic field
FMchain.field([5 10 15]);
gmat = ones(3) + diag([1 2 3]);
FMchain.addmatrix('label', 'g0', 'value', gmat);
FMchain.addg('g0')

% create spectrum and plot
FMspec = FMchain.spinwave({[0 0 0] [1 0 0]},'hermit',false,'gtensor',true);
sw_plotspec(FMspec)
```
"""
import sys

import numpy as np

from examples.raw_calculations.utils import run_example, plot, py_classes

def ferromagnet_gtensor(n_q = 100, classes = py_classes):
    """Basic ferromagnet with magnetic field and g-tensor."""
    rust_kw = {'dtype':complex, 'order':'F'}
    Coupling = classes.coupling
    MagneticField = classes.field

    q_mags = np.linspace(0, 1, n_q).reshape(-1, 1)
    q_vectors = np.array([0, 1, 0]).reshape(1, 3) * q_mags
    positions = [np.array([0., 0., 0.])]

    # external field
    g_tensor = [np.array([[2., 1., 1.],
                          [1., 3., 1.],
                          [1., 1., 4.]], **rust_kw)]
    ext_field = MagneticField(vector=np.array([5., 10., 15.], **rust_kw), g_tensors=g_tensor)

    # Single site
    rotations = [np.eye(3, **rust_kw)]
    magnitudes = np.array([1.0])  # spin-1
    couplings = [Coupling(0, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([0., 1., 0.])),
                 Coupling(0, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([0., -1., 0.])), ]

    return rotations, magnitudes, q_vectors, couplings, positions, np.eye(3, order='F'), ext_field

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    if len(sys.argv) > 1:
        use_rust = "py" not in sys.argv[1]
    else:
        use_rust = True

    q_mags = np.linspace(0, 1, 100)

    _, energies, sqw = run_example(ferromagnet_gtensor, use_rust)

    # Note: we get complex data types with real part zero
    plot(q_mags, energies, sqw)

    # Compare with tutorial 1, 3rd last figure (https://spinw.org/tutorial1_05.png)
    plt.show()
