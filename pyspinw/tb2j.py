""" Reads TB2J output files and returns Hamiltonian object """
import pickle
import numpy as np
import ase
import re

CELL_EXPR = ''.join([r'Cell \(Angstrom\):\s*'] + [r'([\-\d\.]*)\s*']*9)
ATOMS_EXPR = ''.join([r'^(\w*)\s*'] + [r'([\-\d\.]*)\s*']*7)
JPOS_EXPR = r'----\s*\n\s*([A-Za-z0-9]*)\s*([A-Za-z0-9]*)\s*\(\s*([\-\d]*),\s*([\-\d]*),\s*([\-\d]*)\)' \
                        r'\s*([\d\-\.]*)\s*\(\s*([\-\d\.]*),\s*([\-\d\.]*),\s*([\-\d\.]*)\)\s*([\d\-\.]*)'
JISO_EXPR = r'J_iso:\s*([\d\-\.]*)\s*$'
JDMI_EXPR = r'DMI:\s*\(\s*([\d\-\.]*)\s*([\d\-\.]*)\s*([\d\-\.]*)\)\s*$'
JANI_EXPR = r'J_ani:\s*$\s*\[\[\s*([\d\-\.]*)\s*([\d\-\.]*)\s*([\d\-\.]*)\s*\]\s*$' \
                        r'\s*\[\s*([\d\-\.]*)\s*([\d\-\.]*)\s*([\d\-\.]*)\s*\]\s*$' \
                        r'\s*\[\s*([\d\-\.]*)\s*([\d\-\.]*)\s*([\d\-\.]*)\s*\]\]\s*$'

class _sym_index:
    def __init__(self):
        self.syms = {}
    def __call__(self, sym):
        self.syms[sym] = (self.syms[sym] + 1) if sym in self.syms else 1
        return sym + str(self.syms[sym])

class TB2J_Input():

    def __init__(self, filename=None):
        if filename is not None:
            self.file = filename
            if filename.endswith('pickle'):
                with open(filename, 'rb') as f:
                    self._parse_pickle(pickle.load(f))
            else:
                with open(filename, 'r') as f:
                    self._parse_text(f.read())

    def _parse_text(self, text):
        self.type, self.data = ('text', text)
        if not all([kw in self.data for kw in ['Cell (Angstrom):', 'Atoms:', 'Exchange:']]):
            raise RuntimeError('Could not find one of the required keys, "Cell", "Atoms", or "Exchange"')
        cell_re = re.search(CELL_EXPR, self.data)
        cell = np.array([float(cell_re.group(i)) for i in range(1,10)]).reshape(3,3)
        atoms_re = re.finditer(ATOMS_EXPR, self.data.split('Atoms:')[1].split('Total')[0], re.MULTILINE)
        self.atoms, syms, pos, magmoms = [], [], [], []
        for at in atoms_re:
            if at.group(5):
                magmoms.append(list(map(float, at.groups()[5:8])) if at.group(8) else float(at.group(5)))
                self.atoms.append(at.group(1))
                syms.append(re.match(r'([A-Za-z]{1,2})', at.group(1)).group(1))
                pos.append(tuple(map(float, at.groups()[1:4])))
        self.struct = ase.Atoms(cell=cell, symbols=syms, positions=pos, magmoms=magmoms)
        exch_txt = self.data.split('Exchange:')[1]
        j_pos = [[v[0], v[1], list(map(int, v[2:5])), float(v[6]), list(map(float, v[7:10]))] for v in re.findall(JPOS_EXPR, exch_txt)]
        j_iso = np.array(list(map(float, re.findall(JISO_EXPR, exch_txt, re.MULTILINE))))
        j_dmi = np.array([list(map(float, v)) for v in re.findall(JDMI_EXPR, exch_txt, re.MULTILINE)])
        j_ani = np.array([np.array(list(map(float, v))).reshape(3,3) for v in re.findall(JANI_EXPR, exch_txt, re.MULTILINE)])
        if len(np.shape(magmoms)) > 1:
            self.noncolinear = True
            if j_iso.shape[0] != j_dmi.shape[0] or j_iso.shape[0] != j_ani.shape[0]:
                raise RuntimeError("Inconsistent number of exchanges")
        assert len(j_pos) == len(j_iso), "Inconsistent number of exchanges"
        self.exchanges = {'j_pos':j_pos, 'j_iso':j_iso, 'j_dmi':j_dmi, 'j_ani':j_ani}

    def _parse_pickle(self, data):
        self.type, self.data = ('pickle', data)
        self.struct = self.data['atoms']
        self.struct.set_initial_magnetic_moments(self.data['magmoms'])
        self.noncolinear = not self.data['colinear']
        sym_index = _sym_index()
        self.atoms = [sym_index(sym) for sym in self.struct.get_chemical_symbols()]
        d_exch, d_dist = (self.data[k] for k in ['exchange_Jdict', 'distance_dict'])
        ord_k = [k for k, v in d_dist.items() if k in d_exch]
        ord_k = [ord_k[i] for i in np.argsort([v[1] for k, v in d_dist.items() if k in d_exch])]
        j_pos = [[self.atoms[k[1]], self.atoms[k[2]], k[0], 1000*d_exch[k], d_dist[k][0], d_dist[k][1]] for k in ord_k]
        j_iso, j_dmi, j_ani = ([1000*d_exch[k] for k in ord_k], [], [])
        if self.noncolinear:
            j_dmi = [1000*self.data['dmi_ddict'][k] for k in ord_k]
            j_ani = [1000*self.data['Jani_dict'][k] for k in ord_k]
        self.exchanges = {'j_pos':j_pos, 'j_iso':j_iso, 'j_dmi':j_dmi, 'j_ani':j_ani}

    @classmethod
    def from_exchange_out(cls, filename):
        """ Read input from an exchange.out TB2J file """
        with open(filename, 'r') as f:
            return cls()._parse_text(f.read())

    @classmethod
    def from_text(cls, text):
        """ Read input from text string containing the contents of an exchange.out TB2J file """
        return cls()._parse_text(text)

    @classmethod
    def from_pickle(cls, filename):
        """ Read input from a TB2J.pickle file """
        with open(filename, 'rb') as f:
            return cls()._parse_pickle(pickle.load(f))

    def to_hamiltonian(self):
        """ Create a SpinW Hamiltonian object from the TB2J inputs """
        pass

if __name__ == '__main__':
    import sys
    sp = TB2J_Input(sys.argv[1])
