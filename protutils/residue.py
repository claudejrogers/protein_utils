from itertools import groupby
from operator import attrgetter
import matplotlib.pyplot as plt

import atom


class Residue(object):
    """Container for a residue
    """

    def __init__(self, res, nres, chain, atoms):
        self.res = res
        self.nres = nres
        self.chain = chain
        self.atoms = atoms
        self._dict = {a.atom: a for a in atoms}

    @classmethod
    def from_atoms(cls, atoms):
        g = groupby(atoms, attrgetter('nres', 'chain', 'res'))
        resi = [(key, it) for key, it in g]
        if len(resi) != 1:
            raise ValueError(
                '{0} can only contain one residue.'.format(
                    cls.__name
                )
            )
        (nres, chain, res), it = resi[0]
        return cls(res, nres, chain, list(it))

    def __getitem__(self, name):
        """Return *Atom object with selected name
        """
        return self._dict[name]

    def get_next(self):
        """Returns key for 'next' residue
        """
        return '{0}_{1}'.format(self.nres + 1, self.chain)

    def get_prev(self):
        """Returns key for 'previous' residue
        """
        return '{0}_{1}'.format(self.nres - 1, self.chain)


class Residues(object):

    def __init__(self, atoms):
        g = groupby(atoms, attrgetter('nres', 'chain', 'res'))
        self._dict = {
            '{0}_{1}'.format(nres, chain): Residue(
                res, nres, chain, list(it)
            ) for (nres, chain, res), it in g}
        self._list = [v for v in sorted(self._dict.itervalues(),
                                        key=attrgetter('nres', 'chain'))]

    def __getitem__(self, key):
        return self._dict[key]

    def _get_phi_atoms(self, residue):

        if isinstance(residue, Residue):
            this = residue
        elif isinstance(residue, str):
            this = self[residue]
        else:
            raise TypeError(
                'residue must be a string or Residue object,'
                ' not {0}'.format(type(residue))
            )
        prev = self[this.get_prev()]
        C_ = prev['C']
        N = this['N']
        CA = this['CA']
        C = this['C']
        return C_, N, CA, C

    def _get_psi_atoms(self, residue):
        if isinstance(residue, Residue):
            this = residue
        elif isinstance(residue, str):
            this = self[residue]
        else:
            raise TypeError(
                'residue must be a string or Residue object,'
                ' not {0}'.format(type(residue))
            )
        next_ = self[this.get_next()]
        N = this['N']
        CA = this['CA']
        C = this['C']
        N_ = next_['N']
        return N, CA, C, N_

    def phi(self, residue):
        atoms = self._get_phi_atoms(residue)
        return atom.Atom.dihedral(*atoms)

    def psi(self, residue):
        atoms = self._get_psi_atoms(residue)
        return atom.Atom.dihedral(*atoms)

    def ramachandran_plot(self):
        x, y = [], []
        for r in self._list:
            try:
                phi = self.phi(r)
                psi = self.psi(r)
                x.append(phi)
                y.append(psi)
            except:
                continue
        plt.plot(x, y, 'ko')
        plt.xlim((-180, 180))
        plt.ylim((-180, 180))
        plt.xlabel(r'$\phi$', fontsize=22)
        plt.ylabel(r'$\psi$', fontsize=22)
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.show()
