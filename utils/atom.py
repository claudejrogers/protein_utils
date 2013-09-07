import operator
from copy import deepcopy
import numpy as np

from .kabsch import aligned_rmsd, align


class Atom(object):
    """
    Base class for an atom description
    """
    def __init__(self, record, natom, atom, res, chain, nres, x, y, z, charge):
        self.record = record
        self.natom = natom
        self.atom = atom
        self.res = res
        self.chain = chain
        self.nres = nres
        self.x = x
        self.y = y
        self.z = z
        self.charge = charge
        self.coord = np.array([x, y, z])
        self.connections = []

    def add_connections_from_list(self, connections):
        self.connections.extend(connections)

    def distance(self, other):
        return np.linalg.norm(self.coord - other.coord)

    def __repr__(self):
        if self.record == 'ATOM' and not self.atom.startswith('H'):
            atom = ' {0}'.format(self.atom)
        else:
            atom = self.atom
        return "<{0}: {1:6s}{2:5d} {3:4s} {4:3s} {5:5d} ...>".format(
            self.__class__.__name__, self.record, atom, self.res, self.chian,
            self.nres
        )

    def __key(self):
        return "{0}{1}{2}{3}{4}".format(self.record, self.natom, self.atom,
                                        self.res, self.chain, self.nres)

    def __eq__(self, other):
        return self.__key() == other.__key()

    def __hash__(self):
        return hash(self.__key())


class AtomCollection(object):
    """
    Base class for an Atom collection
    """
    def __init__(self, atoms=None):
        self.atoms = atoms or []
        self.matrix = np.array([atom.coord for atom in self.atoms])

    def __len__(self):
        return len(self.atoms)

    def __iter__(self):
        return iter(self.atoms)

    def __getitem__(self, key):
        return self.atoms[key]

    def __eq__(self, other):
        if len(self) != len(other):
            return False
        return all(operator.eq(s, o) for s, o in zip(
            self.atoms, other.atoms
        ))

    def __contains__(self, other):
        if len(self) < len(other):
            return False
        atoms_ = set(self.atoms)
        for a in other.atoms:
            if not a in atoms_:
                return False
        return True

    def __and__(self, other):
        """Returns the common atoms between self and other

        >>> a = AtomCollection(atoms_1)
        >>> b = AtomCollection(atoms_2)
        >>> common_atoms = a & b

        """
        atom_set = set(self.atoms) & set(other.atoms)
        atoms = sorted(atom_set, key=operator.attrgetter('natom'))
        return self._init_with_atoms(atoms)

    def __or__(self, other):
        """Returns the intersection of atoms in self and other

        >>> a = AtomCollection(atoms_1)
        >>> b = AtomCollection(atoms_2)
        >>> all_atoms = a & b

        """
        atom_set = set(self.atoms) | set(other.atoms)
        atoms = sorted(atom_set, key=operator.attrgetter('natom'))
        return self._init_with_atoms(atoms)

    def __xor__(self, other):
        """Returns atoms in self that are not in other and the atoms in
        other not in self.

        >>> all_atoms = AtomCollection(atoms)
        >>> protons = all_atoms.select(atom__startswith='H')
        >>> removed_hydrogens = all_atoms ^ protons

        """
        atom_set = set(self.atoms) ^ set(other.atoms)
        atoms = sorted(atom_set, key=operator.attrgetter('natom'))
        return self._init_with_atoms(atoms)

    def __iand__(self, other):
        atom_set = set(self.atoms) & set(other.atoms)
        self.atoms = sorted(atom_set, key=operator.attrgetter('natom'))
        self.matrix = np.array([atom.coord for atom in self.atoms])
        return self

    def __ior__(self, other):
        atom_set = set(self.atoms) | set(other.atoms)
        self.atoms = sorted(atom_set, key=operator.attrgetter('natom'))
        self.matrix = np.array([atom.coord for atom in self.atoms])
        return self

    def __ixor__(self, other):
        atom_set = set(self.atoms) ^ set(other.atoms)
        self.atoms = sorted(atom_set, key=operator.attrgetter('natom'))
        self.matrix = np.array([atom.coord for atom in self.atoms])
        return self

    def _init_with_atoms(self, atoms, deepcopy_=True):
        if deepcopy_:
            atoms = [deepcopy(a) for a in atoms]
        return self.__class__(atoms)

    def _attr_set(self, attr):
        attr_set = set()
        attrs = []
        for atom in self.atoms:
            at = getattr(atom, attr)
            if at not in attr_set:
                attrs.append(at)
                attr_set.add(at)
        return attrs

    def _select_record(self, record):
        return [a for a in self.atoms if a.record == record]

    def protein(self):
        """Return a new object containing ATOM records.
        """
        atoms = self._select_record('ATOM')
        return self._init_with_atoms(atoms)

    def ligand(self):
        """Return a new object containing HETATM records.
        """
        atoms = self._select_record('HETATM')
        return self._init_with_atoms(atoms)

    def backbone(self):
        """Return a new object containing all backbone atoms
        """
        bbatm = {'N', 'HN', 'CA', 'HCA', 'C', 'O', 'OXT', 'HA', '1HA', '2HA',
                 'HT1', 'HT2', 'HT3'}
        atoms = [a for a in self.atoms if a.record == 'ATOM' and a.atom in bbatm]
        return self._init_with_atoms(atoms)

    def sidechains(self):
        """Return a new object containing all sidechain atoms.
        """
        bbatm = {'N', 'HN', 'CA', 'HCA', 'C', 'O', 'OXT', 'HA', '1HA', '2HA',
                 'HT1', 'HT2', 'HT3'}
        atoms = [a for a in self.atoms if (a.record == 'ATOM' and
                                           a.atom not in bbatm)]
        return self._init_with_atoms(atoms)

    @property
    def chains(self):
        """Return a list of the chains in this object.
        """
        return self._attr_set('chain')

    @property
    def resi(self):
        """Return a list of residue numbers.
        """
        return self._attr_set('nres')

    def select(self, **kwargs):
        """Select atoms based on logical criteria

        Examples:

        >>> prot = AtomCollection(atoms)
        >>> chain_a = prot.select(chain__eq='A')
        >>> other_chains = prot.select(chain__ne='A')
        >>> residues_1_to_250 = prot.select(nres__le=250)
        >>> hydrogens = prot.select(atom__startswith='H')
        >>> res356_A = prot.select(chain='A', nres=356)  # default is 'eq'
        >>> selection = prot.select(nres__gt=200, nres__lt=245).sidechains()

        """
        atoms = self.atoms[::]
        for key, value in kwargs.iteritems():
            try:
                attr, func = key.split('__')
            except ValueError, exception:
                if len(key.split('__')) == 1:
                    attr = key
                    func = 'eq'
                else:
                    raise ValueError(exception)
            if hasattr(operator, func):
                f = getattr(operator, func)
            elif hasattr(str, func):
                f = lambda obj, value: getattr(obj, func)(value)
            else:
                raise AttributeError("Can't apply {0}".format(func))
            atoms = [a for a in atoms if f(getattr(a, attr), value)]
        if atoms:
            return self._init_with_atoms(atoms)
        return None

    def renumber_residues(self, start=1):
        """Renumber residues beginning with start.
        """
        it = iter(self.atoms)
        atom = next(it)
        delta = atom.nres - start
        atom.nres -= delta
        for atom in it:
            atom.nres -= delta
            atom.connections = [at - delta for at in atom.connections]

    def remove_protons(self):
        """Remove protons from collection
        """
        atoms = [p for p in self.atoms if not p.atom.startswith('H')]
        obj = self._init_with_atoms(atoms)
        obj.renumber_atoms()
        return obj

    def renumber_atoms(self):
        """Renumber atoms
        """
        key = {}
        for i, atom in enumerate(self.atoms, start=1):
            key[atom.natom] = i
            atom.natom = i
        for atom in self.atoms:
            atom.connections = [
                key[num] for num in atom.connections if num in key
            ]

    @property
    def sequence(self):
        """Returns one-letter code sequence of collection.
        Missing residues are filled in with -'s
        """
        key = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'CYX': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HID': 'H',
            'HIE': 'H', 'HIP': 'H', 'HIS': 'H', 'HSE': 'H', 'HSP': 'H',
            'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
            'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y',
            'VAL': 'V'
        }
        chain_seq = []
        for chain in self.chains:
            seq = {
                a.nres: key[a.res] for a in self.atoms if (
                    a.record == 'ATOM' and a.atom == 'CA' and a.chain == chain
                )
            }
            MIN = min(seq)
            MAX = max(seq) + 1
            chain_seq.append(''.join(seq.get(i, '-') for i in xrange(MIN, MAX)))
        return ''.join(chain_seq)

    def rmsd(self, other):
        """Return the root-mean-square deviation (RMSD) between this and another
        object.
        """
        factor = np.sqrt(1 / np.float(len(self.atoms)))
        return factor * np.linalg.norm(self.matrix - other.matrix)

    def aligned_rmsd(self, other):
        """Return the RMSD of the aligned structures of this and another object.
        Both objects must have the same number of atoms.
        """
        mol1 = other.matrix
        mol2 = self.matrix
        return aligned_rmsd(mol1, mol2)

    def align(self, other):
        """Align this object to another and return a new object with the
        aligned coordinates.
        """
        mol1 = other.matrix
        mol2 = self.matrix
        new_coords = align(mol1, mol2)
        atoms = [deepcopy(a) for a in self.atoms]
        for at, new_coord in zip(atoms, new_coords):
            at.coord = new_coord
            at.x, at.y, at.z = new_coord
        return self._init_with_atoms(atoms, deepcopy_=False)
