import operator
from copy import deepcopy
from itertools import groupby, product
import numpy as np

import residue
from .cealign.kabsch import aligned_rmsd, align, align_substructure
from .cealign.cealign import distance_matrix, similarity_matrix, find_path


class Atom(object):
    """Base class for an atom description.

    Parameters
    ----------
    record : string {ATOM, HETATM}
        The type of atom. ATOM refers to atoms in amino acids, HETATM refers
        to atoms in other molecules.

    natom : int
        The atom number.

    atom : string
        The identifier for this atom.

    res : string
        The residue of this atom.

    chain : string
        The chain identifier for this atom.

    nres : int
        The residue number of this atom.

    x : float
        The value for the x coordinate of the atom.

    y : float
        The value for the y coordinate of the atom.

    z : float
        The value for the z coordinate of the atom.

    charge : float
        Charge on this atom.
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
        """Add a list of atom numbers (natom) that connect with this Atom

        Parameters
        ----------
        connections : list
            A list of natoms (ints) that connect to this Atom
        """
        self.connections.extend(connections)

    def distance(self, other):
        """Compute the distance between two Atoms
        """
        return np.linalg.norm(self.coord - other.coord)

    @staticmethod
    def dist(atm1, atm2):
        """Distance between two atoms.

        Parameters
        ----------
        atm1 : Atom
            First atom

        atm2 : Atom
            Second atom

        Returns
        -------
        dist : float
            The distance between atm1 and atm2
        """
        return np.linalg.norm(atm1.coord - atm2.coord)

    @staticmethod
    def dot(atm1, atm2):
        """Find dot product of two Atoms.

        Parameters
        ----------
        atm1 : Atom
            First atom

        atm2 : Atom
            Second atom

        Returns
        -------
        dot product : float
            The dot product,
            :math:`\\mathbf{A1}_{xyz} \\cdot \\mathbf{A2}_{xyz}`
        """
        return np.dot(atm1.coord, atm2.coord)

    @staticmethod
    def cross(atm1, atm2):
        """Find cross product of two Atom coordinates.

        Parameters
        ----------
        atm1 : Atom
            First atom

        atm2 : Atom
            Second atom

        Returns
        -------
        cross product : float
            The cross product,
            :math:`\\mathbf{A1}_{xyz} \\times \\mathbf{A2}_{xyz}`
        """
        return np.cross(atm1.coord, atm2.coord)

    @staticmethod
    def nomalize(atom):
        """Normalize Atom coordinates

        Parameters
        ----------
        atom : Atom
            First atom

        Returns
        -------
        coords : ndarray (1 x 3)
            Numpy array of normalized coordinates (|| coords || == 1)
        """
        return _normalize(atom.coord)

    @staticmethod
    def angle(atm1, atm2, atm3):
        """Find angle between three atoms

        Parameters
        ----------
        atm1 : Atom
            First atom

        atm2 : Atom
            Second atom (vertex)

        atm3 : Atom
            Third angle

        Returns
        -------
        angle : float
            Angle between the atoms
        """
        c1 = atm3.coord - atm2.coord
        c2 = atm1.coord - atm2.coord
        norm1 = _normalize(c1)
        norm2 = _normalize(c2)
        dotted = np.dot(norm1, norm2)
        if dotted > 1.0:
            dotted = 1.0
        rad = abs(np.arccos(dotted))
        angle = rad * 180.0 / np.pi
        if angle > 180.0:
            angle = 360.0 - angle
        return angle

    @staticmethod
    def dihedral(atm1, atm2, atm3, atm4):
        """Calculate dihedral angle

        Parameters
        ----------
        atm1 : Atom
            First atom

        atm2 : Atom
            Second atom

        atm3 : Atom
            Third atom

        atm4 : Atom
            Fourth atom

        Returns
        -------
        dihdral angle : float
            The dihedral angle in degrees.
        """
        delta43 = atm4.coord - atm3.coord
        delta32 = atm3.coord - atm2.coord
        delta12 = atm1.coord - atm2.coord

        A = np.cross(delta12, delta32)
        normA = _normalize(A)
        B = np.cross(delta43, delta32)
        normB = _normalize(B)

        scal = np.dot(normA, normB)
        if abs(scal + 1.0) < 1.0e-7:
            value = 180.0
        elif abs(scal - 1.0) < 1.0e-7:
            value = 0.0
        else:
            value = 57.2958 * np.arccos(scal)

        chiral = np.dot(np.cross(normA, normB), delta32)
        if chiral < 0:
            value *= -1.0
        return value

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
    """Base class for an Atom collection
    """
    def __init__(self, atoms=None):
        self.atoms = atoms or []
        self.matrix = np.array([atom.coord for atom in self.atoms])

    def __len__(self):
        """The number of atoms in this collection.
        """
        return len(self.atoms)

    def __iter__(self):
        """Iterate through the atoms in this collection
        """
        return iter(self.atoms)

    def __getitem__(self, index):
        """Get Atom for a given index

        Parameters
        ----------
        index : int
            The index of the Atom

        Returns
        -------
        atom : Atom
            The atom at the supplied index
        """
        return self.atoms[index]

    def __eq__(self, other):
        """Compares two AtomCollections for equality.

        Parameters
        ----------
        other : AtomCollection

        Returns
        -------
        result : boolean
        """
        if len(self) != len(other):
            return False
        return all(operator.eq(s, o) for s, o in zip(
            self.atoms, other.atoms
        ))

    def __contains__(self, other):
        """Returns True if all atoms in other are in this collection

        Parameters
        ----------
        other : AtomCollection

        Returns
        -------
        result : boolean
        """
        if len(self) < len(other):
            return False
        atoms_ = set(self.atoms)
        for a in other.atoms:
            if not a in atoms_:
                return False
        return True

    def __and__(self, other):
        """Returns the common atoms between self and other

        Parameters
        ----------
        other : AtomCollection

        Returns
        -------
        result : AtomCollection
            A new object containing atoms common in this and the other
            collection.

        Examples
        --------
        >>> a = AtomCollection(atoms_1)
        >>> b = AtomCollection(atoms_2)
        >>> common_atoms = a & b

        """
        atom_set = set(self.atoms) & set(other.atoms)
        atoms = sorted(atom_set, key=operator.attrgetter('natom'))
        return self._init_with_atoms(atoms)

    def __or__(self, other):
        """Returns the intersection of atoms in self and other

        Parameters
        ----------
        other : AtomCollection

        Returns
        -------
        result : AtomCollection
            A new object containing the atoms in this and the other
            collection.

        Examples
        --------
        >>> a = AtomCollection(atoms_1)
        >>> b = AtomCollection(atoms_2)
        >>> all_atoms = a | b

        """
        atom_set = set(self.atoms) | set(other.atoms)
        atoms = sorted(atom_set, key=operator.attrgetter('natom'))
        return self._init_with_atoms(atoms)

    def __xor__(self, other):
        """Returns atoms in self that are not in other and the atoms in
        other not in self.

        Parameters
        ----------
        other : AtomCollection

        Returns
        -------
        result : AtomCollection
            A new object that contains atoms in this collection that are not
            in the other collection.

        Examples
        --------
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

    def __repr__(self):
        return '<{0}: ({1} atoms)>'.format(self.__class__.__name__,
                                           len(self.atoms))

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

    def groupby_residue(self):
        """Return dict of residue records grouped by 'nres_chain' keys
        """
        g = groupby(self.atoms, operator.attrgetter('nres', 'chain'))
        return {
            "{0}_{1}".format(*key): self._init_with_atoms(it) for key, it in g
        }

    def protein(self):
        """Return a new object containing ATOM records.

        Returns
        -------
        results : AtomCollection
            A new object containing only protein atoms.

        See Also
        --------
        AtomCollection.ligand, AtomCollection.backbone,
        AtomCollection.sidechains
        """
        atoms = self._select_record('ATOM')
        return self._init_with_atoms(atoms)

    def ligand(self):
        """Return a new object containing HETATM records.

        Returns
        -------
        result : AtomCollection
            A new object containing non-protein atoms.

        See Also
        --------
        AtomCollection.protein, AtomCollection.backbone,
        AtomCollection.sidechains
        """
        atoms = self._select_record('HETATM')
        atoms = [a for a in atoms if a.res != 'HOH']  # no water
        return self._init_with_atoms(atoms)

    def backbone(self):
        """Return a new object containing all backbone atoms

        Returns
        -------
        result : AtomCollection
            A new object containing only backbone atoms

        See Also
        --------
        AtomCollection.protein, AtomCollection.ligands,
        AtomCollection.sidechains
        """
        bbatm = {'N', 'HN', 'CA', 'HCA', 'C', 'O', 'OXT', 'HA', '1HA', '2HA',
                 'HT1', 'HT2', 'HT3'}
        atoms = [a for a in self.atoms if a.record == 'ATOM' and a.atom in bbatm]
        return self._init_with_atoms(atoms)

    def sidechains(self):
        """Return a new object containing all sidechain atoms.

        Returns
        -------
        result : AtomCollection
            A new object containing only sidechain atoms.

        See Also
        --------
        AtomCollection.protein, AtomCollection.ligands,
        AtomCollection.backbone
        """
        bbatm = {'N', 'HN', 'CA', 'HCA', 'C', 'O', 'OXT', 'HA', '1HA', '2HA',
                 'HT1', 'HT2', 'HT3'}
        atoms = [a for a in self.atoms if (a.record == 'ATOM' and
                                           a.atom not in bbatm)]
        return self._init_with_atoms(atoms)

    @property
    def chains(self):
        """Return a list of the chains in this object.

        Returns
        -------
        chains : lists
            A list of the unique chain identifiers in this collection.
        """
        return self._attr_set('chain')

    @property
    def resi(self):
        """Return a list of residue numbers.
        """
        g = groupby(self.atoms, operator.attrgetter('nres', 'chain'))
        return ["{0}_{1}".format(*key) for key, it in g]

    @property
    def residues(self):
        """Returns a list of residues.
        """
        g = groupby(self.atoms, operator.attrgetter('nres', 'chain', 'res'))
        return [key[2] for key, it in g]

    def ramachandran_plot(self):
        """Creates Ramachandran plot for the protein atoms in this collection.
        """
        residues = residue.Residues(self.atoms)
        residues.ramachandran_plot()

    def _select_or_exclude(self, metric, **kwargs):
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
            if metric == 'select':
                atoms = [a for a in atoms if f(getattr(a, attr), value)]
            elif metric == 'exclude':
                atoms = [a for a in atoms if not f(getattr(a, attr), value)]
            else:
                raise ValueError("Unknown metric {0}.".format(metric))
        if atoms:
            return self._init_with_atoms(atoms)
        return None

    def select(self, **kwargs):
        """Select atoms based on logical criteria

        Parameters
        ----------
        kwargs : various
            Keys are composed of an atom attribute (chain, atom, natom, nres,
            record, etc.) and a boolean operator (contains, eq, le, gt,
            startswith), separated by two underscores.
            For example:
            chain__eq, nres__gt, atom__contains.

        Returns
        -------
        AC : AtomCollection
            New object containing the selected subset of Atoms.

        See Also
        --------
        AtomCollection.exclude

        Examples
        --------
        >>> prot = AtomCollection(atoms)
        >>> chain_a = prot.select(chain__eq='A')
        >>> other_chains = prot.select(chain__ne='A')
        >>> residues_1_to_250 = prot.select(nres__le=250)
        >>> hydrogens = prot.select(atom__startswith='H')
        >>> res356_A = prot.select(chain='A', nres=356)  # default is 'eq'
        >>> selection = prot.select(nres__gt=200, nres__lt=245).sidechains()

        """
        return self._select_or_exclude('select', **kwargs)

    def exclude(self, **kwargs):
        """Create new AtomCollection where specified atoms are excluded

        Parameters
        ----------
        kwargs
            Refer to AtomCollection.select documentation

        Returns
        -------
        results : AtomCollection
            A new object containing atoms that are not specified by the
            selection criteria.

        See Also
        --------
        AtomCollection.select
        """
        return self._select_or_exclude('exclude', **kwargs)

    def within(self, distance, other):
        """Get atoms in this instance within some distance of atoms in
        other instance.

        Parameters
        ----------
        distance : float
            The maximum distance between atoms in this object and the other
            object.

        other : AtomCollection
            The other object...

        Returns
        -------
        result : AtomCollection
            New object containing Atoms in this object that are less than
            the selected distance from the other object.
        """
        atoms = {s for s, o in product(
            self.atoms, other.atoms
        ) if s.distance(o) < distance}
        atoms = sorted(atoms, key=operator.attrgetter('natom'))
        return self._init_with_atoms(atoms)

    def renumber_residues(self, start=1):
        """Renumber residues beginning with start.

        Parameters
        ----------
        start : int, optional (default=1)
            The new number of the first residue
        """
        it = iter(self.atoms)
        atom = next(it)
        delta = atom.nres - start
        atom.nres -= delta
        for atom in it:
            atom.nres -= delta
            atom.connections = [at - delta for at in atom.connections]

    def remove_protons(self):
        """Remove protons from collection.
        """
        atoms = [p for p in self.atoms if not p.atom.startswith('H')]
        obj = self._init_with_atoms(atoms)
        obj.renumber_atoms()
        return obj

    def renumber_atoms(self, start=1):
        """Renumber atoms.

        Parameters
        ----------
        start : int, optional (default=1)
            The new number of the first atom.
        """
        key = {}
        for i, atom in enumerate(self.atoms, start=start):
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

        Returns
        -------
        sequence : string
            The sequence of the protein residues in this AtomCollection.
            Gaps in the sequence are indicated by '-'.
        """
        key = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'APP': 'D',
            'CYS': 'C', 'CYX': 'C', 'GLU': 'E', 'GLP': 'E', 'GLN': 'Q',
            'GLY': 'G', 'HID': 'H', 'HIE': 'H', 'HIP': 'H', 'HIS': 'H',
            'HSE': 'H', 'HSP': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
            'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T',
            'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
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

        Parameters
        ----------
        other : AtomCollection
            An atom colection containing the same number of atoms as this
            object.

        Returns
        -------
        rmsd : float
            The root-mean-square deviation between the two AtomCollections.

        See Also
        --------
        AtomCollection.aligned_rmsd
        """
        factor = np.sqrt(1 / np.float(len(self.atoms)))
        return factor * np.linalg.norm(self.matrix - other.matrix)

    def aligned_rmsd(self, other):
        """Return the RMSD of the aligned structures of this and another object.
        Both objects must have the same number of atoms.

        Parameters
        ----------
        other : AtomCollection
            An atom colection containing the same number of atoms as this
            object.

        Returns
        -------
        rmsd : float
            The root-mean-square deviation between the two AtomCollections
            after alignment.

        See Also
        --------
        AtomCollection.rmsd
        """
        mol1 = other.matrix
        mol2 = self.matrix
        return aligned_rmsd(mol1, mol2)

    def new_coordinates(self, matrix):
        """Create new object with new coordinates

        Parameters
        ----------
        matrix : numpy N x 3 matrix
            A matrix with the same shape as self.matrix

        Returns
        -------
        result : AtomCollection
            A new object with the supplied coordinates
        """
        if self.matrix != matrix.shape:
            raise ValueError("New coordinates must have the same dimensions"
                             " as current coordinates.")
        atoms = [deepcopy(a) for a in self.atoms]
        for at, new_coord in zip(atoms, matrix):
            at.coord = new_coord
            at.x, at.y, at.z = new_coord
        return self._init_with_atoms(atoms, deepcopy_=False)

    def orient(self):
        """Orient the protein to be centered at the origin with the
        principle axes aligned properly.

        Returns
        -------
        result : AtomCollection
            A new object with coordinates centered at the origin, with
            principle axes aligned with the Cartesian axes.
        """
        X = self.matrix
        COM = X.mean(axis=0)
        M = X - COM
        V, s, W = np.linalg.svd(M)
        xidx = s.argmax()
        zidx = s.argmin()
        yidx = list(set(range(3)) ^ set([xidx, zidx]))[0]

        A = np.zeros((3, 3))
        A[:, xidx] = [1.0, 0.0, 0.0]
        A[:, yidx] = [0.0, 1.0, 0.0]
        A[:, zidx] = [0.0, 0.0, 1.0]

        V, S, Wt = np.linalg.svd(np.dot(W.T, A))

        reflect = np.linalg.det(V) * np.linalg.det(Wt)

        if np.isclose(reflect, -1.0):
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]

        U = np.dot(V, Wt)

        new_coords = np.dot(M, U)

        return self.new_coordinates(new_coords)

    def align(self, other):
        """Align this object to another and return a new object with the
        aligned coordinates.

        Parameters
        ----------
        other : AtomCollection
            An atom colection containing the same number of atoms as this
            object.

        Returns
        -------
        results : AtomCollection
            New object with coordinates aligned to the other object.

        See Also
        --------
        AtomCollection.cealign
        """
        mol1 = other.matrix
        mol2 = self.matrix
        if mol1.shape != mol2.shape:
            raise ValueError("Objects must have the same number of atoms.")
        new_coords = align(mol1, mol2)
        return self.new_coordinates(new_coords)

    def cealign(self, other, alphas=True):
        """Align using the cealign algorithm from pymol.

        Parameters
        ----------
        other : AtomCollection
            An atom colection to align this object. The number of atoms in the
            two objects can be different.

        alphas : Boolean, optional (default=True)
            If True, align CA carbons only. If False, the full selection is
            aligned.

        Returns
        -------
        result : AtomCollection
            New object with coordinates aligned to the other object.

        See Also
        --------
        AtomCollection.align
        """
        if alphas:
            A = other.select(atom='CA').matrix
            B = self.select(atom='CA').matrix
        else:
            A = other.matrix
            B = self.matrix
        dA = distance_matrix(A)
        dB = distance_matrix(B)
        S = similarity_matrix(dA, dB)
        paths = find_path(S, dA, dB)
        scores = np.zeros(paths.shape[0])
        for i, p in enumerate(paths):
            m1 = A[p[:, 0]]
            m2 = B[p[:, 1]]
            scores[i] = aligned_rmsd(m1, m2)
        best_path = scores.argmin()

        A_MASK = paths[best_path, :, 0]
        B_MASK = paths[best_path, :, 1]

        MOL1 = A[A_MASK]
        MOL2 = B[B_MASK]

        U, COM1, COM2 = align_substructure(MOL1, MOL2)
        new_coords = np.dot(self.matrix - COM2, U) + COM1

        # Make new object
        return self.new_coordinates(new_coords)


def _normalize(coords):
    """Normalize coordinates

    coords must be a numpy array
    """
    dist = np.linalg.norm(coords)
    if dist > 1.0e-7:
        return coords / dist
    return coords
