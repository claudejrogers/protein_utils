from itertools import groupby
from operator import attrgetter
from urllib2 import urlopen

from .residue import Residue, Residues
from .atom import Atom, AtomCollection


class PDBAtom(Atom):
    """PDB file class
    """
    def __init__(self, record, natom, atom, altloc, res, chain, nres, icode,
                 x, y, z, occupancy, bfactor, element, charge):
        super(self.__class__, self).__init__(record, natom, atom, res, chain,
                                             nres, x, y, z, charge)
        self.altloc = altloc
        self.icode = icode
        self.occupancy = occupancy
        self.bfactor = bfactor
        self.element = element

    @classmethod
    def from_line(cls, line):
        """Parse a PDB RECORD line.

        Parameters
        ----------
        line : string
            A RECORD line (ATOM/HETATM) from a PDB file

        Returns
        -------
        result : PDBAtom
        """
        RECORD = slice(0, 6)
        NATOM = slice(6, 11)
        ATOM = slice(12, 16)
        ALTLOC = slice(16, 17)
        RES = slice(17, 20)
        CHAIN = slice(21, 22)
        NRES = slice(22, 26)
        ICODE = slice(26, 27)
        X = slice(30, 38)
        Y = slice(38, 46)
        Z = slice(46, 54)
        OCCUP = slice(54, 60)
        BFACT = slice(60, 66)
        ELEM = slice(76, 78)
        CHARGE = slice(78, 80)
        record = line[RECORD].strip()
        try:
            natom = int(line[NATOM])
        except ValueError:
            natom = int(line[NATOM], base=16)
        atom = line[ATOM].strip()
        altloc = line[ALTLOC].strip()
        res = line[RES].strip()
        chain = line[CHAIN]
        nres = int(line[NRES])
        icode = line[ICODE].strip()
        x = float(line[X])
        y = float(line[Y])
        z = float(line[Z])
        occupancy = float(line[OCCUP])
        bfactor = float(line[BFACT])
        element = line[ELEM].strip() or atom[0]
        try:
            charge = line[CHARGE]
        except ValueError:
            charge = ''
        return cls(record, natom, atom, altloc, res, chain, nres, icode, x, y, z,
                   occupancy, bfactor, element, charge)

    def add_connections_from_line(self, line):
        if not line.startswith('CONECT') and int(line[6:11]) != self.natom:
            return
        line = line.strip()
        it = [iter(line[11:])] * 5
        self.connections.extend(int(''.join(num)) for num in zip(*it))

    def writeline(self):
        """Write object data as formated RECORD and CONECT strings

        Returns
        -------
        line : tuple
            A two-element tuple containing the RECORD (ATOM/HETATM) line and
            the CONECT line (if applicable).
        """
        if self.record == 'ATOM' and not self.atom.startswith('H') and len(self.atom) < 4:
            atom = ' {0}'.format(self.atom)
        else:
            atom = self.atom
        if self.natom > 99999:
            natom = hex(self.natom)[2:]
        else:
            natom = "{0:5d}".format(self.natom)
        aline = (
            "{0:6s}{1} {2:4s}{3:1s}{4:3s} {5:1s}{6:4d}{7:1s}   {8:8.3f}"
            "{9:8.3f}{10:8.3f}{11:6.2f}{12:6.2f}          {13:>2s}{14:>2s}\n"
        ).format(
            self.record, natom, atom, self.altloc, self.res, self.chain,
            self.nres, self.icode, self.x, self.y, self.z, self.occupancy,
            self.bfactor, self.element, self.charge
        )
        cline = ''
        if self.connections:
            cline = "CONECT{0:5d}{1}\n".format(
                self.natom,
                ''.join('{0:5d}'.format(num) for num in self.connections)
            )
        return aline, cline

    def __repr__(self):
        if self.record == 'ATOM' and not self.atom.startswith('H'):
            atom = ' {0}'.format(self.atom)
        else:
            atom = self.atom
        return "<PDB: {0:6s}{1:5d} {2:4s}{3:1s}{4:3s} {5:1s}{6:4d}{7} ...>".format(
            self.record, self.natom, atom, self.altloc, self.res, self.chain,
            self.nres, self.icode
        )


class PDBResidue(Residue):

    def __init__(self, res, nres, chain, icode, atoms):
        self.icode = icode.strip() or ''
        super(self.__class__, self).__init__(res, nres, chain, atoms)

    @classmethod
    def from_atoms(cls, atoms):
        """Instantiate new residue object from sorted atom list

        Parameters
        ----------
        atoms : list
            A list of PDBAtom objects sorted by natom. List must contain only
            one residue.

        Returns
        -------
        result : PDBResidue
        """
        g = groupby(atoms, attrgetter('nres', 'chain', 'icode', 'res'))
        resi = [(key, it) for key, it in g]
        if len(resi) != 1:
            raise ValueError(
                '{0} can only contain one residue.'.format(
                    cls.__name
                )
            )
        (nres, chain, icode, res), it = resi[0]
        return cls(res, nres, chain, icode, list(it))

    def get_next(self):
        """Returns key for 'next' residue
        """
        return '{0}{1}_{2}'.format(self.nres + 1, self.icode, self.chain)

    def get_prev(self):
        """Returns key for 'previous' residue
        """
        return '{0}{1}_{2}'.format(self.nres - 1, self.icode, self.chain)


class PDBResidues(Residues):

    def __init__(self, atoms):
        g = groupby(atoms, attrgetter('nres', 'chain', 'icode', 'res'))
        self._dict = {
            '{0}{1}_{2}'.format(nres, icode.strip(), chain): PDBResidue(
                res, nres, chain, icode, list(it)
            ) for (nres, chain, icode, res), it in g
        }
        self._list = [v for v in sorted(
            self._dict.itervalues(),
            key=attrgetter('nres', 'chain', 'icode')
        )]


def _read_file_iterator(iterable):
    atoms = {}
    for line in iterable:
        if line.startswith(('ATOM', 'HETATM')):
            p = PDBAtom.from_line(line)
            atoms[p.natom] = p
        elif line.startswith('CONECT'):
            natom = int(line[6:11])
            line = line.strip()
            atoms[natom].add_connections_from_line(line)
    return sorted(atoms.values(), key=attrgetter('natom'))


class PDBFile(AtomCollection):

    @classmethod
    def from_file(cls, filename):
        """Initialize from file

        Parameters
        ----------
        filename : path
            Path to a PDB file

        Returns
        -------
        result : PDBFile
        """
        atoms = {}
        with open(filename) as f:
            atoms = _read_file_iterator(f)
        return cls(atoms)

    @classmethod
    def fetch(cls, pdb_code):
        """Fetch pdb from RCSB Protein Data Bank.

        Parameters
        ----------
        pdb_code : str
            Four character PDB code.

        Returns
        -------
        PF : PDBFile
            New PDBFile object.
        """
        url = "http://www.rcsb.org/pdb/files/{0}.pdb".format(pdb_code)
        f = urlopen(url)
        atoms = _read_file_iterator(f)
        f.close()
        return cls(atoms)

    @classmethod
    def fetch_ligand(cls, cci):
        """Fetch ligand from Ligand Expo

        Parameters
        ----------
        cci : string
            Three-letter chemical component identifier

        Returns
        -------
        PF : PDFFile
            New PDBFile object
        """
        url = ('http://ligand-expo.rcsb.org/reports/'
               '{0[0]}/{0}/{0}_model.pdb'.format(cci))
        f = urlopen(url)
        atoms = _read_file_iterator(f)
        f.close()
        return cls(atoms)

    def ramachandran_plot(self):
        """Creates Ramachandran plot of phi and psi backbone angles.
        """
        residues = PDBResidues(self.atoms)
        residues.ramachandran_plot()

    def write_pdb(self, filename):
        """Write object to PDB file

        Parameters
        ----------
        filename : path
            Path to save PDB file.
        """
        atoms = []
        conect = []
        for atom in self.atoms:
            a, c = atom.writeline()
            atoms.append(a)
            if c:
                conect.append(c)
        with open(filename, 'w') as f:
            f.writelines(atoms)
            f.writelines(conect)
            f.write("END\n")
