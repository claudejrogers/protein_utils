from operator import attrgetter
from urllib2 import urlopen

from .atom import Atom, AtomCollection


class PDBLine(Atom):
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
        natom = int(line[NATOM])
        atom = line[ATOM].strip()
        altloc = line[ALTLOC].strip()
        res = line[RES].strip()
        chain = line[CHAIN]
        nres = int(line[NRES])
        icode = line[ICODE]
        x = float(line[X])
        y = float(line[Y])
        z = float(line[Z])
        occupancy = float(line[OCCUP])
        bfactor = float(line[BFACT])
        element = line[ELEM].strip()
        try:
            charge = line[CHARGE]
        except ValueError:
            charge = ''
        return cls(record, natom, atom, altloc, res, chain, nres, icode, x, y, z,
                   occupancy, bfactor, element, charge)

    # def add_connections(self, *connections):
    #     self.conections.extend(connections)

    def add_connections_from_line(self, line):
        if not line.startswith('CONECT') and int(line[6:11]) != self.natom:
            return
        it = [iter(line[11:])] * 5
        self.connections.extend(int(''.join(num)) for num in zip(*it))

    # def distance(self, other):
    #     return np.linalg.norm(self.coord - other.coord)

    def writeline(self):
        if self.record == 'ATOM' and not self.atom.startswith('H'):
            atom = ' {0}'.format(self.atom)
        else:
            atom = self.atom
        aline = (
            "{0:6s}{1:5d} {2:4s}{3:1s}{4:3s} {5:1s}{6:4d}{7:1s}   {8:8.3f}"
            "{9:8.3f}{10:8.3f}{11:6.2f}{12:6.2f}          {13:>2s}{14:>2s}\n"
        ).format(
            self.record, self.natom, atom, self.altloc, self.res, self.chain,
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
        return "<PDB: {0:6s}{1:5d} {2:4s}{3:1s}{4:3s} {5:1s}{6:4d} ...>".format(
            self.record, self.natom, atom, self.altloc, self.res, self.chain,
            self.nres
        )


def _read_file_iterator(iterable):
    atoms = {}
    for line in iterable:
        if line.startswith(('ATOM', 'HETATM')):
            p = PDBLine.from_line(line)
            atoms[p.natom] = p
        elif line.startswith('CONECT'):
            natom = int(line[6:11])
            line = line.strip()
            atoms[natom].add_connections_from_line(line)
    return sorted(atoms.values(), key=attrgetter('natom'))


class PDBFile(AtomCollection):

    @classmethod
    def from_file(cls, filename):
        """
        Initialize from file
        """
        atoms = {}
        with open(filename) as f:
            atoms = _read_file_iterator(f)
        return cls(atoms)

    @classmethod
    def fetch(cls, pdb_code):
        """
        Fetch pdb from RCSB Protein Data Bank
        """
        url = "http://www.rcsb.org/pdb/files/{0}.pdb".format(pdb_code)
        f = urlopen(url)
        atoms = _read_file_iterator(f)
        f.close()
        return cls(atoms)

    def write_pdb(self, filename):
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
