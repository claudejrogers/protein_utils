from .atom import Atom, AtomCollection
from .pdb import PDBAtom, PDBFile


class PQRAtom(Atom):
    """Parse PQR atom format
    """
    def __init__(self, record, natom, atom, res, chain, nres, x, y, z, charge,
                 radius):
        super(self.__class__, self).__init__(record, natom, atom, res, chain,
                                             nres, x, y, z, charge)
        self.radius = radius

    @classmethod
    def from_line(cls, line):
        """PQR files are whitespace delimited
        """
        RECORD = slice(0, 6)
        NATOM = slice(6, 11)
        ATOM = slice(12, 16)
        RES = slice(17, 20)
        CHAIN = slice(21, 22)
        NRES = slice(22, 26)
        X = slice(30, 38)
        Y = slice(38, 46)
        Z = slice(46, 54)
        CHARGE = slice(54, 62)
        RADIUS = slice(62, None)
        record = line[RECORD].strip()
        natom = int(line[NATOM])
        atom = line[ATOM].strip()
        res = line[RES].strip()
        chain = line[CHAIN]
        nres = int(line[NRES])
        x = float(line[X])
        y = float(line[Y])
        z = float(line[Z])
        charge = float(line[CHARGE])
        radius = float(line[RADIUS])
        return cls(record, natom, atom, res, chain, nres, x, y, z, charge,
                   radius)

    def writeline(self):
        if len(self.atom) < 4:
            atom = ' {0}'.format(self.atom)
        else:
            atom = self.atom
        aline = (
            "{0:6s}{1:5d} {2:4s} {3:3s} {4:1s}{5:4d}    {6:8.3f}"
            "{7:8.3f}{8:8.3f} {9:7.4f}{10:7.4f}\n"
        ).format(
            self.record, self.natom, atom, self.res, self.chain,
            self.nres, self.x, self.y, self.z, self.charge,
            self.radius
        )
        return aline

    def __repr__(self):
        if len(self.atom) < 4:
            atom = ' {0}'.format(self.atom)
        else:
            atom = self.atom
        return "<PQR: {0:6s}{1:5d} {2:4s} {3:3s} {5:1s}{6:4d} ...>".format(
            self.record, self.natom, atom, self.res, self.chain, self.nres
        )


class PQRFile(AtomCollection):

    def __init__(self, ff, atoms=None):
        self.ff = ff
        super(self.__class__, self).__init__(atoms)

    @classmethod
    def from_file(cls, filename):
        ff = ''
        atoms = []
        with open(filename) as f:
            for line in f:
                if line.startswith('REMARK'):
                    if 'Forcefield Used:' in line:
                        ff = line.strip().split(': ')[-1]
                elif line.startswith(('ATOM', 'HETATM')):
                    p = PQRAtom.from_line(line)
                    atoms.append(p)
                else:
                    continue
        return cls(ff, atoms)

    def write_pqr(self, filename):
        body = ["REMARK   1 Forcefield Used: {0}\n".format(self.ff)]
        for a in self.atoms:
            line = a.writeline()
            body.append(line)
        body.append("TER\nEND\n")
        with open(filename, 'w') as f:
            f.writelines(body)

    def to_pdb(self):
        pdb_atoms = []
        for a in self.atoms:
            record = a.record
            natom = a.natom
            atom = a.atom
            altloc = ''
            res = a.res
            chain = a.chain
            nres = a.nres
            icode = ''
            x = a.x
            y = a.y
            z = a.z
            occupancy = 0.00
            bfactor = 0.00
            element = atom[0]
            charge = ''
            p = PDBAtom(record, natom, atom, altloc, res, chain, nres, icode,
                        x, y, z, occupancy, bfactor, element, charge)
            pdb_atoms.append(p)
        return PDBFile(pdb_atoms)

    def write_pdb(self, filename):
        pdb = self.to_pdb()
        pdb.write_pdb(filename)
