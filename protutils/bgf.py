from copy import deepcopy
from operator import attrgetter

from .atom import Atom, AtomCollection
from .pdb import PDBAtom, PDBFile


class BGFAtom(Atom):
    """Parse, store, and manipulate RECORD line data from a BGF file.
    """
    def __init__(self, record, natom, atom, res, chain, nres, x, y, z, fftype,
                 nbond, nlonepair, charge, fixed):
        super(self.__class__, self).__init__(record, natom, atom, res, chain,
                                             nres, x, y, z, charge)
        self.fftype = fftype
        self.nbond = nbond
        self.nlonepair = nlonepair
        self.charge = charge
        self.fixed = fixed

    @classmethod
    def from_line(cls, line):
        """Parse a RECORD line from a BGF file.

        Parameters
        ----------
        line : string
            BGF-formatted line

        Returns
        -------
        result : BGFAtom
            Instantiates a BGFAtom object
        """
        # Define slices
        RECORD = slice(0, 6)
        NATOM = slice(7, 12)
        ATOM = slice(13, 18)
        RES = slice(19, 22)
        CHAIN = slice(23, 24)
        NRES = slice(24, 29)
        X = slice(30, 40)
        Y = slice(40, 50)
        Z = slice(50, 60)
        TYPE = slice(61, 66)
        NBOND = slice(66, 69)
        NLP = slice(70, 71)
        CHARGE = slice(71, 80)
        FIXED = slice(81, 82)
        record = line[RECORD].strip()
        natom = int(line[NATOM])
        atom = line[ATOM].strip()
        res = line[RES].strip()
        chain = line[CHAIN].strip()
        nres = int(line[NRES])
        x = float(line[X])
        y = float(line[Y])
        z = float(line[Z])
        fftype = line[TYPE].strip()
        nbond = int(line[NBOND])
        nlonepair = int(line[NLP])
        charge = float(line[CHARGE])
        try:
            fixed = int(line[FIXED])
        except IndexError:
            fixed = 0
        return cls(record, natom, atom, res, chain, nres, x, y, z, fftype,
                   nbond, nlonepair, charge, fixed)

    def add_connections_from_line(self, line):
        """Add connection data
        """
        if not line.startswith('CONECT') or int(line[6:12]) != self.natom:
            return
        it = [iter(line[12:])] * 6  # group conect records
        self.connections.extend(int(''.join(num)) for num in zip(*it))

    def writeline(self):
        """Write BGF ATOM/HETATM line and CONECT info.

        Returns
        -------
        result : tuple
            A two-element tuple containing the ATOM/HETATOM record, and the
            CONECT information.
        """
        if self.record == 'ATOM' and not self.atom.startswith('H'):
            atom = ' {0}'.format(self.atom)
        else:
            atom = self.atom
        aline = (
            "{0:6s} {1:5d} {2:5s} {3} {4}{5:5d} {6:10.5f}{7:10.5f}"
            "{8:10.5f} {9:5s}{10:3d}{11:2d} {12:8.5f} {13}   0\n"
        ).format(
            self.record, self.natom, atom, self.res, self.chain, self.nres,
            self.x, self.y, self.z, self.fftype, self.nbond, self.nlonepair,
            self.charge, self.fixed
        )
        cline = "CONECT{0:6d}{1}\n".format(
            self.natom,
            ''.join('{0:6d}'.format(num) for num in self.connections)
        )
        return aline, cline

    def __repr__(self):
        if self.record == 'ATOM' and not self.atom.startswith('H'):
            atom = ' {0}'.format(self.atom)
        else:
            atom = self.atom
        return "<BGF: {0:6s} {1:5d} {2:5s} {3} {4}{5:5d} ...>".format(
            self.record, self.natom, atom, self.res, self.chain, self.nres
        )


class BGFFile(AtomCollection):
    """Store and manipulate BGF files.
    """
    def __init__(self, biogrf='332', descrp='', ff='DREIDING', atoms=None):
        """Initialize a BGFFile object.

        Creates a blank BGF object by default.

        Parameters
        ----------
        biogrf : str, optional (default='332')
            Biogrf format

        descrp : str, optional
            BGF description.

        ff : str, optional (default='DREIDING')
            Forcefield used to calculate atom charges and fftype

        atoms : list, optional
            List of BGFAtom objects used to create BGFFile instance

        Returns
        -------
        BF : BGFFile
        """
        self.biogrf = biogrf
        self.descrp = ''
        self.ff = ff
        super(self.__class__, self).__init__(atoms)

    @classmethod
    def from_file(cls, filename):
        """Create a new BGFFile object from a file.

        Parses BGF file header and RECORD lines.

        Parameter
        ---------
        filename : Path to BGF file

        Returns
        -------
        BF : BGFFile
            New BGFFile object containing data from the file.
        """
        biogrf = ''
        ff = ''
        descrp = ''
        atoms = {}
        with open(filename) as f:
            for line in f:
                if line.startswith('BIOGRF'):
                    biogrf = line.strip().split()[1]
                elif line.startswith('DESCRP'):
                    descrp = line.strip().split()[1]
                elif line.startswith('FORCEFIELD'):
                    ff = line.strip().split()[1]
                elif line.startswith(('ATOM', 'HETATM')):
                    b = BGFAtom.from_line(line)
                    atoms[b.natom] = b
                elif line.startswith('CONECT'):
                    natom = int(line[6:12])
                    atoms[natom].add_connections_from_line(line)
                else:
                    continue
        atoms = sorted(atoms.values(), key=attrgetter('natom'))
        return cls(biogrf, descrp, ff, atoms)

    def _init_with_atoms(self, atoms, deepcopy_=True):
        if not all(isinstance(b, BGFAtom) for b in atoms):
            raise ValueError("List objects must be instances of BGFAtom.")
        if deepcopy_:
            atoms = [deepcopy(a) for a in atoms]
        return self.__class__(self.biogrf, self.descrp, self.ff, atoms)

    @property
    def fftypes(self):
        return self._attr_set('fftype')

    def set_movable(self, movable=True):
        m = 0 if movable else 1
        for atom in self.atoms:
            atom.fixed = m

    def write_bgf(self, filename):
        """Write object to file

        Parameters
        ----------
        filename : path
            Path to desired file
        """
        body = ["BIOGRF{0:>5s}\n".format(self.biogrf)]
        if self.descrp:
            body.append("DESCRP {0}\n".format(self.descrp))
        body.append("FORCEFIELD {0}\n".format(self.ff))
        body.append("FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5"
                    ",1x,a5,i3,i2,1x,f8.5,i2,i4,f10.5)\n")
        atoms = []
        hetatms = []
        conect = []
        for atom in self.atoms:
            a, c = atom.writeline()
            if atom.record == 'ATOM':
                atoms.append(a)
            elif atom.record == 'HETATM':
                hetatms.append(a)
            conect.append(c)
        body.extend(atoms)
        body.extend(hetatms)
        body.append("FORMAT CONECT (a6,14i6)\nFORMAT ORDER (a6,i6,13f6.3)\n")
        body.extend(conect)
        body.append("END\n")
        with open(filename, 'w') as f:
            f.writelines(body)

    def to_pdb(self):
        """Create PDBFile object.

        Returns
        -------
        result : PDBFile
        """
        pdb_atoms = []
        for _atom in self.atoms:
            record = _atom.record
            natom = _atom.natom
            atom = _atom.atom
            altloc = ''
            res = _atom.res
            chain = _atom.chain
            nres = _atom.nres
            icode = ''
            x = _atom.x
            y = _atom.y
            z = _atom.z
            occupancy = 0.00
            bfactor = 0.00
            element = atom[0]
            charge = ''
            p = PDBAtom(record, natom, atom, altloc, res, chain, nres, icode,
                        x, y, z, occupancy, bfactor, element, charge)
            pdb_atoms.append(p)
        return PDBFile(pdb_atoms)

    def write_pdb(self, filename):
        """Write object contents as a PDB file.
        """
        pdb = self.to_pdb()
        pdb.write_pdb(filename)
