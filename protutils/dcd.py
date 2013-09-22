"""Parse DCD trajectory files.

Based on http://www.csb.pitt.edu/prody/_modules/prody/trajectory/dcdfile.html
"""
from os.path import getsize
from struct import calcsize, unpack
import numpy as np

from .pdb import PDBFile
from .cealign.kabsch import align


class DCDFile(object):
    """Parse DCD trajectory files.
    """
    def __init__(self, filename):
        self._filename = filename
        self._nfi = 0
        self._file = open(filename, 'rb')
        dcd = self._file
        endian = ''
        charmm = None
        dcdcordmagic = unpack('i', 'CORD')[0]

        bits = dcd.read(calcsize('ii'))

        temp = unpack('ii', bits)

        if temp[0] + temp[1] == 84:
            raise IOError('Unsupported DCD format.')
        else:
            if unpack('>ii', bits) == temp:
                endian = '>'
            else:
                endian = '<'

            temp = unpack(endian + 'ii', bits)
            if temp[0] == 84 and temp[1] == dcdcordmagic:
                print "Standard CHARMM 32-bit DCD file"
            else:
                raise IOError('Unsupported DCD format.')

        bits = dcd.read(80)
        temp = unpack(endian + 'i' * 20, bits)

        if temp[-1] != 0:
            charmm = True

        if not charmm:
            raise NotImplementedError('This format is not yet supported.')
        temp = unpack(endian + 'i' * 9 + 'f' + 'i' * 10, bits)
        self._n_csets = temp[0]
        self._first_ts = temp[1]
        self._framefreq = temp[2]
        self._n_fixed = temp[8]
        if self._n_fixed > 0:
            raise NotImplementedError(
                'DCD files with fixed atoms are not yet supported.'
            )
        self._timestep = temp[9]
        self._unitcell = temp[10] == 1

        if unpack(endian + 'i', dcd.read(calcsize('i')))[0] != 84:
            raise IOError('Unrecognized DCD format.')

        temp = unpack(endian + 'i', dcd.read(calcsize('i')))

        if temp[0] != 164:
            raise IOError('Unrecognized DCD format.')

        temp = unpack(endian + 'i', dcd.read(calcsize('i')))

        self._dcdtitle = dcd.read(80)
        self._remarks = dcd.read(80)

        temp = unpack(endian + 'i', dcd.read(calcsize('i')))

        if temp[0] != 164:
            raise IOError('Unrecognized DCD format.')

        if unpack(endian + 'i', dcd.read(calcsize('i')))[0] != 4:
            raise IOError('Unrecognized DCD format.')

        self._n_atoms = unpack(endian + 'i', dcd.read(calcsize('i')))[0]

        if unpack(endian + 'i', dcd.read(calcsize('i')))[0] != 4:
            raise IOError('Bad DCD format.')

        self._endian = endian
        self._n_floats = (self._n_atoms + 2) * 3

        bytes_per_frame = self._n_floats * 4

        if self._unitcell:
            bytes_per_frame += 56
        self._dtype = np.float32
        self._itemsize = 4
        self._bytes_per_frame = bytes_per_frame
        self._first_byte = self._file.tell()
        n_csets = ((getsize(self._filename) - self._first_byte) /
                   self._bytes_per_frame)
        if self._n_csets != n_csets:
            print """\
WARNING: DCD header claims {0} frames, file size indicates
there are actually {1} frames\
            """.format(self._n_csets, n_csets)
            self._n_csets = n_csets

        self._coords = self.next_coordset()
        self._file.seek(self._first_byte)
        self._nfi = 0

    def __next__(self):
        return self.next_coordset()

    def __iter__(self):
        while self._nfi < self._n_csets:
            yield self.__next__()

    def __getitem__(self, index):
        return self.get_coordset(index)

    def next_coordset(self):
        if self._file.closed:
            raise ValueError('IO operation on closed file')
        if self._nfi >= self._n_csets:
            raise IndexError('No more frames')
        if self._unitcell:
            self._file.seek(56, 1)
        return self._next_coordset()

    def _next_coordset(self):
        n_floats = self._n_floats
        n_atoms = self._n_atoms
        coords = np.fromstring(self._file.read(
            self._itemsize * n_floats
        ), self._dtype)
        if len(coords) != n_floats:
            return None
        coords = coords.reshape((3, n_atoms + 2)).T[1:-1, :]
        coords = coords.reshape((n_atoms, 3))
        self._nfi += 1
        return coords

    def get_coordset(self, index):
        """Get the coordinates for a frame at the specified index.

        Parameters
        ----------
        index : int
            The desired frame of dynamics
        """
        nfi = self._nfi
        self.goto(index)
        coords = self.next_coordset()
        self.goto(nfi)
        return coords

    def goto(self, n):
        """Advances the DCD file to a specified frame

        Parameters
        ----------
        n : int
            The desired frame of dynamics
        """
        if self._file.closed:
            raise ValueError('IO operation on closed file')
        n_csets = self._n_csets
        # support negative indexes
        if n < 0:
            n = n_csets + n
            if n < 0:
                n = 0
        elif n > n_csets:
            raise ValueError('File contains {0} frames'.format(n_csets))
        self._file.seek(self._first_byte + n * self._bytes_per_frame)
        self._nfi = n

    def avg_structure(self, pdbfile, selector=None):
        """Get structure closest to the average RMSD for the trajectory.

        Parameters
        ----------
        pdbfile : path
            Path to the input PDB structure for dynamics

        selector : {'CA', None} (optional, default=None)
            Specifies which atoms to consider when determining the average
            structure.
            'CA' considers only the alpha carbons
            None (default) considers any atoms with an occupancy value greater
            than 0.0.

        Returns
        -------
        result : PDBFile
            Creates a PDBFile object with the coordinates of the frame closest
            to the average pose during dynamics.
        """
        pdb = PDBFile.from_file(pdbfile)
        if len(pdb) != self._n_atoms:
            raise ValueError('Number of atoms in the PDB file does not match'
                             ' number of atoms in DCD file.')
        if selector == 'CA':
            indices = [i for i, a in enumerate(pdb.atoms) if a.atom == 'CA']
        else:
            indices = [i for i, a in enumerate(pdb.atoms) if a.occupancy > 0.0]
        self.goto(0)
        top = self.__next__()[indices]
        coords = [top]
        print "Aligning coordinates..."
        for c in self.__iter__():
            xyz = align(top, c[indices])
            coords.append(xyz)
        coords = np.array(coords)
        avg = coords.mean(axis=0)
        rmsds = []
        print "Calculating RMSDs"
        for c in coords:
            rmsd = np.linalg.norm(c - avg)
            rmsds.append(rmsd)
        rmsds = np.array(rmsds) / np.sqrt(avg.shape[0])
        lowest = rmsd.argmin()
        print "Extracting frame {0}".format(lowest)
        self.goto(lowest)
        new_coords = self.__next__()
        return pdb.new_coordinates(new_coords)

    def close(self):
        """Close DCD file object.
        """
        self._file.close()
