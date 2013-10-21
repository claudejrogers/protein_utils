import os
from glob import glob
from itertools import product

import numpy as np
try:
    from numpy import isclose
except ImportError:
    from ._operator import isclose
from sklearn.cluster import KMeans


class Sphere(object):
    """Parse and manipulate sph records
    """
    def __init__(self, iatom, x, y, z, radius, jatom, cluster, color):
        self.iatom = iatom
        self.x = x
        self.y = y
        self.z = z
        self.radius = radius
        self.jatom = jatom
        self.cluster = cluster
        self.color = color
        self.coord = np.array([x, y, z])

    @classmethod
    def from_line(cls, line):
        """Initialize Sphere object from line
        """
        IATOM = slice(0, 5)
        X = slice(5, 15)
        Y = slice(15, 25)
        Z = slice(25, 35)
        RAD = slice(35, 43)
        JATOM = slice(43, 48)
        CLUST = slice(48, 50)
        COLOR = slice(50, 53)
        iatom = int(line[IATOM])
        x = float(line[X])
        y = float(line[Y])
        z = float(line[Z])
        radius = float(line[RAD])
        jatom = int(line[JATOM])
        cluster = int(line[CLUST])
        color = int(line[COLOR])
        return cls(iatom, x, y, z, radius, jatom, cluster, color)

    def __eq__(self, other):
        """Test equality between two Spheres

        Parameters
        ----------
        other : Sphere
            Another Sphere instance to compare.

        Returns
        -------
        result : bool
            Evaluates to true if all the coordinates are close.
        """
        return all(isclose(self.coord, other.coord))

    def __str__(self):
        return "Sphere({0:.5f}, {1:.5f}, {2:.5f})".format(
            self.x, self.y, self.z
        )

    def __repr__(self):
        return "<{0}>".format(self.__str__())

    def __hash__(self):
        return hash(self.__str__())

    def sph_line(self):
        return (
            "{0:5d}{1:10.5f}{2:10.5f}{3:10.5f}{4:8.3f}{5:5d}{6:2d}{7:3d}\n"
        ).format(self.iatom, self.x, self.y, self.z, self.radius,
                 self.jatom, self.cluster, self.color)

    def pdb_line(self):
        return (
            "ATOM   {0:>4d}  C   SPH          {1:8.3f}{2:8.3f}{3:8.3f}\nTER\n"
        ).format(self.iatom, self.x, self.y, self.z)


class SphereFile(object):
    """Parse and manipulate sph files.
    """
    def __init__(self, filename, header, lines):
        self.filename = filename
        self.header = header
        self.lines = lines
        self.matrix = np.array([s.coord for s in lines])

    @classmethod
    def from_file(cls, filename):
        with open(filename) as f:
            header = next(f)
            lines = [Sphere.from_line(line) for line in f]
        return cls(filename, header, lines)

    def write_sph(self, outfile=None):
        if not outfile:
            prefix, ext = os.path.splitext(self.filename)
            outfile = prefix + '.extended' + ext
        with open(outfile, 'w') as f:
            f.write("cluster 1 number of spheres in cluster {0}".format(
                len(self.lines)
            ))
            f.writelines(s.sph_line() for s in self.lines)

    def write_pdb(self, pdbfile=None):
        if not pdbfile:
            prefix, ext = os.path.splitext(self.filename)
            pdbfile = prefix + '.pdb'
        with open(pdbfile, 'w') as f:
            f.writelines(s.pdb_line() for s in self.lines)

    def _centroid(self):
        return self.matrix.mean(axis=0)

    def bounding_box(self, dist=2.8):
        delta = np.array([-dist, dist])
        min_ = self.matrix.min(axis=0)[:, np.newaxis]
        max_ = self.matrix.max(axis=0)[:, np.newaxis]
        minmax = np.concatenate((min_, max_), axis=1) + delta
        box = []
        for (x, y, z) in product(*minmax):
            box.append([x, y, z])
        return np.array(box)

    def __repr__(self):
        return "<SphereFile: {0}".format(self.filename or "BLANK")


def sph2pdb(filename):
    """Convert sphere file to ``pdb``

    Parameters
    ----------
    filename : path
        Path to sphere file

    Returns
    -------
    Writes ``pdb`` to disk
    """
    s = SphereFile.from_file(filename)
    s.write_pdb()


def allsph2pdb():
    """Convert all sphere files in ``cwd`` to ``pdb``
    """
    sphs = glob('*.sph')
    for sph in sphs:
        sph2pdb(sph)


def sphlist2pdb(spheres):
    """Convert a list of sphere files to ``pdb``s

    Parameters
    ----------
    spheres : list
        A list of paths to sphere files

    Returns
    -------
    A ``pdb`` file for each ``sph`` file in spheres
    """
    for sph in spheres:
        sph2pdb(sph)


def split_spheres(filename, n=4):
    """Split sphere into n regions

    Parameters
    ----------
    filename : string
        Path to sphere file to split.

    n : int
        The number of new files to split the input.

    Returns
    -------
    Writes new sphere files to disk.
    """
    s = SphereFile.from_file(filename)
    header = s.header
    X = s.matrix
    km = KMeans(n_clusters=n).fit(X)
    labels = km.labels_
    prefix, ext = os.path.splitext(filename)
    for i in range(n):
        indexes = np.where(labels == i)[0]
        lines = [s.lines[idx] for idx in indexes]
        sphfile = "{0}.{1}{2}".format(prefix, i, ext)
        if os.path.isfile(sphfile):
            raise IOError("File already exists.")
        ns = SphereFile(sphfile, header, lines)
        ns.write_sph(sphfile)


def combine_spheres(outfile, *spheres):
    """Combine multiple sphere files into a single file

    Parameters
    ----------
    outfile : path
        filename to give to combined sphere file

    spheres : paths
        Paths to sphere files to combine

    Returns
    -------
    Writes combined sphere file to disk.
    """
    if not all(os.path.isfile(fn) for fn in spheres):
        raise IOError("Input files must exist.")
    files = [SphereFile.from_file(fn) for fn in spheres]
    s, rest = files[0], files[1:]
    for so in rest:
        s.lines.extend(so.lines)
    s.write_sph(outfile)
