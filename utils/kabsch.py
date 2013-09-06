import numpy as np


def _get_centered_coordinates(coord_array):
    """
    Translate molecule to be centered at 0, 0, 0.
    """
    COM1 = np.average(coord_array, axis=0)
    return (coord_array - COM1, COM1)


def _kabsch(mol1, mol2, rmsd=False):
    """
    Align other molecule to the top. Returns new coordinates for other.
    """
    L = len(mol1)

    top, COM1 = _get_centered_coordinates(mol1)
    other, COM2 = _get_centered_coordinates(mol2)

    assert L == len(other), "molecules must have the same number of atoms"

    # initial residual
    E0 = (np.sum(np.sum(top * top, axis=0), axis=0) +
          np.sum(np.sum(other * other, axis=0), axis=0))

    # Kabsch
    V, S, Wt = np.linalg.svd(np.dot(other.T, top))

    reflect = float(np.linalg.det(V) * np.linalg.det(Wt))

    if reflect == -1:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    RMSD = E0 - (2.0 * sum(S))
    RMSD = np.sqrt(abs(RMSD / L))
    if rmsd:
        return RMSD

    print "RMSD =", RMSD

    U = np.dot(V, Wt)

    # rotate and translate the molecule
    mol2 = np.dot((other), U)
    return mol2 + COM1


def aligned_rmsd(mol1, mol2):
    return _kabsch(mol1, mol2, rmsd=True)


def align(mol1, mol2):
    return _kabsch(mol1, mol2, rmsd=False)
