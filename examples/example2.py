from protutils.pdb import PDBFile


# write a small helper function
def residues_for_selection(selected, pdbfile):
    """Get all residues atoms for a selection
    """
    residues = {atm.nres for atm in selected}
    return pdbfile.select(nres__isin=residues)


pdb = PDBFile.fetch('1HPV')

# Select ligand
ligand = pdb.ligand()

# select protein atoms with 5 Angstroms of the ligand
atoms = pdb.protein().within(5.0, ligand)

prot = residues_for_selection(atoms, pdb)

ligand.write_pdb('ligand.pdb')
prot.write_pdb('prot.pdb')
