from protutils.pdb import PDBFile

pdb = PDBFile.fetch('4K5Y')
protein = pdb.select(chain='A', nres__lt=1000).protein()

orient = protein.orient()

protein.write_pdb('4K5YA.pdb')
orient.write_pdb('4K5YA_0.pdb')
