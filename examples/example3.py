from protutils.pdb import PDBFile


pdb = PDBFile.fetch('4K5Y')
protein = pdb.protein()  # remove HETATM records
protein.ramachandran_plot()
