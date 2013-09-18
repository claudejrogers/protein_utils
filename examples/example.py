from protutils.pdb import PDBFile


# Download 4k5Y.pdb from the RCSB Protein Data Bank
pdb = PDBFile.fetch('4K5Y')

# Select residues less than 263 on chain A, then chain B
chain_a = pdb.select(chain__eq='A', nres__lt=263)
chain_b = pdb.select(chain__eq='B', nres__lt=263)

# Align selection with another selection that contains the same number of atoms
aligned_chain_b = chain_b.align(chain_a)
# prints: RMSD = 0.983960869568

# write structures to pdb files
chain_a.write_pdb('chain_a.pdb')
chain_b.write_pdb('chain_b.pdb')
aligned_chain_b.write_pdb('aligned_chain_b.pdb')
