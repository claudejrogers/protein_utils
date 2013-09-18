from protutils.pdb import PDBFile


pdb = PDBFile.fetch('4K5Y')
chain_a = pdb.select(chain='A')  # equivalent to pdb.select(chain__eq='A')
chain_b = pdb.select(chain='B')

# compare the number of atoms in the selection
len(chain_a) == len(chain_b)
# returns False

# the align method would fail for these selections
aligned_chain_b = chain_b.cealign(chain_a)
# prints RMSD = 0.874903919378  (RMSD of alpha carbons)

# write structures
chain_a.write_pdb('chain_a.pdb')
chain_b.write_pdb('chain_b.pdb')
aligned_chain_b.write_pdb('aligned_chain_b.pdb')
