from protutils.pdb import PDBFile
from protutils.ncbi.blastp import BLASTPDBRecord

pdb = PDBFile.fetch('4K5Y')
chain_a = pdb.select(chain='A', nres__lt=1000)
# Get sequence
sequence = chain_a.sequence.replace('-', '')

# search for similar sequences to 4K5Y_A
query = BLASTPDBRecord(sequence)

# Get top hit
PDB = query.get_best()['pdb']
similar = PDBFile.fetch(PDB)
similar.select(chain='A', nres__lt=1000)

# Compare aligned structures
aligned = similar.cealign(chain_a)

chain_a.write_pdb('4k5y_A.mod.pdb')
aligned.write_pdb('{0}_A.mod.pdb'.format(PDB))
