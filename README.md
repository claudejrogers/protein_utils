# Protein Utilities

Utilities to manipulate PDB, BGF, (mol2, mae, etc.) files.

## Examples

### Selections and alignment

```python
>>> from utils.pdb import PDBFile

# Download 4K5Y.pdb from RCSB Protein Data Bank
>>> pdb = PDBFile.fetch('4K5Y')
# Select residues less than 263 on chain A
>>> chain_a = pdb.select(chain__eq='A', nres__lt=263)
>>> chain_b = pdb.select(chain__eq='B', nres__lt=263)
>>> aligned_chain_b = chain_b.align(chain_a)
RMSD = 0.983960869568
>>> chain_a.write_pdb('chain_a.pdb')
>>> chain_b.write_pdb('chain_b.pdb')
>>> aligned_chain_b = chain_b.write_pdb('aligned_chain_b.pdb')
```

Now, we can visualize our selections:

```
# example.pml
load chain_a.pdb
load chain_b.pdb
load aligned_chain_b.pdb

zoom
show cartoon
ray
png example.png
```

to get:

![pymol img](examples/example.png)
