User guide
==========


Examples
========

Here are some examples...

Download structure, select, align, and write
--------------------------------------------

Download a ``pdb``, select desired atoms, align the structures, and save the
structures to ``pdb`` files.

    .. literalinclude:: ../examples/example.py

Now, we can visualize the structures in ``pymol``.

    .. literalinclude:: ../examples/example.pml

This produces the following image

    .. image:: ../examples/example.png


The ``align`` method works selections of equal size (e.g. different conformers
of the same molecules).
It will fail for molecules with different number of atoms.
To align different molecules, use the ``cealign`` method.

    .. literalinclude:: ../examples/example1.py

Visualizing the saved structures gives

    .. image:: ../examples/example1.png


Select binding site residues
----------------------------

Selection helper functions allow quick selections of protein, ligand,
backbone, and side chain atoms.
Distance selections between two objects can be made.
Most methods return a new object, allowing chaining methods together.

    .. literalinclude:: ../examples/example2.py

Visualizing gives

    .. image:: ../examples/example2.png


Structure analysis: Ramachandran plot
-------------------------------------

Create a Ramachandran plot of a protein.

    .. literalinclude:: ../examples/example3.py

    .. image:: ../examples/example3.png


Orient principle axes of a protein along Cartesian axes
-------------------------------------------------------

It is sometimes useful to orient structures such that the centroid lies at
the origin and the principle axes of the protein align with the *x*, *y*, and
*z* axes.

    .. literalinclude:: ../examples/example4.py

Visulaize the transformed structure in ``pymol`` (Cartesian axes indicated by
``rgb`` lines).

    .. image:: ../examples/example4.png
