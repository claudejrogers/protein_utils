import unittest

from utils.pdb import PDBAtom, PDBFile


class TestPDBAtomMethods(unittest.TestCase):

    def setUp(self):
        self.line = ("ATOM    158  C   LEU A 135     -42.160 -11.386  54.208"
                     "  1.00 43.54           C  \n")
        self.het = ("HETATM 8488  C13 1Q5 A 401     -41.176 -19.302  34.303 "
                    " 1.00 48.59           C  \n")
        self.conect = ("CONECT 8488 8481 8487 8503                          "
                       "                            \n")
        self.atom = PDBAtom.from_line(self.line)
        self.hetatm = PDBAtom.from_line(self.het)
        self.hetatm.add_connections_from_line(self.conect)

    def test_pdb_init_from_line(self):
        self.assertIsInstance(self.atom, PDBAtom,
                              msg="PDBAtom.from_line failed")

    def test_pdb_parse_line(self):
        pdb = self.atom
        self.assertEqual(pdb.record, 'ATOM', msg="Parse RECORD failed")
        self.assertEqual(pdb.natom, 158, msg="Parse ATOM ID failed")
        self.assertEqual(pdb.atom, 'C', msg="Parse ATOM failed")
        self.assertEqual(pdb.altloc, ' ', msg="Parse ALTLOC failed")
        self.assertEqual(pdb.res, 'LEU', msg="Parse RES failed")
        self.assertEqual(pdb.chain, 'A', msg="Parse CHAIN failed")
        self.assertEqual(pdb.nres, 135, msg="Parse RES ID failed")
        self.assertEqual(pdb.icode, ' ', msg="Parse ICODE failed")
        self.assertAlmostEqual(pdb.x, -42.160, places=3, msg="Parse X failed")
        self.assertAlmostEqual(pdb.y, -11.386, places=3, msg="Parse Y failed")
        self.assertAlmostEqual(pdb.z, 54.208, places=3, msg="Parse Z failed")
        self.assertAlmostEqual(pdb.occupancy, 1.00, places=2,
                               msg="Parse OCCUPANCY failed")
        self.assertAlmostEqual(pdb.bfactor, 43.54, places=2,
                               msg="Parse BFACTOR failed")
        self.assertEqual(pdb.element, 'C', msg="Parse ELEMENT failed")
        self.assertEqual(pdb.charge, '  ', msg="Parse CHARGE failed")

    def test_pdb_add_connections_from_line(self):
        pdb = self.hetatm
        self.assertListEqual(pdb.connections, [8481, 8487, 8503],
                             msg="Parse CONECT failed")

    def test_pdb_atom_distance_calculation(self):
        dist = self.atom.distance(self.hetatm)
        distance = sum(
            (x - y) ** 2 for x, y in zip([-42.160, -11.386, 54.208],
                                         [-41.176, -19.302, 34.303])
        ) ** 0.5
        self.assertAlmostEqual(dist, distance, places=5,
                               msg="Distance calculation failed")

    def test_pdb_writeline(self):
        line, conect = self.atom.writeline()
        self.assertEqual(line, self.line, msg="PDBAtom.writeline failed")
        self.assertEqual(conect, '', msg="PDBAtom.writeline failed")

    def test_pdb_writeline_with_conect(self):
        line, conect = self.hetatm.writeline()
        test_conect = self.conect.strip() + '\n'
        self.assertEqual(line, self.het, msg="PDBAtom.writeline failed")
        self.assertEqual(conect, test_conect,
                         msg="CONECT PDBAtom.writeline failed")

    def test_pdb_truncated_line_ok(self):
        line = self.line.strip()
        pdb = PDBAtom.from_line(line)
        self.assertEqual(pdb.charge, '',
                         msg="Parse CHARGE failed on truncated line")

    def test_pdb_repr(self):
        r = repr(self.atom)
        self.assertEqual(r, "<PDB: ATOM    158  C   LEU A 135 ...>",
                         msg="PDBAtom.__repr__ failed")


class TestPDBFileMethods(unittest.TestCase):

    def setUp(self):
        self.line = ("ATOM    158  C   LEU A 135     -42.160 -11.386  54.208"
                     "  1.00 43.54           C  \n")
        self.web = PDBFile.fetch('4K5Y')
        self.file_ = PDBFile.from_file('tests/resources/4K5Y.pdb')

    def test_pdb_file_from_file_inits(self):
        self.assertIsInstance(self.file_, PDBFile, msg="web fetch failed")

    def test_pdb_file_fetch_inits(self):
        self.assertIsInstance(self.web, PDBFile, msg="web fetch failed")

    def test_pdb_file_parse(self):
        line = self.file_[157].writeline()[0]
        print repr(line)
        print repr(self.line)
        self.assertEqual(line, self.line, msg="File parse failed")

    def test_len_of_pdb(self):
        length = len(self.file_)
        self.assertEqual(length, 8823, msg="Length calculation failed")

    def test_pdb_chain_property(self):
        chains = self.file_.chains
        self.assertListEqual(chains, ['A', 'B', 'C'],
                             msg="Chain property failed")


class TestPDBFileBitwiseOperations(unittest.TestCase):

    def setUp(self):
        self.pdb = PDBFile.from_file('tests/resources/4K5Y.pdb')
        self.chain_a = self.pdb.select(chain__eq='A', nres__lt=263)

    def test_select_chain(self):
        chain_a = self.pdb.select(chain__eq='A')
        self.assertIsInstance(chain_a, PDBFile,
                              msg="PDBFile not returned by select")
        self.assertNotEqual(chain_a, self.pdb,
                            msg="New object not created by select")
        self.assertListEqual(chain_a.chains, ['A'], msg="Selection failed")

    def test_dunder_contains(self):
        self.assertIn(self.chain_a, self.pdb)

    def test_dunder_and_method(self):
        cha = self.pdb & self.chain_a
        self.assertEqual(cha, self.chain_a)

    def test_dunder_or_method(self):
        chain_a = self.pdb.select(chain__eq='A')
        chain_bc = self.pdb.select(chain__gt='A')
        all_ = chain_a | chain_bc
        self.assertEqual(all_, self.pdb)

    def test_dunder_xor_method(self):
        chain_a = self.pdb.select(chain__eq='A')
        chain_bc = self.pdb.select(chain__gt='A')
        cbc = self.pdb ^ chain_a
        self.assertEqual(cbc, chain_bc)

    def test_dunder_iand(self):
        chain_a = self.pdb.select(chain__eq='A', nres__lt=263)
        chain_a &= self.pdb
        self.assertEqual(chain_a, self.chain_a)

    def test_dunder_ior(self):
        chain_a = self.pdb.select(chain='A')
        chain_bc = self.pdb.select(chain__gt='A')
        chain_a |= chain_bc
        self.assertEqual(chain_a, self.pdb)

    def test_dunder_ixor(self):
        chain_a = self.pdb.select(chain='A')
        chain_bc = self.pdb.select(chain__gt='A')
        chain_a ^= self.pdb
        self.assertEqual(chain_a, chain_bc)


class TestWithinMethod(unittest.TestCase):

    def setUp(self):
        self.pdb = PDBFile.from_file('tests/resources/1hpv.pdb')
        self.lig = self.pdb.ligand()
        self.prot = self.pdb.protein()

    def test_distance_method(self):
        binding_site = self.prot.within(5.0, self.lig)
        self.assertEqual(len(binding_site), 100)
