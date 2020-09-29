from find_intersecting_snps import generate_reads
import unittest


class TestErrors(unittest.TestCase):

    def test_missing_position(self):
        with self.assertRaises(TypeError):
            generated_reads = generate_reads(
                read_seq='AAAAAAAA',
                read_pos=[1],
                ref_alleles=[b'A', b'G'],
                alt_alleles=[b'G', b'C']
            )

    def test_missing_ref(self):
        with self.assertRaises(TypeError):
            generated_reads = generate_reads(
                read_seq='AAAAAAAA',
                read_pos=[1, 8],
                ref_alleles=[b'A'],
                alt_alleles=[b'G', b'C']
            )

    def test_missing_alt(self):
        with self.assertRaises(TypeError):
            generated_reads = generate_reads(
                read_seq='AAAAAAAA',
                read_pos=[1, 8],
                ref_alleles=[b'A', b'G'],
                alt_alleles=[b'G']
            )

    def test_low_position(self):
        with self.assertRaises(AssertionError):
            generated_reads = generate_reads(
                read_seq='AAAAAAAA',
                read_pos=[0, 1],
                ref_alleles=[b'A', b'G'],
                alt_alleles=[b'G', b'C']
            )
            print(generated_reads)

    def test_high_position(self):
        with self.assertRaises(AssertionError):
            generated_reads = generate_reads(
                read_seq='AAAAAAAA',
                read_pos=[1, 9],
                ref_alleles=[b'A', b'G'],
                alt_alleles=[b'G', b'C']
            )
            print(generated_reads)


class TestSnpMethods(unittest.TestCase):

    def test_single_reference_snp(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_pos=[3],
            ref_alleles=[b'A'],
            alt_alleles=[b'T']
        )
        expected_reads = {
            'AAAAAAAA',
            'AATAAAAA'
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_nonreference_snp(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_pos=[3],
            ref_alleles=[b'C'],
            alt_alleles=[b'T']
        )
        expected_reads = {
            'AAAAAAAA',
            'AATAAAAA',
            'AACAAAAA'
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_double_snp(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_pos=[2, 6],
            ref_alleles=[b'C', b'A'],
            alt_alleles=[b'T', b'G']
        )
        expected_reads = {
            'AAAAAAAA',
            'ATAAAAAA',
            'ACAAAAAA',
            'AAAAAGAA',
            'ATAAAGAA',
            'ACAAAGAA',
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_double_terminal_snp(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_pos=[1, 8],
            ref_alleles=[b'A', b'G'],
            alt_alleles=[b'G', b'C']
        )
        expected_reads = {
            'AAAAAAAA',
            'GAAAAAAA',
            'AAAAAAAG',
            'AAAAAAAC',
            'GAAAAAAG',
            'GAAAAAAC',
        }
        self.assertEqual(generated_reads, expected_reads)


class TestAdjacentSnpMethods(unittest.TestCase):

    def test_single_reference_variant(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_pos=[3],
            ref_alleles=[b'AA'],
            alt_alleles=[b'TT']
        )
        expected_reads = {
            'AAAAAAAA',
            'AATTAAAA'
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_nonreference_variant(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_pos=[7],
            ref_alleles=[b'CG'],
            alt_alleles=[b'TT']
        )
        expected_reads = {
            'AAAAAAAA',
            'AAAAAACG',
            'AAAAAATT',
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_nonreference_variant(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_pos=[7],
            ref_alleles=[b'CG'],
            alt_alleles=[b'TT']
        )
        expected_reads = {
            'AAAAAAAA',
            'AAAAAACG',
            'AAAAAATT',
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_multiple_variants(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_pos=[2, 7],
            ref_alleles=[b'AAA', b'CG'],
            alt_alleles=[b'GCT', b'TT']
        )
        expected_reads = {
            'AAAAAAAA',
            'AGCTAAAA',
            'AAAAAACG',
            'AAAAAATT',
            'AGCTAACG',
            'AGCTAATT',
        }
        self.assertEqual(generated_reads, expected_reads)

if __name__ == '__main__':
    unittest.main()
