from find_intersecting_snps import generate_reads
import unittest


class MissTestErrors(unittest.TestCase):

    def test_missing_position(self):
        with self.assertRaises(AssertionError):
            generate_reads(
                read_seq='AAAAAAAA',
                read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
                read_pos=[1],
                ref_alleles=[b'A', b'G'],
                alt_alleles=[b'G', b'C']
            )

    def test_missing_ref(self):
        with self.assertRaises(AssertionError):
            generate_reads(
                read_seq='AAAAAAAA',
                read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
                read_pos=[1, 8],
                ref_alleles=[b'A'],
                alt_alleles=[b'G', b'C']
            )

    def test_missing_alt(self):
        with self.assertRaises(AssertionError):
            generate_reads(
                read_seq='AAAAAAAA',
                read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
                read_pos=[1, 8],
                ref_alleles=[b'A', b'G'],
                alt_alleles=[b'G']
            )

    def test_low_position(self):
        with self.assertRaises(AssertionError):
            generate_reads(
                read_seq='AAAAAAAA',
                read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
                read_pos=[0],
                ref_alleles=[b'A'],
                alt_alleles=[b'G']
            )

    def test_high_position(self):
        with self.assertRaises(AssertionError):
            generate_reads(
                read_seq='AAAAAAAA',
                read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
                read_pos=[9],
                ref_alleles=[b'A'],
                alt_alleles=[b'G']
            )

    def test_extend_position(self):
        with self.assertRaises(AssertionError):
            generate_reads(
                read_seq='AAAAAAAA',
                read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
                read_pos=[8],
                ref_alleles=[b'AA'],
                alt_alleles=[b'G']
            )

    def test_overlapping_variants(self):
        with self.assertRaises(AssertionError):
            generated_reads = generate_reads(
                read_seq='AAAAAAAA',
                read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
                read_pos=[5, 7],
                ref_alleles=[b'AAA', b'C'],
                alt_alleles=[b'GCT', b'T']
            )


class TestSnpMethods(unittest.TestCase):

    def test_single_reference_snp(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
            read_pos=[3],
            ref_alleles=[b'A'],
            alt_alleles=[b'T']
        )
        expected_reads = {
            ('AATAAAAA', (10, 12, 14, 16, 18, 20, 22, 24))
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_nonreference_snp(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
            read_pos=[3],
            ref_alleles=[b'C'],
            alt_alleles=[b'T']
        )
        expected_reads = {
            ('AATAAAAA', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('AACAAAAA', (10, 12, 14, 16, 18, 20, 22, 24))
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_double_snp(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
            read_pos=[2, 6],
            ref_alleles=[b'C', b'A'],
            alt_alleles=[b'T', b'G']
        )
        expected_reads = {
            ('ATAAAAAA', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('ACAAAAAA', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('AAAAAGAA', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('ATAAAGAA', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('ACAAAGAA', (10, 12, 14, 16, 18, 20, 22, 24))
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_double_terminal_snp(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
            read_pos=[1, 8],
            ref_alleles=[b'A', b'G'],
            alt_alleles=[b'G', b'C']
        )
        expected_reads = {
            ('GAAAAAAA', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('AAAAAAAG', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('AAAAAAAC', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('GAAAAAAG', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('GAAAAAAC', (10, 12, 14, 16, 18, 20, 22, 24))
        }
        self.assertEqual(generated_reads, expected_reads)


class TestLargerNonIndelMethods(unittest.TestCase):
    def test_single_reference_variant(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
            read_pos=[3],
            ref_alleles=[b'AA'],
            alt_alleles=[b'TT']
        )
        expected_reads = {
            ('AATTAAAA', (10, 12, 14, 16, 18, 20, 22, 24))
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_nonreference_variant(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
            read_pos=[7],
            ref_alleles=[b'CG'],
            alt_alleles=[b'TT']
        )
        expected_reads = {
            ('AAAAAACG', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('AAAAAATT', (10, 12, 14, 16, 18, 20, 22, 24))
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_multiple_variants(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
            read_pos=[2, 7],
            ref_alleles=[b'AAA', b'CG'],
            alt_alleles=[b'GCT', b'TT']
        )
        expected_reads = {
            ('AGCTAAAA', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('AAAAAACG', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('AAAAAATT', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('AGCTAACG', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('AGCTAATT', (10, 12, 14, 16, 18, 20, 22, 24))
        }
        self.assertEqual(generated_reads, expected_reads)


class TestDeletionMethods(unittest.TestCase):

    def test_single_reference_deletion(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
            read_pos=[3],
            ref_alleles=[b'AA'],
            alt_alleles=[b'A']
        )
        expected_reads = {
            ('AAAAAAA', (10, 12, 15, 18, 20, 22, 24))
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_nonreference_deletion(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
            read_pos=[1],
            ref_alleles=[b'TAT'],
            alt_alleles=[b'G']
        )
        expected_reads = {
            ('TATAAAAA', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('GAAAAA', (12, 16, 18, 20, 22, 24))
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_multiple_deletions(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
            read_pos=[2, 6],
            ref_alleles=[b'AAA', b'TAC'],
            alt_alleles=[b'TG', b'T']
        )
        expected_reads = {
            ('ATGAAAA', '0224567'),
            ('AAAAATAC', '01234567'),
            ('AAAAAT', '012346'),
            ('ATGATAC', '0224567'),
            ('ATGAT', '02246')
        }
        expected_reads = {
            ('ATGAAAA', (10, 14, 14, 18, 20, 22, 24)),
            ('AAAAATAC', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('AAAAAT', (10, 12, 14, 16, 18, 22)),
            ('ATGATAC', (10, 14, 14, 18, 20, 22, 24)),
            ('ATGAT', (10, 14, 14, 18, 22))
        }
        self.assertEqual(generated_reads, expected_reads)

class TestInsertionMethods(unittest.TestCase):
    
    def test_single_reference_insertion(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
            read_pos=[3],
            ref_alleles=[b'A'],
            alt_alleles=[b'TGAT']
        )
        expected_reads = {
            ('AATGATAAAAA', (10, 12, 14, 14, 14, 14, 16, 18, 20, 22, 24))
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_nonreference_insertion(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
            read_pos=[1],
            ref_alleles=[b'TAT'],
            alt_alleles=[b'GGGCC']
        )
        expected_reads = {
            ('TATAAAAA', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('GGGCCAAAAA', (12, 12, 12, 12, 12, 16, 18, 20, 22, 24))
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_multiple_insertions(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual=[10, 12, 14, 16, 18, 20, 22, 24],
            read_pos=[2, 6],
            ref_alleles=[b'A', b'TAC'],
            alt_alleles=[b'TGG', b'GGCCC']
        )
        expected_reads = {
            ('ATGGAAAAAA', (10, 12, 12, 12, 14, 16, 18, 20, 22, 24)),
            ('AAAAATAC', (10, 12, 14, 16, 18, 20, 22, 24)),
            ('AAAAAGGCCC', (10, 12, 14, 16, 18, 22, 22, 22, 22, 22)),
            ('ATGGAAATAC', (10, 12, 12, 12, 14, 16, 18, 20, 22, 24)),
            ('ATGGAAAGGCCC', (10, 12, 12, 12, 14, 16, 18, 22, 22, 22, 22, 22)),
        }
        self.assertEqual(generated_reads, expected_reads)

if __name__ == '__main__':
    unittest.main()
