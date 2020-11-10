from find_intersecting_snps import generate_reads
import unittest


class MissTestErrors(unittest.TestCase):

    def test_missing_position(self):
        with self.assertRaises(AssertionError):
            generate_reads(
                read_seq='AAAAAAAA',
                read_qual='01234567',
                read_pos=[1],
                ref_alleles=[b'A', b'G'],
                alt_alleles=[b'G', b'C']
            )

    def test_missing_ref(self):
        with self.assertRaises(AssertionError):
            generate_reads(
                read_seq='AAAAAAAA',
                read_qual='01234567',
                read_pos=[1, 8],
                ref_alleles=[b'A'],
                alt_alleles=[b'G', b'C']
            )

    def test_missing_alt(self):
        with self.assertRaises(AssertionError):
            generate_reads(
                read_seq='AAAAAAAA',
                read_qual='01234567',
                read_pos=[1, 8],
                ref_alleles=[b'A', b'G'],
                alt_alleles=[b'G']
            )

    def test_low_position(self):
        with self.assertRaises(AssertionError):
            generate_reads(
                read_seq='AAAAAAAA',
                read_qual='01234567',
                read_pos=[0],
                ref_alleles=[b'A'],
                alt_alleles=[b'G']
            )

    def test_high_position(self):
        with self.assertRaises(AssertionError):
            generate_reads(
                read_seq='AAAAAAAA',
                read_qual='01234567',
                read_pos=[9],
                ref_alleles=[b'A'],
                alt_alleles=[b'G']
            )

    def test_extend_position(self):
        with self.assertRaises(AssertionError):
            generate_reads(
                read_seq='AAAAAAAA',
                read_qual='01234567',
                read_pos=[8],
                ref_alleles=[b'AA'],
                alt_alleles=[b'G']
            )

    def test_overlapping_variants(self):
        with self.assertRaises(AssertionError):
            generated_reads = generate_reads(
                read_seq='AAAAAAAA',
                read_qual='01234567',
                read_pos=[5, 7],
                ref_alleles=[b'AAA', b'C'],
                alt_alleles=[b'GCT', b'T']
            )


class TestSnpMethods(unittest.TestCase):

    def test_single_reference_snp(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual='01234567',
            read_pos=[3],
            ref_alleles=[b'A'],
            alt_alleles=[b'T']
        )
        expected_reads = {
            ('AATAAAAA', '01234567')
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_nonreference_snp(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual='01234567',
            read_pos=[3],
            ref_alleles=[b'C'],
            alt_alleles=[b'T']
        )
        expected_reads = {
            ('AATAAAAA', '01234567'),
            ('AACAAAAA', '01234567')
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_double_snp(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual='01234567',
            read_pos=[2, 6],
            ref_alleles=[b'C', b'A'],
            alt_alleles=[b'T', b'G']
        )
        expected_reads = {
            ('ATAAAAAA', '01234567'),
            ('ACAAAAAA', '01234567'),
            ('AAAAAGAA', '01234567'),
            ('ATAAAGAA', '01234567'),
            ('ACAAAGAA', '01234567')
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_double_terminal_snp(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual='01234567',
            read_pos=[1, 8],
            ref_alleles=[b'A', b'G'],
            alt_alleles=[b'G', b'C']
        )
        expected_reads = {
            ('GAAAAAAA', '01234567'),
            ('AAAAAAAG', '01234567'),
            ('AAAAAAAC', '01234567'),
            ('GAAAAAAG', '01234567'),
            ('GAAAAAAC', '01234567')
        }
        self.assertEqual(generated_reads, expected_reads)


class TestLargerNonIndelMethods(unittest.TestCase):
    def test_single_reference_variant(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual='01234567',
            read_pos=[3],
            ref_alleles=[b'AA'],
            alt_alleles=[b'TT']
        )
        expected_reads = {
            ('AATTAAAA', '01234567')
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_nonreference_variant(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual='01234567',
            read_pos=[7],
            ref_alleles=[b'CG'],
            alt_alleles=[b'TT']
        )
        expected_reads = {
            ('AAAAAACG', '01234567'),
            ('AAAAAATT', '01234567')
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_multiple_variants(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual='01234567',
            read_pos=[2, 7],
            ref_alleles=[b'AAA', b'CG'],
            alt_alleles=[b'GCT', b'TT']
        )
        expected_reads = {
            ('AGCTAAAA', '01234567'),
            ('AAAAAACG', '01234567'),
            ('AAAAAATT', '01234567'),
            ('AGCTAACG', '01234567'),
            ('AGCTAATT', '01234567')
        }
        self.assertEqual(generated_reads, expected_reads)


class TestDeletionMethods(unittest.TestCase):
    
    def test_single_reference_deletion(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual='01234567',
            read_pos=[3],
            ref_alleles=[b'AA'],
            alt_alleles=[b'A']
        )
        expected_reads = {
            ('AAAAAAA', '0124567')
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_nonreference_deletion(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual='01234567',
            read_pos=[1],
            ref_alleles=[b'TAT'],
            alt_alleles=[b'G']
        )
        expected_reads = {
            ('TATAAAAA', '01234567'),
            ('GAAAAA', '134567')
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_multiple_deletions(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual='01234567',
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
        self.assertEqual(generated_reads, expected_reads)

class TestInsertionMethods(unittest.TestCase):
    
    def test_single_reference_insertion(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual='01234567',
            read_pos=[3],
            ref_alleles=[b'A'],
            alt_alleles=[b'TGAT']
        )
        expected_reads = {
            ('AATGATAAAAA', '01222234567')
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_nonreference_insertion(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual='01234567',
            read_pos=[1],
            ref_alleles=[b'TAT'],
            alt_alleles=[b'GGGCC']
        )
        expected_reads = {
            ('TATAAAAA', '01234567'),
            ('GGGCCAAAAA', '1111134567')
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_multiple_insertions(self):
        generated_reads = generate_reads(
            read_seq='AAAAAAAA',
            read_qual='01234567',
            read_pos=[2, 6],
            ref_alleles=[b'A', b'TAC'],
            alt_alleles=[b'TGG', b'GGCCC']
        )
        expected_reads = {
            ('ATGGAAAAAA', '0111234567'),
            ('AAAAATAC', '01234567'),
            ('AAAAAGGCCC', '0123466666'),
            ('ATGGAAATAC', '0111234567'),
            ('ATGGAAAGGCCC', '011123466666'),
        }
        self.assertEqual(generated_reads, expected_reads)

if __name__ == '__main__':
    unittest.main()
