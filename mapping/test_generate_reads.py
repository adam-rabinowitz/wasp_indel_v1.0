from find_intersecting_snps import generate_reads
import unittest


# class MissTestErrors(unittest.TestCase):

#     def test_low_position(self):
#         with self.assertRaises(AssertionError):
#             generate_reads(
#                 read_seq='AAAAAAAA',
#                 read_qual=[20, 22, 24, 26, 28, 30, 32, 34],
#                 read_pos=[0],
#                 ref_alleles=[b'A'],
#                 alt_alleles=[b'G']
#             )

#     def test_high_position(self):
#         with self.assertRaises(AssertionError):
#             generate_reads(
#                 read_seq='AAAAAAAA',
#                 read_qual=[20, 22, 24, 26, 28, 30, 32, 34],
#                 read_pos=[9],
#                 ref_alleles=[b'A'],
#                 alt_alleles=[b'G']
#             )

#     def test_extend_position(self):
#         with self.assertRaises(AssertionError):
#             generate_reads(
#                 read_seq='AAAAAAAA',
#                 read_qual=[20, 22, 24, 26, 28, 30, 32, 34],
#                 read_pos=[8],
#                 ref_alleles=[b'AA'],
#                 alt_alleles=[b'G']
#             )

#     def test_overlapping_variants(self):
#         with self.assertRaises(AssertionError):
#             generated_reads = generate_reads(
#                 read_seq='AAAAAAAA',
#                 read_qual=[20, 22, 24, 26, 28, 30, 32, 34],
#                 read_pos=[5, 7],
#                 ref_alleles=[b'AAA', b'C'],
#                 alt_alleles=[b'GCT', b'T']
#             )


class TestSnpMethods(unittest.TestCase):

    def test_single_reference_snp(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(2, 3, 'A', 'G')]
        )
        expected_reads = {
            'AAGAAAAA': [20, 22, 24, 26, 28, 30, 32, 34]
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_nonreference_snp(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(6, 7, 'T', 'C,G,A')]
        )
        expected_reads = {
            'AAAAAATA': [20, 22, 24, 26, 28, 30, 32, 34],
            'AAAAAACA': [20, 22, 24, 26, 28, 30, 32, 34],
            'AAAAAAGA': [20, 22, 24, 26, 28, 30, 32, 34]
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_start_snp(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(0, 1, 'C', 'G')]
        )
        expected_reads = {
            'CAAAAAAA': [20, 22, 24, 26, 28, 30, 32, 34],
            'GAAAAAAA': [20, 22, 24, 26, 28, 30, 32, 34]
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_end_snp(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(7, 8, 'A', 'G,T')]
        )
        expected_reads = {
            'AAAAAAAG': [20, 22, 24, 26, 28, 30, 32, 34],
            'AAAAAAAT': [20, 22, 24, 26, 28, 30, 32, 34]
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_double_snp(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(1, 2, 'C', 'T'), (5, 6, 'A', 'G')]
        )
        expected_reads = {
            'ATAAAAAA': [20, 22, 24, 26, 28, 30, 32, 34],
            'ACAAAAAA': [20, 22, 24, 26, 28, 30, 32, 34],
            'AAAAAGAA': [20, 22, 24, 26, 28, 30, 32, 34],
            'ATAAAGAA': [20, 22, 24, 26, 28, 30, 32, 34],
            'ACAAAGAA': [20, 22, 24, 26, 28, 30, 32, 34]
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_double_neighbouring_snp(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(3, 4, 'C', 'T'), (4, 5, 'A', 'G,T')]
        )
        expected_reads = {
            'AAATAAAA': [20, 22, 24, 26, 28, 30, 32, 34],
            'AAACAAAA': [20, 22, 24, 26, 28, 30, 32, 34],
            'AAAAGAAA': [20, 22, 24, 26, 28, 30, 32, 34],
            'AAAATAAA': [20, 22, 24, 26, 28, 30, 32, 34],
            'AAATGAAA': [20, 22, 24, 26, 28, 30, 32, 34],
            'AAATTAAA': [20, 22, 24, 26, 28, 30, 32, 34],
            'AAACGAAA': [20, 22, 24, 26, 28, 30, 32, 34],
            'AAACTAAA': [20, 22, 24, 26, 28, 30, 32, 34]
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_double_terminal_snp(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(0, 1, 'A', 'G'), (7, 8, 'G', 'C')]
        )
        expected_reads = {
            'GAAAAAAA': [20, 22, 24, 26, 28, 30, 32, 34],
            'AAAAAAAG': [20, 22, 24, 26, 28, 30, 32, 34],
            'AAAAAAAC': [20, 22, 24, 26, 28, 30, 32, 34],
            'GAAAAAAG': [20, 22, 24, 26, 28, 30, 32, 34],
            'GAAAAAAC': [20, 22, 24, 26, 28, 30, 32, 34]
        }
        self.assertEqual(generated_reads, expected_reads)


class TestLargerNonIndelMethods(unittest.TestCase):

    def test_single_reference_variant(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(2, 4, 'AA', 'TT,TG')]
        )
        expected_reads = {
            'AATTAAAA': [20, 22, 25, 25, 28, 30, 32, 34],
            'AATGAAAA': [20, 22, 25, 25, 28, 30, 32, 34]
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_end_nonreference_variant(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(6, 8, 'CG', 'TT')]
        )
        expected_reads = {
            'AAAAAACG': [20, 22, 24, 26, 28, 30, 33, 33],
            'AAAAAATT': [20, 22, 24, 26, 28, 30, 33, 33]
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_multiple_variants(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(0, 3, 'AAA', 'GCT,TCT'), (6, 8, 'CG', 'TT')]
        )
        expected_reads = {
            'GCTAAAAA': [22, 22, 22, 26, 28, 30, 32, 34],
            'TCTAAAAA': [22, 22, 22, 26, 28, 30, 32, 34],
            'AAAAAACG': [20, 22, 24, 26, 28, 30, 33, 33],
            'AAAAAATT': [20, 22, 24, 26, 28, 30, 33, 33],
            'GCTAAACG': [22, 22, 22, 26, 28, 30, 33, 33],
            'TCTAAACG': [22, 22, 22, 26, 28, 30, 33, 33],
            'GCTAAATT': [22, 22, 22, 26, 28, 30, 33, 33],
            'TCTAAATT': [22, 22, 22, 26, 28, 30, 33, 33]
        }
        self.assertEqual(generated_reads, expected_reads)


class TestDeletionMethods(unittest.TestCase):

    def test_single_reference_deletion(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(2, 4, 'AA', 'A')]
        )
        expected_reads = {
            'AAAAAAA': [20, 22, 25, 28, 30, 32, 34]
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_nonreference_deletion(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(0, 3, 'TAT', 'G,CC,')]
        )
        expected_reads = {
            'TATAAAAA': [22, 22, 22, 26, 28, 30, 32, 34],
            'GAAAAA': [22, 26, 28, 30, 32, 34],
            'CCAAAAA': [22, 22, 26, 28, 30, 32, 34],
            'AAAAA': [26, 28, 30, 32, 34]
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_multiple_deletions(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(1, 3, 'TG', 'T'), (5, 8, 'AAA', 'TC,GC')]
        )
        expected_reads = {
            'ATGAAAAA': [20, 23, 23, 26, 28, 30, 32, 34],
            'ATAAAAA': [20, 23, 26, 28, 30, 32, 34],
            'AAAAATC': [20, 22, 24, 26, 28, 32, 32],
            'AAAAAGC': [20, 22, 24, 26, 28, 32, 32],
            'ATGAATC': [20, 23, 23, 26, 28, 32, 32],
            'ATGAAGC': [20, 23, 23, 26, 28, 32, 32],
            'ATAATC': [20, 23, 26, 28, 32, 32],
            'ATAAGC': [20, 23, 26, 28, 32, 32],
        }
        self.assertEqual(generated_reads, expected_reads)


class TestInsertionMethods(unittest.TestCase):

    def test_single_reference_insertion(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(2, 3, 'A', 'TGAT,TG')]
        )
        expected_reads = {
            'AATGATAAAAA': [20, 22, 24, 24, 24, 24, 26, 28, 30, 32, 34],
            'AATGAAAAA': [20, 22, 24, 24, 26, 28, 30, 32, 34]
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_nonreference_insertion(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(0, 3, 'TAT', 'GGGCC')]
        )
        expected_reads = {
            'TATAAAAA': [22, 22, 22, 26, 28, 30, 32, 34],
            'GGGCCAAAAA': [22, 22, 22, 22, 22, 26, 28, 30, 32, 34]
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_absent_reference_start(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(0, 0, 'T', 'CC')]
        )
        expected_reads = {
            'TAAAAAAAA': [20, 20, 22, 24, 26, 28, 30, 32, 34],
            'CCAAAAAAAA': [20, 20, 20, 22, 24, 26, 28, 30, 32, 34]
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_single_absent_reference_end(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(8, 8, 'T', 'CC')]
        )
        expected_reads = {
            'AAAAAAAAT': [20, 22, 24, 26, 28, 30, 32, 34, 34],
            'AAAAAAAACC': [20, 22, 24, 26, 28, 30, 32, 34, 34, 34]
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_terminal_insertion(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(5, 8, 'TAC', 'GGCCC')]
        )
        expected_reads = {
            'AAAAATAC': [20, 22, 24, 26, 28, 32, 32, 32],
            'AAAAAGGCCC': [20, 22, 24, 26, 28, 32, 32, 32, 32, 32]
        }
        self.assertEqual(generated_reads, expected_reads)

    def test_multiple_insertions(self):
        generated_reads = generate_reads(
            sequence='AAAAAAAA',
            quality=[20, 22, 24, 26, 28, 30, 32, 34],
            variants=[(0, 1, 'A', 'TGG'), (5, 8, 'TAC', 'GGCCC')]
        )
        expected_reads = {
            'TGGAAAAAAA': [20, 20, 20, 22, 24, 26, 28, 30, 32, 34],
            'AAAAATAC': [20, 22, 24, 26, 28, 32, 32, 32],
            'AAAAAGGCCC': [20, 22, 24, 26, 28, 32, 32, 32, 32, 32],
            'TGGAAAATAC': [20, 20, 20, 22, 24, 26, 28, 32, 32, 32],
            'TGGAAAAGGCCC': [20, 20, 20, 22, 24, 26, 28, 32, 32, 32, 32, 32]
        }
        self.assertEqual(generated_reads, expected_reads)


if __name__ == '__main__':
    unittest.main()
