import subprocess
import unittest
from test_classes import TestGenerateVariantReads


class TestAlignmentFilter(TestGenerateVariantReads):

    def test_read_filter_se(self):
        # Add variant
        self.input.add_variant(
            name='V01', pos=5, ref='A', alt='T'
        )
        # Add unmapped single-end read
        self.input.add_alignment(
            name='R01', pos=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ',
            cigar='*', flag=4, mapq=30
        )
        # Add secondary alignment
        self.input.add_alignment(
            name='R02', pos=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ',
            cigar='10M', flag=256, mapq=30
        )
        # Add suplementary alingment
        self.input.add_alignment(
            name='R03', pos=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ',
            cigar='10M', flag=2048, mapq=30
        )
        # Add read just below mapping quality threshold
        self.input.add_alignment(
            name='R04', pos=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ',
            cigar='10M', flag=0, mapq=9
        )
        # Add read at mapping quality threshold
        self.input.add_alignment(
            name='R05', pos=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ',
            cigar='10M', flag=0, mapq=10
        )
        # Generate input files and run script
        arguments = [
            '--min_mapq', 10, '--max_vars', 4,
            '--max_seqs', 8, '--min_len', 10
        ]
        log = self.run_script(arguments)[2]
        # Extract log file and check alignment filter
        align_filter = log['alignment filter']
        self.assertEqual(align_filter['unmapped'], 1)
        self.assertEqual(align_filter['secondary'], 1)
        self.assertEqual(align_filter['supplementary'], 1)
        self.assertEqual(align_filter['low mapping quality'], 1)
        self.assertEqual(align_filter['passed'], 1)

    def test_read_filter_pe(self):
        # Add variant
        self.input.add_variant(
            name='V01', pos=5, ref='A', alt='T'
        )
        # Add paired-end read with one end unmapped
        self.input.add_alignment(
            name='R01', chrom='chr1', pos=1, cigar='10M', flag=73, mapq=30,
            rnext='*', pnext=0, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ'
        )
        self.input.add_alignment(
            name='R01', chrom='chr1', pos=1, cigar='*', flag=133, mapq=30,
            rnext='chr1', pnext=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ'
        )
        # Add secondary alignment for first read
        self.input.add_alignment(
            name='R01', chrom='chr1', pos=1, cigar='10M', flag=385, mapq=30,
            rnext='chr1', pnext=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ'
        )
        # Add suplementary alignment for first read
        self.input.add_alignment(
            name='R01', chrom='chr1', pos=1, cigar='10M', flag=2177, mapq=30,
            rnext='chr1', pnext=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ'
        )
        # Add read pair mapped to different chromosomes
        self.input.add_alignment(
            name='R02', chrom='chr1', pos=1, cigar='10M', flag=65, mapq=30,
            rnext='chr2', pnext=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ'
        )
        self.input.add_alignment(
            name='R02', chrom='chr2', pos=1, cigar='10M', flag=129, mapq=30,
            rnext='chr1', pnext=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ'
        )
        # Add improper read pair
        self.input.add_alignment(
            name='R03', chrom='chr1', pos=1, cigar='10M', flag=65, mapq=30,
            rnext='chr1', pnext=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ'
        )
        self.input.add_alignment(
            name='R03', chrom='chr1', pos=1, cigar='10M', flag=129, mapq=30,
            rnext='chr1', pnext=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ'
        )
        # Add read pair with one read just below mapping quality threshold
        self.input.add_alignment(
            name='R04', chrom='chr1', pos=1, cigar='10M', flag=67, mapq=30,
            rnext='chr1', pnext=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ'
        )
        self.input.add_alignment(
            name='R04', chrom='chr1', pos=1, cigar='10M', flag=131, mapq=9,
            rnext='chr1', pnext=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ'
        )
        # Add read pair with one read just at mapping quality threshold
        self.input.add_alignment(
            name='R05', chrom='chr1', pos=1, cigar='10M', flag=67, mapq=30,
            rnext='chr1', pnext=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ'
        )
        self.input.add_alignment(
            name='R05', chrom='chr1', pos=1, cigar='10M', flag=131, mapq=10,
            rnext='chr1', pnext=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ'
        )
        # Generate input files and run script
        arguments = [
            '--min_mapq', 10, '--max_vars', 4,
            '--max_seqs', 8, '--min_len', 10
        ]
        log = self.run_script(arguments)[2]
        # Extract log file and check alignment filter
        align_filter = log['alignment filter']
        self.assertEqual(align_filter['unmapped'], 1)
        self.assertEqual(align_filter['mate unmapped'], 1)
        self.assertEqual(align_filter['secondary'], 1)
        self.assertEqual(align_filter['supplementary'], 1)
        self.assertEqual(align_filter['different chromosomes'], 2)
        self.assertEqual(align_filter['improper pair'], 2)
        self.assertEqual(align_filter['low mapping quality'], 2)
        self.assertEqual(align_filter['passed'], 2)


if __name__ == '__main__':
    unittest.main(verbosity=2)
