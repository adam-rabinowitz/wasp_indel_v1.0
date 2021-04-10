import unittest
from test_classes import TestGenerateVariantReads


class TestSNV(TestGenerateVariantReads):

    def test01_single_snv(self):
        # Add variants
        self.add_variant(
            name='V01', chrom='chr1', pos=6, ref='A', alt='T'
        )
        # Add mapped single-end read
        self.add_alignment(
            name='R01', chrom='chr1', pos=1, sequence='AAAAAAAAAA',
            quality='ABCDEFGHIJ', cigar='10M', flag=67, mapq=30, rnext='chr1',
            pnext=4
        )
        self.add_alignment(
            name='R01', chrom='chr1', pos=4, sequence='AAAAAAAAAA',
            quality='HIJKLMNOPQ', cigar='10M', flag=131, mapq=30, rnext='chr1',
            pnext=4
        )
        # Run script and extract output
        arguments = [
            '--min_mapq', 10, '--max_vars', '3', '--max_seqs', 18,
            '--min_len', 8
        ]
        bam, fastq, log = self.run_script(arguments, paired=True)
        # Test output
        self.assertEqual(bam, {})
        self.assertEqual(
            fastq,
            {
                ('R01', '0-10-20'): set([
                    (('AAAAAAAAAA', 'ABCDEFGHIJ'), ('AAAAAAAAAA', 'HIJKLMNOPQ')),
            }
        )

if __name__ == '__main__':
    unittest.main(verbosity=2)