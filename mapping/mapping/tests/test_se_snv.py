import unittest
from test_classes import TestGenerateVariantReads


class TestSNV(TestGenerateVariantReads):

    def test01_multiple_snv(self):
        # Add variants
        self.add_variant(
            name='V01', chrom='chr1', pos=3, ref='A', alt='T'
        )
        self.add_variant(
            name='V02', chrom='chr1', pos=8, ref='T', alt='G'
        )
        self.add_variant(
            name='V03', chrom='chr2', pos=11, ref='T', alt='A'
        )
        self.add_variant(
            name='V04', chrom='chr2', pos=19, ref='C', alt='G'
        )
        self.add_variant(
            name='V05', chrom='chr2', pos=20, ref='A', alt='G,C'
        )
        # Add mapped single-end read
        self.add_alignment(
            name='R01', chrom='chr1', pos=1, sequence='AAAAAAAAAA',
            quality='ABCDEFGHIJ', cigar='10M', flag=0, mapq=30
        )
        self.add_alignment(
            name='R02', chrom='chr1', pos=3, sequence='AAAAAAAAAA',
            quality='ABCDEFGHIJ', cigar='10M', flag=16, mapq=30
        )
        self.add_alignment(
            name='R03', chrom='chr2', pos=1, sequence='AAAAAAAAAA',
            quality='ABCDEFGHIJ', cigar='10M', flag=0, mapq=30
        )
        self.add_alignment(
            name='R04', chrom='chr2', pos=11, sequence='AAAAAAAAAA',
            quality='ABCDEFGHIJ', cigar='10M', flag=0, mapq=30
        )
        self.add_alignment(
            name='R05', chrom='chr2', pos=21, sequence='AAAAAAAAAA',
            quality='ABCDEFGHIJ', cigar='10M', flag=0, mapq=30
        )
        # Generate input files and run script
        arguments = [
            '--min_mapq', 10, '--max_vars', 3,
            '--max_seqs', 18, '--min_len', 10
        ]
        bam, fastq, log = self.run_script(arguments)
        # Extract log file and check alignment filter
        self.assertEqual(bam, {'R03': 1, 'R05': 1})
        self.assertEqual(
            fastq,
            {
                ('R01', '0-0-10'): set([
                    ('AAAAAAAAAA', 'ABCDEFGHIJ'),
                    ('AAAAAAATAA', 'ABCDEFGHIJ'),
                    ('AAAAAAAGAA', 'ABCDEFGHIJ'),
                    ('AATAAAAAAA', 'ABCDEFGHIJ'),
                    ('AATAAAATAA', 'ABCDEFGHIJ'),
                    ('AATAAAAGAA', 'ABCDEFGHIJ'),
                ]),
                ('R02', '0-2-12'): set([
                    ('TTTTTTTTTT', 'JIHGFEDCBA'),
                    ('TTTTATTTTT', 'JIHGFEDCBA'),
                    ('TTTTCTTTTT', 'JIHGFEDCBA'),
                    ('TTTTTTTTTA', 'JIHGFEDCBA'),
                    ('TTTTATTTTA', 'JIHGFEDCBA'),
                    ('TTTTCTTTTA', 'JIHGFEDCBA')
                ]),
                ('R04', '1-10-20'): set([
                    ('AAAAAAAAAA', 'ABCDEFGHIJ'),
                    ('AAAAAAAAAG', 'ABCDEFGHIJ'),
                    ('AAAAAAAAAC', 'ABCDEFGHIJ'),
                    ('AAAAAAAACA', 'ABCDEFGHIJ'),
                    ('AAAAAAAACG', 'ABCDEFGHIJ'),
                    ('AAAAAAAACC', 'ABCDEFGHIJ'),
                    ('AAAAAAAAGA', 'ABCDEFGHIJ'),
                    ('AAAAAAAAGG', 'ABCDEFGHIJ'),
                    ('AAAAAAAAGC', 'ABCDEFGHIJ'),
                    ('TAAAAAAAAA', 'ABCDEFGHIJ'),
                    ('TAAAAAAAAG', 'ABCDEFGHIJ'),
                    ('TAAAAAAAAC', 'ABCDEFGHIJ'),
                    ('TAAAAAAACA', 'ABCDEFGHIJ'),
                    ('TAAAAAAACG', 'ABCDEFGHIJ'),
                    ('TAAAAAAACC', 'ABCDEFGHIJ'),
                    ('TAAAAAAAGA', 'ABCDEFGHIJ'),
                    ('TAAAAAAAGG', 'ABCDEFGHIJ'),
                    ('TAAAAAAAGC', 'ABCDEFGHIJ')
                ])
            }
        )
        self.assertEqual(log['read variant types']['none'], 2)
        self.assertEqual(log['read variant types']['only biallelic snv'], 2)
        self.assertEqual(log['read variant types']['other'], 1)
        self.assertEqual(log['read allele counts']['reference'], 3)
        self.assertEqual(log['read allele counts']['alternative'], 1)
        self.assertEqual(log['read allele counts']['other'], 3)

    def test02_snv_over_indels(self):
        # Add variants
        self.add_variant(
            name='V01', pos=13, ref='T', alt='A'
        )
        self.add_variant(
            name='V02', pos=18, ref='C', alt='G'
        )
        # Add mapped single-end read
        self.add_alignment(
            name='R01', pos=11, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ',
            cigar='10M', flag=0, mapq=30
        )
        self.add_alignment(
            name='R02', pos=11, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ',
            cigar='2M1D5M2I1M', flag=0, mapq=30
        )
        # Generate input files and run script
        arguments = [
            '--min_mapq', 10, '--max_vars', 3,
            '--max_seqs', 18, '--min_len', 8
        ]
        bam, fastq, log = self.run_script(arguments)
        # Extract log file and check alignment filter
        self.assertEqual(bam, {})
        self.assertEqual(
            fastq,
            {
                ('R01', '0-10-20'): set([
                    ('AAAAAAAAAA', 'ABCDEFGHIJ'),
                    ('AAAAAAACAA', 'ABCDEFGHIJ'),
                    ('AAAAAAAGAA', 'ABCDEFGHIJ'),
                    ('AATAAAAAAA', 'ABCDEFGHIJ'),
                    ('AATAAAACAA', 'ABCDEFGHIJ'),
                    ('AATAAAAGAA', 'ABCDEFGHIJ')
                ]),
                ('R02', '0-10-19'): set([
                    ('AAAAAAAAAA', 'ABCDEFGHIJ'),
                    ('AAAAAACA', 'ABCDEFHJ'),
                    ('AAAAAAGA', 'ABCDEFHJ')
                ])
            }
        )
        self.assertEqual(log['read variant types']['only biallelic snv'], 2)
        self.assertEqual(log['read allele counts']['reference'], 0)
        self.assertEqual(log['read allele counts']['alternative'], 1)
        self.assertEqual(log['read allele counts']['other'], 2)

    def test03_excess_seqs(self):
        # Add variants
        self.add_variant(
            name='V01', pos=11, ref='T', alt='A'
        )
        self.add_variant(
            name='V02', pos=19, ref='C', alt='G'
        )
        self.add_variant(
            name='V03', pos=20, ref='A', alt='G,C'
        )
        # Add mapped single-end read
        self.add_alignment(
            name='R01', pos=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ',
            cigar='10M', flag=0, mapq=30
        )
        self.add_alignment(
            name='R02', pos=11, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ',
            cigar='10M', flag=0, mapq=30
        )
        self.add_alignment(
            name='R03', pos=21, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ',
            cigar='10M', flag=0, mapq=30
        )
        # Generate input files and run script
        arguments = [
            '--min_mapq', 10, '--max_vars', 3,
            '--max_seqs', 17, '--min_len', 10
        ]
        bam, fastq, log = self.run_script(arguments)
        # Extract log file and check alignment filter
        self.assertEqual(bam, {'R01': 1, 'R03': 1})
        self.assertEqual(fastq, {})
        self.assertEqual(log['read variant types']['none'], 2)
        self.assertEqual(log['read variant types']['other'], 1)
        self.assertEqual(log['variant filter']['excess reads'], 1)
        self.assertEqual(log['read allele counts']['reference'], 1)
        self.assertEqual(log['read allele counts']['alternative'], 1)
        self.assertEqual(log['read allele counts']['other'], 1)

    def test04_excess_vars(self):
        # Add variants
        self.add_variant(
            name='V01', pos=11, ref='T', alt='A'
        )
        self.add_variant(
            name='V02', pos=19, ref='C', alt='G'
        )
        self.add_variant(
            name='V03', pos=20, ref='A', alt='G,C'
        )
        # Add mapped single-end read
        self.add_alignment(
            name='R01', pos=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ',
            cigar='10M', flag=0, mapq=30
        )
        self.add_alignment(
            name='R02', pos=11, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ',
            cigar='10M', flag=0, mapq=30
        )
        self.add_alignment(
            name='R03', pos=21, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ',
            cigar='10M', flag=0, mapq=30
        )
        # Generate input files and run script
        arguments = [
            '--min_mapq', 10, '--max_vars', 2,
            '--max_seqs', 18, '--min_len', 10
        ]
        bam, fastq, log = self.run_script(arguments)
        # Extract log file and check alignment filter
        self.assertEqual(bam, {'R01': 1, 'R03': 1})
        self.assertEqual(fastq, {})
        self.assertEqual(log['read variant types']['none'], 2)
        self.assertEqual(log['read variant types']['other'], 1)
        self.assertEqual(log['variant filter']['excess variants'], 1)
        self.assertEqual(log['read allele counts']['reference'], 1)
        self.assertEqual(log['read allele counts']['alternative'], 1)
        self.assertEqual(log['read allele counts']['other'], 1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
