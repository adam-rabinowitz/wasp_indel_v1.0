import unittest
from test_classes import TestGenerateVariantReads


class TestSNV(TestGenerateVariantReads):

    def test01_multiple_snv_1(self):
        # Add variants
        self.input.add_variant(
            name='V01', pos=3, ref='AAA', alt='T,TTTT'
        )
        self.input.add_variant(
            name='V02', pos=8, ref='AT', alt='GAT'
        )
        # Add mapped single-end read
        self.input.add_alignment(
            name='R01', pos=1, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ',
            cigar='10M', flag=0, mapq=30
        )
        self.input.add_alignment(
            name='R02', pos=3, sequence='AAAAAAAAAA', quality='ABCDEFGHIJ',
            cigar='10M', flag=16, mapq=30
        )
        # Run script and extract output
        arguments = [
            '--min_mapq', 10, '--max_vars', '3', '--max_seqs', 18,
            '--min_len', 8
        ]
        bam, fastq, log = self.run_script(arguments)
        # Extract
        self.assertEqual(bam, {})
        self.assertEqual(
            fastq,
            {
                ('R01', '0-0-10'): set([
                    ('AAAAAAAAAA', 'ABCDEFGHIJ'),
                    ('AAAAAAAATA', 'ABCDEFGHIJ'),
                    ('AAAAAAAGATA', 'ABCDEFGHHHJ'),
                    ('AATAAAAA', 'ABDFGHIJ'),
                    ('AATAAATA', 'ABDFGHIJ'),
                    ('AATAAGATA', 'ABDFGHHHJ'),
                    ('AATTTTAAAAA', 'ABDDDDFGHIJ'),
                    ('AATTTTAAATA', 'ABDDDDFGHIJ'),
                    ('AATTTTAAGATA', 'ABDDDDFGHHHJ')
                ]),
                ('R02', '0-2-12'): set([
                    ('TTTTTTTTTT', 'JIHGFEDCBA'),
                    ('TTTATTTTTT', 'JIHGFEDCBA'),
                    ('TTTATCTTTTT', 'JIHFFFEDCBA'),
                    ('TTTTTTTA', 'JIHGFEDB'),
                    ('TTTATTTA', 'JIHGFEDB'),
                    ('TTTATCTTA', 'JIHFFFEDB'),
                    ('TTTTTTTAAAA', 'JIHGFEDBBBB'),
                    ('TTTATTTAAAA', 'JIHGFEDBBBB'),
                    ('TTTATCTTAAAA', 'JIHFFFEDBBBB')
                ])
            }
        )
        self.assertEqual(log['read variant types']['mixed'], 2)
        self.assertEqual(log['read allele counts']['reference'], 2)
        self.assertEqual(log['read allele counts']['alternative'], 0)
        self.assertEqual(log['read allele counts']['other'], 2)


if __name__ == '__main__':
    unittest.main(verbosity=2)
