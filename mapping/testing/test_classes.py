import os
import pysam
import unittest
import yaml


class InputFiles(object):

    def __init__(self, prefix):
        # Store prefix
        if prefix.endswith('.'):
            prefix = prefix[:-1]
        self.prefix = prefix
        # Create alignment paths
        self.sam = prefix + '.sam'
        self.bam = prefix + '.bam'
        self.bai = prefix + '.bai'
        # Create vcf paths
        self.vcf = prefix + '.vcf'
        self.vcf_gzip = prefix + '.vcf.gz'
        self.tbi = prefix + '.vcf.gz.tbi'
        # Create sam file with header
        with open(self.sam, 'wt') as outfile:
            outfile.write(
                '@HD\tVN:1.6\tSO:coordinate\n'
                '@SQ\tSN:chr1\tLN:100\n'
                '@SQ\tSN:chr2\tLN:100\n'
            )
        # Create vcf file with header
        with open(self.vcf, 'wt') as outfile:
            outfile.write(
                '##fileformat=VCFv4.2\n'
                '##FORMAT=<ID=GT,Number=1,Type=String>\n'
                '##FORMAT=<ID=PL,Number=G,Type=Integer>\n'
                '##contig=<ID=chr1,length=100>\n'
                '##contig=<ID=chr2,length=100>\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
                'S001\tS002\n'
            )

    def add_alignment(
        self, name, pos, sequence, quality, cigar, chrom='chr1', flag=0,
        mapq=60, rnext='*', pnext=0, tlen=0
    ):
        # Create string
        alignment_string = (
            '{name}\t{flag}\t{chrom}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t'
            '{pnext}\t{tlen}\t{sequence}\t{quality}\n'
        ).format(
            name=name, flag=flag, chrom=chrom, pos=pos, mapq=mapq,
            cigar=cigar, rnext=rnext, pnext=pnext, tlen=tlen,
            sequence=sequence, quality=quality
        )
        # Write string to file
        with open(self.sam, 'at') as outfile:
            outfile.write(alignment_string)

    def add_variant(
        self, name, pos, ref, alt, chrom='chr1', quality=100, filter='.',
        info='.', S001='0/0:0,30,300', S002='0/0:0,30,300'
    ):
        # Create string
        variant_string = (
            '{chrom}\t{pos}\t{name}\t{ref}\t{alt}\t{quality}\t{filter}\t'
            '{info}\t{S001}\t{S002}\n'
        ).format(
            chrom=chrom, pos=pos, name=name, ref=ref, alt=alt,
            quality=quality, filter=filter, info=info, S001=S001, S002=S002
        )
        # Write string to file
        with open(self.vcf, 'at') as outfile:
            outfile.write(variant_string)

    def prepare_input(self):
        # Compress and index SAM file
        pysam.sort('-o', self.bam, '-O', 'bam', self.sam)
        pysam.index(self.bam, self.bai)
        # Compress and index VCF file
        pysam.tabix_compress(self.vcf, self.vcf_gzip, force=True)
        pysam.tabix_index(
            self.vcf_gzip, force=True, preset='vcf', index=self.tbi
        )
        # Return indexed files
        return(self.bam, self.vcf_gzip)

    def delete(self):
        for path in (
            self.sam, self.bam, self.bai, self.vcf, self.vcf_gzip, self.tbi
        ):
            if os.path.isfile(path):
                os.remove(path)


class OutputFiles(object):

    def __init__(self, prefix):
        # Store prefix
        if prefix.endswith('.'):
            prefix = prefix[:-1]
        self.prefix = prefix
        # Create paths
        self.bam = prefix + '.no_variants.bam'
        self.fastq = prefix + '.allele_flipped.fq.gz'
        self.log = prefix + '.first_alignment_log.txt'

    def parse_bam(self):
        bam_dict = {}
        with pysam.AlignmentFile(self.bam) as bam:
            for read in bam:
                if read.query_name in bam_dict:
                    bam_dict[read.query_name] += 1
                else:
                    bam_dict[read.query_name] = 1
        return(bam_dict)

    def parse_fastq(self):
        fastq_dict = {}
        with pysam.FastxFile(self.fastq) as fastq:
            for read in fastq:
                name, location, total, count = read.name.split('.')
                read_id = (name, location)
                if read_id not in fastq_dict:
                    fastq_dict[read_id] = [None] * int(total)
                fastq_dict[read_id][int(count)] = (read.sequence, read.quality)
        return(fastq_dict)

    def parse_log(self):
        with open(self.log, 'rt') as log:
            log_dict = yaml.safe_load(log)
        return(log_dict)

    def parse_all(self):
        bam_dict = self.parse_bam()
        fastq_dict = self.parse_fastq()
        log_dict = self.parse_log()
        return(bam_dict, fastq_dict, log_dict)

    def delete(self):
        os.remove(self.bam)
        os.remove(self.fastq)
        os.remove(self.log)


class TestGenerateVariantReads(unittest.TestCase):

    def setUp(self):
        self.script = '../generate_variant_reads.py'
        self.prefix = 'unit_test'
        self.input = InputFiles(self.prefix)
        self.output = OutputFiles(self.prefix)

    def tearDown(self):
        self.input.delete()
        self.output.delete()
