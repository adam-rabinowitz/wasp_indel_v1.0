import os
import pysam
import subprocess
import unittest
import yaml


class FlipAlleles(object):

    def create_paths(self):
        # Set script path
        test_dir = os.path.dirname(os.path.realpath(__file__))
        # Get script directory
        self.script = os.path.join(
            os.path.dirname(os.path.dirname(test_dir)),
            'generate_variant_reads.py'
        )
        # Create prefix
        self.prefix = os.path.join(
            test_dir, 'unit_test'
        )
        # Create input alignment paths
        self.sam = self.prefix + '.sam'
        self.bam = self.prefix + '.bam'
        self.bai = self.prefix + '.bai'
        # Create input vcf paths
        self.vcf = self.prefix + '.vcf'
        self.vcf_gzip = self.prefix + '.vcf.gz'
        self.tbi = self.prefix + '.vcf.gz.tbi'
        # Create output paths
        self.out_bam = self.prefix + '.no_variants.bam'
        self.out_fastq = self.prefix + '.allele_flipped.fq.gz'
        self.out_log = self.prefix + '.first_alignment_log.txt'
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

    def parse_bam(self):
        bam_dict = {}
        with pysam.AlignmentFile(self.out_bam) as bam:
            for read in bam:
                if read.query_name in bam_dict:
                    bam_dict[read.query_name] += 1
                else:
                    bam_dict[read.query_name] = 1
        return(bam_dict)

    def parse_fastq(self, paired):
        fastq_dict = {}
        cached_read = None
        with pysam.FastxFile(self.out_fastq) as fastq:
            for i, read in enumerate(fastq):
                # Process paired end read
                if paired:
                    # Store first read in pair and continue
                    if i % 2 == 0:
                        cached_read = read
                        continue
                    # Extract data for second read in pair
                    else:
                        read1, read2 = cached_read, read
                        assert(read1.name == read2.name)
                        chached_read = None
                        # Extract read data
                        name, location = read1.name.split('.')[0:2]
                        read_id = (name, location)
                        read_data = (
                            (read1.sequence, read1.quality),
                            (read2.sequence, read2.quality)
                        )
                # Process signle end read
                else:
                    # Extract data for unpaired read
                    name, location = read.name.split('.')[0:2]
                    read_id = (name, location)
                    read_data = (read.sequence, read.quality)
                # Store read data
                if read_id not in fastq_dict:
                    fastq_dict[read_id] = set()
                fastq_dict[read_id].add(read_data)
        return(fastq_dict)

    def parse_log(self):
        with open(self.out_log, 'rt') as log:
            log_dict = yaml.safe_load(log)
        return(log_dict)

    def parse_all(self, paired):
        bam_dict = self.parse_bam()
        fastq_dict = self.parse_fastq(paired=paired)
        log_dict = self.parse_log()
        return(bam_dict, fastq_dict, log_dict)

    def run_script(self, arguments, paired=False):
        # Prepare input
        self.prepare_input()
        # Create command
        arguments = list(map(str, arguments))
        command = [
            'python', self.script, '--bam', self.bam, '--vcf',
            self.vcf_gzip, '--out_prefix', self.prefix
        ] + arguments
        # Run command
        try:
            subprocess.check_call(command)
        except subprocess.CalledProcessError:
            print('\n' + ' '.join(command))
            raise
        # Parse and return files
        return(self.parse_all(paired=paired))

    def delete(self):
        for path in (
            self.sam, self.bam, self.bai, self.vcf, self.vcf_gzip, self.tbi,
            self.out_bam, self.out_fastq, self.out_log
        ):
            if os.path.isfile(path):
                os.remove(path)


class TestGenerateVariantReads(FlipAlleles, unittest.TestCase):

    def setUp(self):
        self.create_paths()

    def tearDown(self):
        self.delete()
