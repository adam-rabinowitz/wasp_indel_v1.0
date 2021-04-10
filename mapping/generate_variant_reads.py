import argparse
import collections
import gzip
import pysam
# Import classes from modules
from mapping.bams import FirstBamGenerator
from mapping.flipvar import FlipVar
from mapping.vartree import VarTree


class Log(object):

    def __init__(self, path):
        self.path = path
        self.counts = collections.Counter()

    def write_counts(self):
        # Get alignment counts log
        alignment_log = (
            'total reads: {total_align}\n'
            'alignment filter:\n'
            '  secondary: {secondary_align}\n'
            '  supplementary: {supplementary_align} \n'
            '  unmapped: {unmapped_align} \n'
            '  mate unmapped: {mate_unmapped_align}\n'
            '  different chromosomes: {different_chrom_align}\n'
            '  improper pair: {improper_pair_align}\n'
            '  low mapping quality: {low_mapq_align}\n'
            '  abnormal cigar: {abnormal_cigar_align}\n'
            '  passed: {passed_align}\n'
        ).format(
            total_align=self.counts['total_align'],
            secondary_align=self.counts['secondary_align'],
            supplementary_align=self.counts['supplementary_align'],
            unmapped_align=self.counts['unmapped_align'],
            mate_unmapped_align=self.counts['mate_unampped_align'],
            different_chrom_align=self.counts['different_chrom_align'],
            improper_pair_align=self.counts['improper_pair_align'],
            low_mapq_align=self.counts['low_mapq_align'],
            abnormal_cigar_align=self.counts['abnormal_cigar_align'],
            passed_align=self.counts['passed_align']
        )
        # Get variant counts log
        variant_count_log = (
            'read variant types:\n'
            '  none: {no_variants}\n'
            '  only biallelic snv: {biallelic_snv_variants}\n'
            '  other: {other_variants}\n'
        ).format(
            no_variants=self.counts['no_variants'],
            biallelic_snv_variants=self.counts['biallelic_snv_variants'],
            other_variants=self.counts['other_variants']
        )
        # Get allele count log
        allele_count_log = (
            'read allele counts:\n'
            '  reference: {ref_allele}\n'
            '  alternative: {alt_allele}\n'
            '  other: {other_allele}\n'
        ).format(
            ref_allele=self.counts['ref_allele'],
            alt_allele=self.counts['alt_allele'],
            other_allele=self.counts['other_allele']
        )
        # Get variant filter log
        variant_filter_log = (
            'variant filter:\n'
            '  unwanted variants: {unwanted_variants}\n'
            '  excess variants: {excess_variants}\n'
            '  overlapping variants: {overlapping_variants}\n'
            '  conflicting alleles: {conflicting_alleles}\n'
            '  incomplete haplotype: {incomplete_haplotype}\n'
            '  excess reads: {excess_reads}\n'
            '  truncated reads: {too_short}\n'
            '  to remap: {to_remap}\n'
        ).format(
            unwanted_variants=self.counts['unwanted_variants'],
            excess_variants=self.counts['excess_variants'],
            overlapping_variants=self.counts['overlapping_variants'],
            conflicting_alleles=self.counts['conflicting_alleles'],
            incomplete_haplotype=self.counts['incomplete_haplotype'],
            excess_reads=self.counts['excess_reads'],
            too_short=self.counts['too_short'],
            to_remap=self.counts['to_remap']
        )
        # Write logs to file
        with open(self.path, 'wt') as outfile:
            outfile.write(alignment_log)
            outfile.write(variant_count_log)
            outfile.write(allele_count_log)
            outfile.write(variant_filter_log)


class ProcessAlignments(object):

    def __init__(
        self, bam, vcf, out_prefix, min_mapq, max_vars, max_seqs, min_len,
        only_biallelic_snv, samples, check_phase
    ):
        # Store input parameters
        self.inbam_path = bam
        self.vcf_path = vcf
        self.out_prefix = out_prefix
        self.min_mapq = min_mapq
        self.max_vars = max_vars
        self.max_seqs = max_seqs
        self.min_len = min_len
        self.only_biallelic_snv = only_biallelic_snv
        self.samples = samples
        self.check_phase = check_phase
        # Create output paths
        self.outbam_path = self.out_prefix + '.no_variants.bam'
        self.fastq_path = self.out_prefix + '.allele_flipped.fq.gz'
        self.log_path = self.out_prefix + '.first_alignment_log.txt'
        # Initialise objects and open files
        self.bam_generator = FirstBamGenerator(
            self.inbam_path, min_mapq=self.min_mapq
        )
        self.vartree = VarTree(
            self.vcf_path, samples=self.samples, check_phase=self.check_phase
        )
        self.outbam = pysam.AlignmentFile(
            self.outbam_path, mode='wb', template=self.bam_generator.bam
        )
        self.fastq = gzip.open(self.fastq_path, 'wt')
        # Create object to log processing
        self.log = Log(self.log_path)

    def close(
        self
    ):
        # Close all open objects and files
        self.bam_generator.close()
        self.outbam.close()
        self.fastq.close()

    def filter_reads(
        self
    ):
        # Get variants for each chromosome
        for chrom_stats in self.bam_generator.index_statistics:
            # Get chromosome name and total mapped reads
            chromosome = chrom_stats.contig
            total_reads = chrom_stats.total
            # Skip chromosomes without mapped reads
            if total_reads == 0:
                continue
            # Extract variants for chromosome
            self.vartree.read_vcf(chromosome)
            # Loop through reads on chromosome
            for reads in self.bam_generator.get_reads(chromosome):
                read_no = len(reads)
                # Get variants for each read and create object to flip reads
                variants = [
                    self.vartree.get_read_variants(read, partial=True)
                    for read in reads
                ]
                flipvar = FlipVar(reads, variants)
                # Save reads with no variants
                if not flipvar.unique_variants:
                    self.log.counts['no_variants'] += read_no
                    for read in reads:
                        self.outbam.write(read)
                    continue
                # Update allele counts
                self.log.counts.update(flipvar.count_variant_alleles())
                # Count variant types
                biallelic_snv = all(
                    [v.is_biallelic_snv() for v in flipvar.unique_variants]
                )
                if biallelic_snv:
                    self.log.counts['biallelic_snv_variants'] += read_no
                else:
                    self.log.counts['other_variants'] += read_no
                # Count reads with unwanted variants and continue
                if not biallelic_snv and self.only_biallelic_snv:
                    self.log.counts['unwanted_variants'] += read_no
                    continue
                # Count reads with excess variants and continue
                if flipvar.max_variant_count() > self.max_vars:
                    self.log.counts['excess_variants'] += read_no
                    continue
                # Count reads containing overlapping variants and continue
                if flipvar.variants_overlap():
                    self.log.counts['overlapping_variants'] += read_no
                    continue
                # Count paired reads with conflicting alleles and continue
                if flipvar.conflicting_alleles():
                    self.log.counts['conflicting_alleles'] += read_no
                    continue
                # Get flipped read haplotypes
                if self.samples and flipvar.incomplete_haplotype():
                    self.log.counts['incomplete_haplotype'] += read_no
                    continue
                # Count flipped read haplotypes
                flipvar.get_haplotypes(self.samples)
                if flipvar.count_haplotypes() > self.max_seqs:
                    self.log.counts['excess_reads'] += read_no
                    continue
                # Generate flipped reads and check lengths
                flipvar.generate_flipped_reads()
                if flipvar.min_length() < self.min_len:
                    self.log.counts['too_short'] += read_no
                    continue
                # Generate fastq and write to file
                fastq_entry = flipvar.generate_fastq()
                self.fastq.write(fastq_entry)
                self.log.counts['to_remap'] += read_no
        # Get alignment counts
        alignment_counts = self.bam_generator.counts
        alignment_counts['total'] += self.bam_generator.nocoordinate
        alignment_counts['unmapped'] += self.bam_generator.nocoordinate
        # Add alignment counts to counter and write to file
        self.log.counts.update(alignment_counts)
        self.log.write_counts()


# Run script
if __name__ == '__main__':
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="Looks for variants (SNPs & indels) overlapping reads. "
        "Three output files are generated with the provided prefix and the "
        "following suffixes. 1) '.invariant.bam' - BAM file containing reads "
        "not overlapping variants. 2) '.remap.fq.gz' - a gzipped FASTQ file "
        "containing the original and variant flipped reads for realignment. "
        "If the input BAM file contains paired end reads then the pairs are "
        "interlevead in the FASTQ file. 3) '.variant_log.txt' text file "
        "containing read processing metrics."
    )
    parser.add_argument(
        "--bam", required=True, help=(
            "Coordinate-sorted and indexed input BAM file."
        )
    )
    parser.add_argument(
        "--vcf", required=True, help=(
            "Coordinate sorted and tabix indexed VCF file."
        )
    )
    parser.add_argument(
        "--out_prefix", required=True, help=(
            "Prefix of output files."
        )
    )
    parser.add_argument(
        "--min_mapq", type=int, default=0, help=(
            "Minimum read mapping quality (default=0)."
        )
    )
    parser.add_argument(
        "--max_vars", type=int, default=6, help=(
            "The maximum number of variants allowed to overlap a read before "
            "discarding the read/read-pair (default=6)."
        )
    )
    parser.add_argument(
        "--max_seqs", type=int, default=64, help=(
            "The maximum number of allele flipped reads/read-pairs to "
            " consider remapping. Reads with more allelic combinations "
            "than the specified value are discarded (defult=64)."
        )
    )
    parser.add_argument(
        "--min_len", type=int, default=20, help=(
            "The minimum length of allele flipped reads for remapping "
            "(default=20)"
        )
    )
    parser.add_argument(
        "--biallelic_snv_only", action='store_true', help=(
            "Discard reads overlapping non biallelic-snv variants "
            "(default=False)."
        )
    )
    parser.add_argument(
        "--samples", nargs='+', help=(
            "Sample names for generation of phased allele flipped reads."
        )
    )
    parser.add_argument(
        "--assume_phased", action='store_true', help=(
            "Assume sample genotypes in VCF are phased Only considered "
            "when '--samples' are specified (default=False)."
        )
    )
    args = parser.parse_args()
    # Run program
    check_phase = not args.assume_phased
    process_alignments = ProcessAlignments(
        bam=args.bam, vcf=args.vcf, out_prefix=args.out_prefix,
        min_mapq=args.min_mapq, max_vars=args.max_vars, max_seqs=args.max_seqs,
        min_len=args.min_len, only_biallelic_snv=args.biallelic_snv_only,
        samples=args.samples, check_phase=check_phase
    )
    process_alignments.filter_reads()
    process_alignments.close()
