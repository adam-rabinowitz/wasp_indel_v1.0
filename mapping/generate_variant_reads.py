import argparse
import collections
import gzip
import pysam
# Import classes from modules
from flipvar import FlipVar, HaplotypeError
from vartree import VarTree, CigarError


class BamGenerator():

    def __init__(self, bam, min_mapq):
        # Store parameters and open BAM file
        self.bam_path = bam
        self.bam = pysam.AlignmentFile(self.bam_path)
        self.min_mapq = min_mapq
        # Check pairing
        paired = [read.is_paired for read in self.bam.head(100)]
        if any(paired):
            assert(all(paired))
            self.paired = True
        else:
            self.paired = False
        # Get genome and alignment metrics
        self.index_statistics = self.bam.get_index_statistics()
        self.nocoordinate = self.bam.nocoordinate
        # Create counter
        self.counter = {
            'total': 0,
            'secondary': 0,
            'supplementary': 0,
            'unmapped': 0,
            'mate_unmapped': 0,
            'different_chromosomes': 0,
            'improper_pair': 0,
            'low_mapq': 0,
            'passed': 0
        }

    def close(self):
        self.bam.close()

    def get_reads(self, chromosome):
        # Create cache to store unfound reads
        read_pair_cache = {}
        # Loop though reads on chromosome
        for read in self.bam.fetch(contig=chromosome):
            # Count total reads
            self.counter['total'] += 1
            # Count and skip secondary reads
            if read.is_secondary:
                self.counter['secondary'] += 1
                continue
            # Count and skip supplementary reads
            if read.is_supplementary:
                self.counter['supplementary'] += 1
                continue
            # Count and skip unmapped reads
            if read.is_unmapped:
                self.counter['unmapped'] += 1
                continue
            # Count and skip mate unmapped reads
            if self.paired and read.mate_is_unmapped:
                self.counter['mate_unmapped'] += 1
                continue
            # Count and skip mates mapped to different chromosome
            if self.paired and read.reference_id != read.next_reference_id:
                self.counter['different_chromosomes'] += 1
                continue
            # Count and skip imporperly paired mates
            if self.paired and not read.is_proper_pair:
                self.counter['improper_pair'] += 1
                continue
            # Process reads passing initial filters
            if self.paired:
                # Get paired reads if mate has been cached...
                if read.query_name in read_pair_cache:
                    cached_read = read_pair_cache.pop(read.query_name)
                    if read.is_read1:
                        assert(cached_read.is_read2)
                        mapped_reads = [read, cached_read]
                    else:
                        assert(cached_read.is_read1)
                        assert(read.is_read2)
                        mapped_reads = [cached_read, read]
                # or store first in pair
                else:
                    read_pair_cache[read.query_name] = read
                    mapped_reads = []
            # Process single end reads
            else:
                assert(not read.is_paired)
                mapped_reads = [read]
            # Further filter properly mapped reads
            if mapped_reads:
                # Count and skip reads with low mapping quality
                low_mapq = [
                    r.mapping_quality < self.min_mapq for r in mapped_reads
                ]
                if any(low_mapq):
                    self.counter['low_mapq'] += len(mapped_reads)
                else:
                    self.counter['passed'] += len(mapped_reads)
                    yield(mapped_reads)
        # Check read pair cache is empty
        assert(len(read_pair_cache) == 0)


class CreateLog(object):

    def __init__(self, path):
        self.path = path

    def write_counts(
        self, alignment_counts, variant_counts
    ):
        # Get alignment counts log
        alignment_log = (
            'total reads: {total}\n'
            'alignment filter:\n'
            '  secondary: {secondary}\n'
            '  supplementary: {supplementary} \n'
            '  unmapped: {unmapped} \n'
            '  mate unmapped: {mate_unmapped}\n'
            '  different chromosomes: {different_chromosomes}\n'
            '  improper pair: {improper_pair}\n'
            '  low mapping quality: {low_mapq}\n'
            '  passed: {passed}\n'
        ).format(**alignment_counts)
        # Get variant counts log
        variant_count_log = (
            'read variant types:\n'
            '  none: {variants_absent}\n'
            '  unknown: {abnormal_alignment}\n'
            '  biallelic snv: {biallelic_snv}\n'
            '  biallelic: {biallelic}\n'
            '  snv: {snv}\n'
            '  mixed: {mixed}\n'
        ).format(**variant_counts)
        # Get allele count log
        allele_count_log = (
            'read allele counts:\n'
            '  reference: {ref}\n'
            '  alternative: {alt}\n'
            '  other: {other}\n'
        ).format(**variant_counts['alleles'])
        # Get variant filter log
        variant_filter_log = (
            'variant filter:\n'
            '  unwanted variants: {unwanted_variants}\n'
            '  excess variants: {excess_variants}\n'
            '  overlapping variants: {overlapping_variants}\n'
            '  conflicting alleles: {conflicting_alleles}\n'
            '  incomplete haplotype: {incomplete_haplotype}\n'
            '  truncated reads: {too_short}\n'
            '  excess reads: {excess_reads}\n'
            '  to remap: {to_remap}\n'
        ).format(**variant_counts)
        # Write logs to file
        with open(self.path, 'wt') as outfile:
            outfile.write(alignment_log)
            outfile.write(variant_count_log)
            outfile.write(allele_count_log)
            outfile.write(variant_filter_log)


class ProcessAlignments(object):

    def __init__(
        self, bam, vcf, out_prefix, min_mapq, max_vars, max_seqs, min_len,
        only_snv, only_biallelic, samples, check_phase
    ):
        # Store input parameters
        self.inbam_path = bam
        self.vcf_path = vcf
        self.out_prefix = out_prefix
        self.min_mapq = min_mapq
        self.max_vars = max_vars
        self.max_seqs = max_seqs
        self.min_len = min_len
        self.only_snv = only_snv
        self.only_biallelic = only_biallelic
        self.samples = samples
        self.check_phase = check_phase
        # Create output paths
        self.outbam_path = self.out_prefix + '.no_variants.bam'
        self.fastq_path = self.out_prefix + '.allele_flipped.fq.gz'
        self.log_path = self.out_prefix + '.first_alignment_log.txt'
        # Initialise objects and open files
        self.bam_generator = BamGenerator(
            self.inbam_path, min_mapq=self.min_mapq
        )
        self.vartree = VarTree(
            self.vcf_path, samples=self.samples, check_phase=self.check_phase
        )
        self.outbam = pysam.AlignmentFile(
            self.outbam_path, mode='wb', template=self.bam_generator.bam
        )
        self.fastq = gzip.open(self.fastq_path, 'wt')
        self.log = CreateLog(self.log_path)
        # Create counter
        self.counter = {
            'abnormal_alignment': 0,
            'variants_absent': 0,
            'biallelic_snv': 0,
            'biallelic': 0,
            'snv': 0,
            'mixed': 0,
            'unwanted_variants': 0,
            'excess_variants': 0,
            'overlapping_variants': 0,
            'conflicting_alleles': 0,
            'incomplete_haplotype': 0,
            'too_short': 0,
            'excess_reads': 0,
            'to_remap': 0,
            'alleles': collections.Counter(
                {'ref': 0, 'alt': 0, 'other': 0}
            )
        }

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
                # Get variants for each read...
                try:
                    variants = [
                        self.vartree.get_read_variants(read, partial=True) for
                        read in reads
                    ]
                # ...or count abnormal alignments and continue
                except CigarError:
                    self.counter['abnormal_alignment'] += read_no
                    continue
                # Create object to flip reads
                flipvar = FlipVar(reads, variants)
                unique_variants = flipvar.get_unique_variants()
                # Or save reads with no variants
                if not unique_variants:
                    self.counter['variants_absent'] += read_no
                    for read in reads:
                        self.outbam.write(read)
                    continue
                # Update allele counts
                self.counter['alleles'].update(flipvar.count_variant_alleles())
                # Count variant types...
                snv = all([v.is_snv() for v in unique_variants])
                biallelic = all([v.is_biallelic() for v in unique_variants])
                if snv:
                    if biallelic:
                        self.counter['biallelic_snv'] += read_no
                    else:
                        self.counter['snv'] += read_no
                else:
                    if biallelic:
                        self.counter['biallelic'] += read_no
                    else:
                        self.counter['mixed'] += read_no
                # Count reads with unwanted variants and continue
                if not snv and self.only_snv:
                    self.counter['unwanted_variants'] += read_no
                    continue
                if not biallelic and self.only_biallelic:
                    self.counter['unwanted_variants'] += read_no
                    continue
                # Count reads with excess variants and continue
                if flipvar.max_variant_count() > self.max_vars:
                    self.counter['excess_variants'] += read_no
                    continue
                # Count reads containing overlapping variants and continue
                if flipvar.variants_overlap():
                    self.counter['overlapping_variants'] += read_no
                    continue
                # Count paired reads with conflicting alleles and continue
                if flipvar.conflicting_alleles():
                    self.counter['conflicting_alleles'] += read_no
                    continue
                # Get flipped read haplotypes
                try:
                    flipvar.get_haplotypes(self.samples)
                except HaplotypeError:
                    self.counter['incomplete_haplotype'] += read_no
                    continue
                # Count flipped read haplotypes
                if flipvar.count_haplotypes() > self.max_seqs:
                    self.counter['excess_reads'] += read_no
                    continue
                # Generate flipped reads and check lengths
                flipvar.generate_flipped_reads()
                if flipvar.min_length() < self.min_len:
                    self.counter['too_short'] += read_no
                    continue
                # Generate fastq and write to file
                fastq_entry = flipvar.generate_fastq()
                self.fastq.write(fastq_entry)
                self.counter['to_remap'] += read_no
        # Get counts, adjust and write to file
        alignment_counts = self.bam_generator.counter
        alignment_counts['total'] += self.bam_generator.nocoordinate
        alignment_counts['unmapped'] += self.bam_generator.nocoordinate
        variant_counts = self.counter
        self.log.write_counts(
            alignment_counts=alignment_counts, variant_counts=variant_counts
        )


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
            "discarding the read (default=6)."
        )
    )
    parser.add_argument(
        "--max_seqs", type=int, default=64, help=(
            "The maximum number of allele flipped reads (or read pairs) "
            " to consider remapping. Reads with more allelic combinations "
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
        "--snv_only", action='store_true', help=(
            "Discard reads overlapping non snv variants (default=False)."
        )
    )
    parser.add_argument(
        "--biallelic_only", action='store_true', help=(
            "Discard reads overlapping non biallelic variants (default=False)."
        )
    )
    parser.add_argument(
        "--samples", nargs='+', help=(
            "Sample names for generation of phased allele flipped reads."
        )
    )
    parser.add_argument(
        "--assume_phased", action='store_true', help=(
            "Assume sample genotypes in VCF are phased (default=False)."
        )
    )
    args = parser.parse_args()
    # Run program
    check_phase = not args.assume_phased
    process_alignments = ProcessAlignments(
        bam=args.bam, vcf=args.vcf, out_prefix=args.out_prefix,
        min_mapq=args.min_mapq, max_vars=args.max_vars, max_seqs=args.max_seqs,
        min_len=args.min_len, only_snv=args.snv_only,
        only_biallelic=args.biallelic_only, samples=args.samples,
        check_phase=check_phase
    )
    process_alignments.filter_reads()
    process_alignments.close()
