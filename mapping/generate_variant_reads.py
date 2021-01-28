import argparse
import collections
import gzip
import itertools
import pysam
import sys
import os
import util
import vartree


class DataFiles(object):
    """Object to hold names and filehandles for all input / output
        datafiles"""

    def __init__(
        self, bam, vcf, out_prefix, is_paired
    ):
        # Process supplied arguments
        self.in_bam_path = bam
        self.vcf = vcf
        self.prefix = out_prefix
        self.is_paired = is_paired
        # Print input path and open files
        sys.stderr.write(
            "reading reads from:\n  {}\n"
            "reading_variants from:\n  {}\n".format(
                self.in_bam_path, self.vcf
            )
        )
        self.in_bam = pysam.AlignmentFile(self.in_bam_path)
        # Create output directory
        output_dir = os.path.dirname(self.prefix)
        if output_dir and not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        # Create output paths
        self.invariant_bam_path = self.prefix + ".invariant.bam"
        self.remap_fastq_path = self.prefix + ".remap.fq.gz"
        self.log_path = self.prefix + '.variant_log.txt'
        # Print output paths and open files
        sys.stderr.write(
            "writing output files to:\n  {}\n  {}\n  {}\n\n".format(
                self.invariant_bam_path, self.remap_fastq_path,
                self.log_path
            )
        )
        self.invariant_bam = pysam.AlignmentFile(
            self.invariant_bam_path, 'wb', template=self.in_bam
        )
        self.remap_fastq = gzip.open(self.remap_fastq_path, 'wt')
        self.log = open(self.log_path, 'wt')

    def close(self):
        """close open filehandles"""
        for filehandle in (
            self.in_bam, self.invariant_bam, self.remap_fastq, self.log
        ):
            if filehandle:
                filehandle.close()


class ReadStats(object):
    """Track information about reads and SNPs that they overlap"""

    def __init__(self):
        # number of read matches to reference allele
        self.ref_count = 0
        # number of read matches to alternative allele
        self.alt_count = 0
        # number of reads that overlap SNP but match neither allele
        self.other_count = 0
        # total number of reads
        self.total_reads = 0
        # number of reads discarded because secondary match
        self.discard_secondary = 0
        # number of chimeric reads discarded
        self.discard_supplementary = 0
        # number of reads discarded becaused not mapped
        self.discard_unmapped = 0
        # number of reads discarded because mate unmapped
        self.discard_mate_unmapped = 0
        # paired reads map to different chromosomes
        self.discard_different_chromosome = 0
        # number of reads discarded due to low mapping quality
        self.discard_low_mapq = 0
        # number of reads discarded because not proper pair
        self.discard_improper_pair = 0
        # number of reads with abnormal alignments
        self.discard_abnormal_alignment = 0
        # read encompases overlapping alleles
        self.discard_overlapping_alleles = 0
        # number of reads discarded because of too many overlapping alleles
        self.discard_excess_alleles = 0
        # number of reads discarded because too many allelic combinations
        self.discard_excess_reads = 0
        # when read pairs share SNP locations but have different alleles there
        self.discard_discordant_shared_var = 0
        # number of single reads kept
        self.invariant_single = 0
        # number of read pairs kept
        self.invariant_pair = 0
        # number of single reads that need remapping
        self.remap_single = 0
        # number of read pairs kept
        self.remap_pair = 0

    def write(self, filehandle):
        # Write read processing counts
        counts = (
            "Total reads: {total}\n"
            "Discard reads:\n"
            "  secondary alignment: {secondary}\n"
            "  supplementary alignment: {supplementary}\n"
            "  unmapped: {unmapped}\n"
            "  mate unmapped: {mate_unampped}\n"
            "  different chromosomes: {different_chromosome}\n"
            "  low mapping quality: {low_mapq}\n"
            "  improper pair: {improper_pair}\n"
            "  abnormal alignment: {abnornmal_alignment}\n"
            "  overlapping alleles: {overlapping_alleles}\n"
            "  excess read alleles: {excess_alleles}\n"
            "  excess allelic combinations: {excess_reads}\n"
            "  discordant shared variants: {discordant_shared_var}\n"
            "Invariant reads:\n"
            "  single-end: {invariant_single}\n"
            "  paired: {invariant_pair}\n"
            "Remap reads:\n"
            "  single-end: {remap_single}\n"
            "  paired: {remap_pair}\n"
            "Allele counts:\n"
            "  reference count: {ref_count}\n"
            "  alternative count: {alt_count}\n"
            "  other count: {other_count}\n"
        ).format(
            total=self.total_reads,
            secondary=self.discard_secondary,
            supplementary=self.discard_supplementary,
            unmapped=self.discard_unmapped,
            mate_unampped=self.discard_mate_unmapped,
            different_chromosome=self.discard_different_chromosome,
            low_mapq=self.discard_low_mapq,
            improper_pair=self.discard_improper_pair,
            abnornmal_alignment=self.discard_abnormal_alignment,
            overlapping_alleles=self.discard_overlapping_alleles,
            excess_alleles=self.discard_excess_alleles,
            excess_reads=self.discard_excess_reads,
            discordant_shared_var=self.discard_discordant_shared_var,
            invariant_single=self.invariant_single,
            invariant_pair=self.invariant_pair,
            remap_single=self.remap_single,
            remap_pair=self.remap_pair,
            ref_count=self.ref_count,
            alt_count=self.alt_count,
            other_count=self.other_count
        )
        filehandle.write(counts)

    def check_vcf(self):
        # Calculate percentage of mismatched 'other' alleles
        total = self.ref_count + self.alt_count + self.other_count
        if total > 0:
            mismatch_pct = 100.0 * float(self.other_count) / total
            # Print warning message if mismatch rate is above 10%
            if mismatch_pct > 10.0:
                sys.stderr.write(
                    "WARNING: many read variants ({:.1f}%) do not match "
                    "either reference or alternative alleles. Coordinates "
                    "in input vcf may be incorrect.\n".format(mismatch_pct)
                )


# Function to count allele matches
def count_ref_alt_matches(
    read, variants, read_stats
):
    # Loop through varaints and extract read allele
    for variant in variants.values():
        read_allele = read.query_sequence[variant.start:variant.end]
        # Count if read allele matches, ref, alt or other
        if read_allele == variant.ref:
            read_stats.ref_count += 1
        elif read_allele in variant.alts:
            read_stats.alt_count += 1
        else:
            read_stats.other_count += 1


# Function generates alternative reads for supplied alleles
def generate_reads(
    sequence, quality, variants
):
    """Generate set of reads with all possible combinations"""
    # Create set to hold all current and novel versions of the read
    assert(len(sequence) == len(quality))
    Read = collections.namedtuple('Read', ['sequence', 'quality', 'edits'])
    initial_read = Read(sequence=sequence, quality=quality, edits={})
    current_reads = [initial_read]
    new_reads = []
    # Loop though variants end to start
    positions = list(variants.keys())
    positions.sort(reverse=True)
    for position in positions:
        variant = variants[position]
        # Calculate mean quality across reference
        if variant.start == variant.end:
            ref_quality = quality[(variant.start - 1):(variant.end + 1)]
        else:
            ref_quality = quality[variant.start:variant.end]
        mean_ref_quality = [sum(ref_quality) // len(ref_quality)]
        # Merge possible reference and alternative alleles
        allele_list = [variant.ref] + list(variant.alts)
        for new_allele in allele_list:
            for old_read in current_reads:
                # Skip identical alleles
                old_allele = old_read.sequence[variant.start:variant.end]
                if old_allele == new_allele:
                    continue
                # Skip '*' marking deletions spanning variant position
                if new_allele == '*':
                    continue
                # Create new sequence and quality
                new_sequence = (
                    old_read.sequence[:variant.start] +
                    new_allele +
                    old_read.sequence[variant.end:]
                )
                new_quality = (
                    old_read.quality[:variant.start] +
                    mean_ref_quality * len(new_allele) +
                    old_read.quality[variant.end:]
                )
                assert(len(new_sequence) == len(new_quality))
                # Create new edits
                new_edits = old_read.edits
                new_edits[position] = new_allele
                # Store modified read
                new_read = Read(
                    sequence=new_sequence, quality=new_quality, edits=new_edits
                )
                new_reads.append(new_read)
        # update current reads with new read versions
        current_reads.extend(new_reads)
        new_reads = []
    # Return all reads
    return(current_reads)


def read_pair_combos(read1_list, read2_list, max_seqs):
    """ Check common position of read pairs """
    # Set output variables
    pair_list = []
    # Loop through all combinations of reads
    for read1, read2 in itertools.product(read1_list, read2_list):
        # Check alleles are identical across common positions
        common = list(read1.edits.keys() & read2.edits.keys())
        if len(common) > 0:
            # Skip read pairs where common positions differ
            read1_common = [read1.edits[c] for c in common]
            read2_common = [read2.edits[c] for c in common]
            if read1_common != read2_common:
                continue
        # Add new alleles and check pait list length
        pair_list.append((read1, read2))
        if len(pair_list) > max_seqs:
            break
    # Return pair list and errors
    return(pair_list)


def write_fastq(
    fastq, read, new_reads, offset=33
):
    # Check reads
    assert(read.query_sequence == new_reads[0].sequence)
    # Get read name
    read_name = read.query_name
    reference_id = read.reference_id
    # Get position of original mapped reads and number of new reads
    position = "{id}-{start}-{end}".format(
        id=reference_id, start=read.reference_start,
        end=read.reference_end
    )
    n_seq = len(new_reads)
    # Process new reads sequentially
    for i, new_read in enumerate(new_reads):
        # Get new name
        fastq_name = "{}.{}.{}.{:06d}".format(
            read_name, position, n_seq, i
        )
        # Get new read sequence
        if read.is_reverse:
            fastq_seq = util.revcomp(new_read.sequence)
            fastq_qual = ''.join(
                [chr(x + offset) for x in new_read.quality]
            )
        else:
            fastq_seq = new_read.sequence
            fastq_qual = ''.join(
                [chr(x + offset) for x in new_read.quality[::-1]]
            )
        # Write read to file
        fastq_entry = '@{name}\n{seq}\n+\n{quality}\n'.format(
            name=fastq_name, seq=fastq_seq, quality=fastq_qual
        )
        fastq.write(fastq_entry)


def write_pair_fastq(
    fastq, read1, read2, new_pairs, offset=33
):
    # Check reads
    assert(read1.is_read1 & read2.is_read2)
    assert(read1.query_name == read2.query_name)
    assert(read1.reference_id == read2.reference_id)
    assert(read1.query_sequence == new_pairs[0][0].sequence)
    assert(read2.query_sequence == new_pairs[0][1].sequence)
    # Check name and reference id of original reads
    read_name = read1.query_name
    reference_id = read1.reference_id
    # Get position of original mapped reads and number of new reads
    position = "{id}-{start1}-{end1}-{start2}-{end2}".format(
        id=reference_id, start1=read1.reference_start,
        end1=read1.reference_end, start2=read2.reference_start,
        end2=read2.reference_end
    )
    n_pair = len(new_pairs)
    # Process new pairs sequentially
    for i, (new_read1, new_read2) in enumerate(new_pairs):
        # Generate read pair name
        fastq_name = "{}.{}.{}.{:06d}".format(read_name, position, n_pair, i)
        # Process each read in pair sequentially
        for old_read, new_read in [(read1, new_read1), (read2, new_read2)]:
            # Extract sequence and quality for read 1
            if old_read.is_reverse:
                fastq_seq = util.revcomp(new_read.sequence)
                fastq_qual = ''.join(
                    [chr(x + offset) for x in new_read.quality[::-1]]
                )
            else:
                fastq_seq = new_read.sequence
                fastq_qual = ''.join(
                    [chr(x + offset) for x in new_read.quality]
                )
            assert(len(fastq_seq) == len(fastq_qual))
            # Write fastq
            fastq_entry = '@{name}\n{seq}\n+\n{quality}\n'.format(
                name=fastq_name, seq=fastq_seq, quality=fastq_qual
            )
            fastq.write(fastq_entry)


def variants_overlap(
    variants
):
    """Function checks if variants in dictionary overlap"""
    overlap = False
    start_position = 0
    for variant in variants.values():
        if variant.start < start_position:
            overlap = True
            break
        start_position = variant.end
    return(overlap)


def process_single_read(
    read, read_stats, files, var_tree, max_seqs, max_vars
):
    """Check if a single read overlaps SNPs or indels, and writes
    this read (or generated read pairs) to appropriate output files"""
    # check if read overlaps variants
    try:
        read_variants = var_tree.get_read_variants(read=read, partial=True)
    except AssertionError:
        read_variants = None
    # Count and skip abnormal reads
    if read_variants is None:
        read_stats.discard_abnormal_alignment += 1
    # Keep reads without variants
    elif len(read_variants) == 0:
        files.invariant_bam.write(read)
        read_stats.invariant_single += 1
    # Count and skip overlapping variants
    if variants_overlap(read_variants):
        read_stats.discard_overlapping_alleles += 1
    # count and skip excess variants
    elif len(read_variants) > max_vars:
        read_stats.discard_excess_alleles += 1
    # Process passed variants
    else:
        # Count reference and alternative allele matches
        count_ref_alt_matches(
            read=read, variants=read_variants, read_stats=read_stats
        )
        # Generate new reads
        new_reads = generate_reads(
            sequence=read.query_sequence,
            quality=list(read.query_qualities),
            variants=read_variants
        )
        # Skip reads with too many variants
        if len(new_reads) > max_seqs:
            read_stats.discard_excess_reads += 1
        # Write passed pairs to file
        else:
            write_fastq(
                fastq=files.remap_fastq, read=read, new_reads=new_reads
            )
            read_stats.remap_single += 1


def process_paired_read(
    read1, read2, read_stats, files, var_tree, max_vars, max_seqs
):
    """Checks if either end of read pair overlaps SNPs or indels
    and writes read pair (or generated read pairs) to appropriate
    output files"""
    # Find variants for each read and total variant count
    try:
        read1_variants, read2_variants, identical_variants = (
            var_tree.get_paired_read_variants(
                read1=read1, read2=read2, partial=True
            )
        )
        variant_count = len(read1_variants.keys() | read2_variants.keys())
    except AssertionError:
        variant_count = None
    # Count and skip abnormal reads
    if variant_count is None:
        read_stats.discard_abnormal_alignment += 2
    # Keep reads without variants
    elif variant_count == 0:
        files.invariant_bam.write(read1)
        files.invariant_bam.write(read2)
        read_stats.invariant_pair += 2
    # Count and skip mismatched variants
    elif not identical_variants:
        read_stats.discard_discordant_shared_var += 2
    # Count and skip overlapping variants
    elif variants_overlap(read1_variants) or variants_overlap(read2_variants):
        read_stats.discard_overlapping_alleles += 2
    # count and discard excess variants or...
    elif len(read1_variants) > max_vars or len(read2_variants) > max_vars:
        read_stats.discard_excess_alleles += 2
    # count and discard excess reads...
    elif (2 ** variant_count) > max_seqs:
        read_stats.discard_excess_reads += 2
    # or process variant reads
    else:
        # Count reference and alternative allele matches
        count_ref_alt_matches(
            read=read1, variants=read1_variants, read_stats=read_stats
        )
        count_ref_alt_matches(
            read=read2, variants=read2_variants, read_stats=read_stats
        )
        # Generate new reads
        new_read1 = generate_reads(
            sequence=read1.query_sequence, quality=list(read1.query_qualities),
            variants=read1_variants
        )
        new_read2 = generate_reads(
            sequence=read2.query_sequence, quality=list(read2.query_qualities),
            variants=read2_variants
        )
        # Get all paired combinations of reads
        new_pairs = read_pair_combos(
            read1_list=new_read1, read2_list=new_read2, max_seqs=max_seqs
        )
        # Discard excess pairs
        if len(new_pairs) > max_seqs:
            read_stats.discard_excess_reads += 2
        # Write passed pairs to file
        else:
            write_pair_fastq(
                fastq=files.remap_fastq, read1=read1, read2=read2,
                new_pairs=new_pairs
            )
            read_stats.remap_pair += 2


def filter_reads(
    files, max_seqs, max_vars, min_mapq
):
    '''Main function to filter reads within sorted BAM file'''
    # Set variables to process read chromosome
    cur_chrom = None
    seen_chrom = set([])
    # Set variables for processing reads
    var_tree = vartree.VarTree(files.vcf)
    read_stats = ReadStats()
    read_pair_cache = {}
    # Loop through and count read in bam file
    for read in files.in_bam:
        read_stats.total_reads += 1
        # Skip and count secondary reads
        if read.is_secondary:
            read_stats.discard_secondary += 1
            continue
        # Skip and count supplementary reads
        if read.is_supplementary:
            read_stats.discard_supplementary += 1
            continue
        # Skip and count unmapped reads
        if read.is_unmapped:
            read_stats.discard_unmapped += 1
            continue
        # Skip and count unmapped pairs
        if read.is_paired and read.mate_is_unmapped:
            read_stats.discard_mate_unmapped += 1
            continue
        # Skip reads mapping to different chromosomes
        if read.is_paired and read.reference_id != read.next_reference_id:
            read_stats.discard_different_chromosome += 1
            continue
        # Process reads on a new chromosome
        if read.reference_name != cur_chrom:
            # Check all pairs have been found
            if len(read_pair_cache) != 0:
                raise ValueError(
                    'unpaired reads for chromosome: {}. '
                    'Has alignment been filtered?'.format(cur_chrom)
                )
            # Check that input bam file is sorted
            cur_chrom = read.reference_name
            if cur_chrom in seen_chrom:
                raise ValueError("BAM file is not sorted")
            seen_chrom.add(cur_chrom)
            # Read variants from file
            sys.stderr.write(
                "reading variants for chromosome {}\n".format(cur_chrom)
            )
            var_tree.read_vcf(cur_chrom)
            sys.stderr.write("processing reads\n")
        # Find single end or complete pairs
        if read.is_paired:
            # Store unpaired reads and continue to next read
            if read.query_name not in read_pair_cache:
                found = None
                read_pair_cache[read.query_name] = read
                continue
            # Return pair if other read has been stored
            else:
                if read.is_read1:
                    read1 = read
                    read2 = read_pair_cache.pop(read.query_name)
                    assert(read2.is_read2)
                elif read.is_read2:
                    read2 = read
                    read1 = read_pair_cache.pop(read.query_name)
                    assert(read1.is_read1)
                found = [read1, read2]
        else:
            found = [read]
        # Count and skip poorly mapped reads
        if any([read.mapping_quality < min_mapq for read in found]):
            read_stats.discard_low_mapq += len(found)
            continue
        # Process paired end reads
        if len(found) == 2:
            # Count and skip improper pairs
            if not all([read.is_proper_pair for read in found]):
                read_stats.discard_improper_pair += 2
                continue
            # Process paired end reads
            process_paired_read(
                read1=found[0], read2=found[1], read_stats=read_stats,
                files=files, var_tree=var_tree, max_seqs=max_seqs,
                max_vars=max_vars
            )
        # Process single end reads
        else:
            process_single_read(
                read=found[0], read_stats=read_stats, files=files,
                var_tree=var_tree, max_seqs=max_seqs, max_vars=max_vars
            )
    # Check all pairs have been found
    if len(read_pair_cache) != 0:
        raise ValueError(
            'unpaired reads for chromosome: {}. '
            'Has alignment been filtered?'.format(cur_chrom)
        )
    # Write read stats to file and check vcf
    read_stats.write(files.log)
    read_stats.check_vcf()


# Run script
if __name__ == '__main__':
    # Print command and software versions
    sys.stderr.write("command line:\n  %s\n" % " ".join(sys.argv))
    sys.stderr.write("python version:\n  %s\n" % sys.version)
    sys.stderr.write("pysam version:\n  %s\n\n" % pysam.__version__)
    # Check versions of programs
    util.check_python_version()
    util.check_pysam_version()
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
        "--paired_end", action='store_true', default=False,
        help=(
            "Reads are paired-end (default is single)."
        )
    )
    parser.add_argument(
        "--min_mapq", type=int, default=0, help=(
            "Minimum read mapping quality (default=0)."
        )
    )
    parser.add_argument(
        "--max_seqs", type=int, default=64, help=(
            "The maximum number of sequences with different allelic "
            "combinations to consider remapping (default=64). Read pairs "
            "with more allelic combinations than MAX_SEQs are discarded"
        )
    )
    parser.add_argument(
        "--max_vars", type=int, default=6, help=(
            "The maximum number of variants allowed to overlap a read before "
            "discarding the read. Allowing higher numbers will decrease speed "
            "and increase memory usage (default=6)."
        )
    )
    parser.add_argument(
        "bam", action='store', help=(
            "Coordinate-sorted input BAM file."
        )
    )
    parser.add_argument(
        "vcf", action='store', help=(
            "Coordinate sorted, and tabix indexed, VCF file."
        )
    )
    parser.add_argument(
        "out_prefix", action='store', help=(
            "Prefix of output files."
        )
    )
    args = parser.parse_args()
    # Run program
    files = DataFiles(
        bam=args.bam, is_paired=args.paired_end,
        out_prefix=args.out_prefix, vcf=args.vcf
    )
    filter_reads(
        files=files, max_seqs=args.max_seqs, max_vars=args.max_vars,
        min_mapq=args.min_mapq
    )
    files.close()
