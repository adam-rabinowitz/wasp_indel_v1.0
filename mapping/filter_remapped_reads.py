import argparse
from pysam import AlignmentFile  # pylint: disable=no-name-in-module


class ReadStats:

    """Track information about what the program is doing with reads"""
    def __init__(self):
        # total number of reads
        self.total_reads = 0
        # number of reads with missing alignments
        self.discard_missing_reads = 0
        # number of reads with unmapped alignments
        self.discard_unmapped_reads = 0
        # number of reads with low mapping quality
        self.discard_low_mapq = 0
        # number of reads with mismapped original
        self.discard_original_mismapped = 0
        # number of reads with mismapped variants
        self.discard_variant_mismapped = 0
        # number of reads with identical mapping
        self.keep_reads = 0

    def write(self, filehandle):
        counts = (
            "Read counts:\n"
            "  total: {total}\n"
            "  missing: {missing}\n"
            "  unmapped: {unmapped}\n"
            "  low mapping quality: {low_mapq}\n"
            "  mismapped original: {original_mismapped}\n"
            "  mismapped variant: {variant_mismapped}\n"
            "  kept reads: {keep}\n"
        ).format(
            total=self.total_reads,
            missing=self.discard_missing_reads,
            unmapped=self.discard_unmapped_reads,
            low_mapq=self.discard_low_mapq,
            original_mismapped=self.discard_original_mismapped,
            variant_mismapped=self.discard_variant_mismapped,
            keep=self.keep_reads
        )
        filehandle.write(counts)


def bam_generator(
    bam
):
    # Set processing variables
    seen_names = set()
    current_name = None
    reads = []
    # Loop through yields in bam file
    for read in bam:
        # Skip secondary and supplementary read
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        # Process query name
        words = read.query_name.split(".")
        if len(words) < 4:
            raise ValueError(
                "expected read names to be formatted like <orig_name>."
                "<coordinate>.<read_number>.<total_read_number> but got "
                "{}".format(read.query_name)
            )
        # Extract original name accouting for that it might contain "."
        original_name = ".".join(words[0:len(words)-3])
        # Check file has been name sorted
        if original_name in seen_names:
            raise IOError('check input file is name sorted')
        # Process stored reads if read has a new name
        if original_name != current_name:
            # Return current reads
            if current_name is not None:
                seen_names.add(current_name)
                yield(current_name, reads)
            # Set current name and reset reads
            current_name = original_name
            reads = []
        # Add reads
        reads.append(read)
    # Return final set od reads
    if current_name is not None:
        yield(current_name, reads)


def differ(positions, max_diff):
    '''determined if a list of integers differ by more than the permitted
    value from the first integer in the list'''
    expected = positions[0]
    do_differ = [abs(x - expected) > max_diff for x in positions]
    any_differ = any(do_differ)
    return(any_differ)


def filter_bam_single(
    inpath, outpath, min_mapq, max_diff
):
    # Create object to collect read stats
    read_stats = ReadStats()
    # Open bam files
    inbam = AlignmentFile(inpath)
    outbam = AlignmentFile(outpath, 'wb', template=inbam)
    # Loop through samples in bam file
    for name, reads in bam_generator(inbam):
        read_stats.total_reads += 1
        # Get read location and expected no
        old_location, expected_reads = reads[0].query_name.split('.')[-3:-1]
        old_location = [int(x) for x in old_location.split('-')]
        expected_reads = int(expected_reads)
        # Count and discrad reads with missing alignments
        if len(reads) != expected_reads:
            assert(len(reads) < expected_reads)
            read_stats.discard_missing_reads += 1
            continue
        # Count and skip reads with unmapped variants
        if any([read.is_unmapped for read in reads]):
            read_stats.discard_unmapped_reads += 1
            continue
        # Count and skip reads with poorly mapped variants
        if any([read.mapping_quality < min_mapq for read in reads]):
            read_stats.discard_low_mapq += 1
            continue
        # Check read numbers
        expected_read_nos = list(range(expected_reads))
        read_nos = [int(read.query_name.split('.')[-1]) for read in reads]
        assert(read_nos == expected_read_nos)
        # Extract mapped intervals
        read_intervals = [
            (read.reference_id, read.reference_start, read.reference_end)
            for read in reads
        ]
        # Check original read maps to the same location
        new_location = list(read_intervals[0])
        if new_location != old_location:
            read_stats.discard_original_mismapped += 1
            continue
        # Extract location of all mapped reads
        read_id_list, read_start_list, read_end_list = zip(*read_intervals)
        if len(set(read_id_list)) > 1:
            read_stats.discard_variant_mismapped += 1
            continue
        if differ(read_start_list, max_diff):
            read_stats.discard_variant_mismapped += 1
            continue
        if differ(read_end_list, max_diff):
            read_stats.discard_variant_mismapped += 1
            continue
        # Keep remaining reads
        out_read = reads[0]
        out_read.query_name = name
        outbam.write(out_read)
        read_stats.keep_reads += 1
    # Return read stats
    return(read_stats)


def filter_bam_pairs(
    inpath, outpath, min_mapq, max_diff
):
    # Create object to collect read stats
    read_stats = ReadStats()
    # Open bam files
    inbam = AlignmentFile(inpath)
    outbam = AlignmentFile(outpath, 'wb', template=inbam)
    # Loop through samples in bam file
    for name, reads in bam_generator(inbam):
        read_stats.total_reads += 2
        # Get read location and expected no
        old_location, expected_no = reads[0].query_name.split('.')[-3:-1]
        old_location = [int(x) for x in old_location.split('-')]
        expected_pairs = int(expected_no)
        expected_reads = 2 * expected_pairs
        # Check expected read number
        if len(reads) != expected_reads:
            assert(len(reads) < expected_reads)
            read_stats.discard_missing_reads += 2
            continue
        # Count and skip reads with unmapped variants
        if any([read.is_unmapped for read in reads]):
            read_stats.discard_unmapped_reads += 2
            continue
        # Count and skip reads with poorly mapped variants
        if any([read.mapping_quality < min_mapq for read in reads]):
            read_stats.discard_low_mapq += 2
            continue
        # Split reads into read1 and read2
        read1 = [read for read in reads if read.is_read1]
        read2 = [read for read in reads if read.is_read2]
        # Check read numbers
        expected_read_nos = list(range(expected_pairs))
        read1_nos = [int(read.query_name.split('.')[-1]) for read in read1]
        read2_nos = [int(read.query_name.split('.')[-1]) for read in read2]
        assert(read1_nos == expected_read_nos)
        assert(read2_nos == expected_read_nos)
        # Extract mapped intervals
        read1_intervals = [
            [read.reference_id, read.reference_start, read.reference_end]
            for read in read1
        ]
        read2_intervals = [
            [read.reference_id, read.reference_start, read.reference_end]
            for read in read2
        ]
        # Check original read maps to the same location
        orig_read1_interval = read1_intervals[0]
        orig_read2_interval = read2_intervals[0]
        if orig_read1_interval[0] != orig_read2_interval[0]:
            read_stats.discard_original_mismapped += 2
            continue
        new_location = orig_read1_interval + orig_read2_interval[1:3]
        if new_location != old_location:
            read_stats.discard_original_mismapped += 2
            continue
        # Transpose list of intervals
        read1_id_list, read1_start_list, read1_end_list = zip(
            *read1_intervals
        )
        read2_id_list, read2_start_list, read2_end_list = zip(
            *read2_intervals
        )
        # Count and skip mismatched variant reads
        if len(set(read1_id_list)) > 1:
            read_stats.discard_variant_mismapped += 2
            continue
        if differ(read1_start_list, max_diff):
            read_stats.discard_variant_mismapped += 2
            continue
        if differ(read1_end_list, max_diff):
            read_stats.discard_variant_mismapped += 2
            continue
        if len(set(read2_id_list)) > 1:
            read_stats.discard_variant_mismapped += 2
            continue
        if differ(read2_start_list, max_diff):
            read_stats.discard_variant_mismapped += 2
            continue
        if differ(read2_end_list, max_diff):
            read_stats.discard_variant_mismapped += 2
            continue
        # Keep remaining reads
        for out_read in (read1[0], read2[0]):
            out_read.query_name = name
            out_read.is_proper_pair = True
            outbam.write(out_read)
        read_stats.keep_reads += 2
    # Return read stats
    return(read_stats)


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="This program checks whether reads that overlap "
        "variants map back to the same location as the original reads "
        "after their alleles are flipped by the find_intersecting_snps.py "
        "script. Reads where one or more allelic versions map to a "
        "different location (or fail to map) are discarded. Reads in "
        "the input bam file are expected to have read names encoding "
        "the original map location and number of allelic variants. "
        "Specifically, the read names should be delimited with the '.' "
        "character and contain the following fields: <original_name>."
        "<coordinate>.<total_read_number>.<read_number>. These read names "
        "are generated by the find_intersecting_snps.py script. Reads with "
        "the same original name should be adjacent within the BAM file. "
        "Two output files are generated with the provided prefix and the "
        "following suffixes. 1) '.consistent.bam' - BAM file containing reads "
        "whose allele flipped variants map to the same location. 2) "
        "'.filter_log.txt' -  text file containing read filtering metrics."
    )
    parser.add_argument(
        "--min_mapq", type=int, default=0, help=(
            "Minimum read mapping quality for remapped reads (default=0)."
        )
    )
    parser.add_argument(
        "--wobble", type=int, default=0, help=(
            "Maximum allowed difference between original mapped position "
            "and allele flipped positions (default=0)."
        )
    )
    parser.add_argument(
        "bam", help=(
            "Input BAM file containing alignments of allele flipped reads."
        )
    )
    parser.add_argument(
        "out_prefix", help=(
            "output BAM file containing alignments where allele flipped "
            "reads map to the same location."
        )
    )
    args = parser.parse_args()
    # Get BAM paiting
    with AlignmentFile(args.bam) as bam:
        paired = [read.is_paired for read in bam.head(100)]
        if any(paired):
            assert(all(paired))
            paired = True
        else:
            paired = False
    # Filter reads
    out_bam = args.out_prefix + '.consistent.bam'
    if paired:
        read_stats = filter_bam_pairs(
            inpath=args.bam, outpath=out_bam, min_mapq=args.min_mapq,
            max_diff=args.wobble
        )
    else:
        read_stats = filter_bam_single(
            inpath=args.bam, outpath=out_bam, min_mapq=args.min_mapq,
            max_diff=args.wobble
        )
    # Write readstats to file
    with open(args.out_prefix + '.second_alignment_log.txt', 'wt') as log:
        read_stats.write(log)
