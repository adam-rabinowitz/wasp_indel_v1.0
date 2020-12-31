import argparse
import collections
import gzip
import pysam
import sys
import vartable


def bam_generator(
    bam, chromosome
):
    # Loop through reads on chromosome
    for read in bam.fetch(chromosome):
        # Skip unmapped, secondary and supplementary alignments
        if read.is_unmapped:
            continue
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        # Yield read
        yield(read)


def paired_bam_generator(
    bam, chromosome
):
    # Set processing variables
    unpaired = {}
    # Loop through reads on chromosome
    for read in bam.fetch(chromosome):
        # Skip unmapped, secondary and supplementary alignments
        if read.is_unmapped:
            continue
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        # Extract read1 and read2 for previously observed reads...
        query_name = read.query_name
        if query_name in unpaired:
            if read.is_read1:
                read1 = read
                read2 = unpaired.pop(query_name)
                assert(read2.is_read2)
            elif read.is_read2:
                read2 = read
                read1 = unpaired.pop(query_name)
                assert(read1.is_read1)
            yield(read1, read2)
        # Or add read to dictionary for new reads
        else:
            unpaired[query_name] = read
    # Print warning if unpaired reads remain
    if len(unpaired) > 0:
        sys.stderr.write(
            'WARNING: {} unpaired reads remain for chromsome {}\n'.format(
                len(unpaired), chromosome
            )
        )


def count_chromosome_variants(
    bam, vcf, chromosome, partial
):
    # Create empty counter
    counter = collections.defaultdict(collections.Counter)
    # Read in variants
    var_table = vartable.VarTable(vcf)
    var_table.read_vcf(chromosome)
    # Loop through bam file
    for read in bam_generator(bam, chromosome):
        # Get read sequence and variants
        read_sequence = read.query_sequence
        read_variants = var_table.get_read_variants(read=read, partial=False)
        # Count alleles
        for position, variant in read_variants.items():
            # Generate counter for position
            variant_tuple = (chromosome, position, variant.ref, variant.alt)
            # Find if read is variant allele
            read_allele = read_sequence[variant.start:variant.end]
            if read_allele == variant.ref:
                counter[variant_tuple]['ref'] += 1
            elif read_allele in variant.alt.split(','):
                counter[variant_tuple]['alt'] += 1
            else:
                counter[variant_tuple]['other'] += 1
    # Return counter
    return(counter)


def count_paired_chromosome_variants(
    bam, vcf, chromosome, partial
):
    # Create empty counter
    counter = collections.defaultdict(collections.Counter)
    # Read in variants
    var_table = vartable.VarTable(vcf)
    var_table.read_vcf(chromosome)
    # Loop through bam file
    for read1, read2 in paired_bam_generator(bam, chromosome):
        # Get sequence variants and check they agree
        read1_variants, read2_variants, identical = (
            var_table.get_paired_read_variants(
                read1=read1, read2=read2, partial=False
            )
        )
        if not identical:
            continue
        # Remove common variants from read2 so as not to count twice
        common_positions = read1_variants.keys() & read2_variants.keys()
        for position in common_positions:
            read2_variants.pop(position)
        # Get read sequence
        read1_sequence = read1.query_sequence
        read2_sequence = read2.query_sequence
        # Loop through reads in pairs
        for read_sequence, read_variants in [
            (read1_sequence, read1_variants), (read2_sequence, read2_variants)
        ]:
            # Loop through read variants
            for position, variant in read_variants.items():
                # Generate counter for position
                variant_tuple = (
                    chromosome, position, variant.ref, variant.alt
                )
                # Find if read is variant allele
                read_allele = read_sequence[variant.start:variant.end]
                if read_allele == variant.ref:
                    counter[variant_tuple]['ref'] += 1
                elif read_allele in variant.alt.split(','):
                    counter[variant_tuple]['alt'] += 1
                else:
                    counter[variant_tuple]['other'] += 1
    # Return counter
    return(counter)


def write_counts(
    outfile, counter
):
    # Extract variants and sort by positions
    variants = list(counter.keys())
    variants = sorted(variants, key=lambda x: x[1])
    # Set file parameters
    for variant in variants:
        variant_counts = counter[variant]
        out_line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            *variant, variant_counts['ref'], variant_counts['alt'],
            variant_counts['other']
        )
        outfile.write(out_line)


def bam_variant_counts(
    bam, vcf, outpath, paired, partial
):
    # Open input files
    inbam = pysam.AlignmentFile(bam)
    if outpath.endswith('.gz'):
        outfile = gzip.open(outpath, 'wt')
    else:
        outfile = open(outpath, 'wt')
    # Loop through chromosomes
    for chromosome in inbam.references:
        # Get counts for each chromosome
        if paired:
            counter = count_paired_chromosome_variants(
                bam=inbam, vcf=vcf, chromosome=chromosome, partial=partial
            )
        else:
            counter = count_chromosome_variants(
                bam=inbam, vcf=vcf, chromosome=chromosome, partial=partial
            )
        # Write count to file
        write_counts(outfile, counter)


# Run script
if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Counts 'reference', 'alternative' and 'other' allele "
        "counts for variants overlapping reads in a BAM file. For paired end "
        "reads, variants overlapping both read1 and read2 are only counted "
        "once. The input BAM file is expected to only contain properly "
        "aligned reads. The output text file contains 7 columns: chromosome, "
        "position, reference allele, alternative allele, reference allele "
        "count, alternative allele count and other allele count."
    )
    parser.add_argument(
        "--paired_end", action='store_true', default=False,
        help=(
            "Indicates that reads are paired-end (default is single)."
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
        "outfile", action='store', help=(
            "Path of output text files."
        )
    )
    args = parser.parse_args()
    # Generate counts
    bam_variant_counts(
        bam=args.bam, vcf=args.vcf, outpath=args.outfile,
        paired=args.paired_end, partial=False
    )
