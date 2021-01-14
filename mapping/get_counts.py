import argparse
import collections
import gzip
import pysam
import sys
import vartree


def bam_generator(
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
        # Process paired end reads
        if read.is_paired:
            # Return pair if other read has been stored
            if read.query_name in unpaired:
                if read.is_read1:
                    read1 = read
                    read2 = unpaired.pop(read.query_name)
                    assert(read2.is_read2)
                elif read.is_read2:
                    read2 = read
                    read1 = unpaired.pop(read.query_name)
                    assert(read1.is_read1)
                yield([read1, read2])
            # Store read if first observed of pair
            else:
                unpaired[read.query_name] = read
        # Return single end reads
        else:
            yield([read])
    # Print warning if unpaired reads remain
    if len(unpaired) > 0:
        sys.stderr.write(
            'WARNING: {} unpaired reads remain for chromsome {}\n'.format(
                len(unpaired), chromosome
            )
        )


def write_counts(
    outfile, counter
):
    # Extract variants and sort by positions
    variants = list(counter.keys())
    variants = sorted(variants, key=lambda x: x[1])
    # Get counts for each variant
    for variant in variants:
        variant_counts = counter[variant]
        # Extract data and readjust position
        chrom, position, ref, alt = variant
        position += 1
        # Process alternative alleles and their counts
        alt_list = alt.split(',')
        alt_list = [x if x != '' else '*' for x in alt_list]
        alt = ','.join(alt_list)
        indv_alt_counts = [variant_counts[i] for i in range(len(alt_list))]
        # Create output line
        out_line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            chrom, position, ref, alt, variant_counts['ref'],
            sum(indv_alt_counts), variant_counts['other'],
            ','.join(map(str, indv_alt_counts)),
        )
        outfile.write(out_line)


def count_vcf_variants(
    bam, vcf, outpath, partial
):
    # Open BAM file
    inbam = pysam.AlignmentFile(bam)
    # Open output file
    if outpath.endswith('.gz'):
        outfile = gzip.open(outpath, 'wt')
    else:
        outfile = open(outpath, 'wt')
    # Sequentially read in variants for each chromosome
    var_tree = vartree.VarTree(vcf)
    for chromosome in var_tree.chromosomes:
        var_tree.read_vcf(chromosome)
        # Create empty counter for chromosome
        counter_keys = [
            (chromosome, interval.begin, interval.data.ref, interval.data.alt)
            for interval in var_tree.tree
        ]
        counter = {key: collections.Counter() for key in counter_keys}
        # Loop through BAM reads for each chromosome
        for reads in bam_generator(inbam, chromosome):
            # Extract read sequences
            read_sequences = [read.query_sequence for read in reads]
            # Get read variants for single end reads
            if len(reads) == 1:
                read_variants = [
                    var_tree.get_read_variants(read=reads[0], partial=partial)
                ]
            # Get read variants for paired end reads
            if len(reads) == 2:
                read1_variants, read2_variants, identical = (
                    var_tree.get_paired_read_variants(
                        read1=reads[0], read2=reads[1], partial=partial
                    )
                )
                assert(identical)
                read_variants = [read1_variants, read2_variants]
            # Loop thorugh reads and alleles
            for sequence, variants in zip(read_sequences, read_variants):
                for position, variant in variants.items():
                    # Generate key for counter
                    counter_key = (
                        chromosome, position, variant.ref, variant.alt
                    )
                    # Get read allele and process alternativ alleles
                    read_allele = sequence[variant.start:variant.end]
                    alt_list = variant.alt.split(',')
                    # Count matches to refrence
                    if read_allele == variant.ref:
                        counter[counter_key]['ref'] += 1
                    # Count matched to alternative alleles individually
                    elif read_allele in alt_list:
                        for i in range(len(alt_list)):
                            if read_allele == alt_list[i]:
                                counter[counter_key][i] += 1
                    # Or count other
                    else:
                        counter[counter_key]['other'] += 1
        # Write counter to file
        write_counts(outfile=outfile, counter=counter)
    # Close output file
    outfile.close()


# Run script
if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Counts 'reference', 'alternative' and 'other' allele "
        "counts for variants overlapping reads in a BAM file. For paired end "
        "reads, variants overlapping both read1 and read2 are only counted "
        "once. The input BAM file is expected to only contain properly "
        "aligned reads. The output text file contains 8 columns: chromosome, "
        "position, reference allele, alternative allele(s), reference allele "
        "count, alternative allele(s) count, other allele(s) count, and comma "
        "seperated list on individual alternative allele counts."
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
    count_vcf_variants(
        bam=args.bam, vcf=args.vcf, outpath=args.outfile, partial=False
    )
