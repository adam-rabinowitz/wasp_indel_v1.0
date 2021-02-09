import argparse
import os
import pysam
import random
import vartree


class GenerateCounts(object):

    def __init__(
        self, bam, vcf, sample, partial
    ):
        # Store input arguments
        self.bam_path = bam
        self.vcf_path = vcf
        self.sample = sample
        self.partial = partial
        # Create VarTree object and open bam file
        self.vartree = vartree.VarTree(
            path=self.vcf_path, sample=self.sample
        )
        self.bam = pysam.AlignmentFile(self.bam_path)
        # Generate counter
        self.counter = {
            'non_biallelic': 0, 'no_haplotype': 0, 'homozygous': 0,
            'heterozygous': 0, 'no_variant': 0, 'biallelic_het_variant': 0,
            'other_variant': 0, 'ref': 0, 'alt': 0, 'other': 0
        }

    def close(
        self
    ):
        self.bam.close()

    def write_metrics(
        self, variant_metrics, outfile
    ):
        # Extract variants and sort by positions
        variants = list(variant_metrics.keys())
        variants = sorted(variants, key=lambda x: x.start)
        # Get counts for each variant
        for variant in variants:
            metrics = variant_metrics[variant]
            # Create output line
            out_line = (
                '{chromosome}\t{position}\t{id}\t{ref}\t{alts}\t{haplotype}\t'
                '{ref_prob}\t{het_prob}\t{alt_prob}\t{ref_count}\t'
                '{alt_count}\t{other_count}\n'
            ).format(
                chromosome=variant.chrom,
                position=variant.start + 1,
                id=variant.id,
                ref=variant.alleles[0],
                alts=','.join(variant.alleles[1:]),
                haplotype=metrics['haplotype'],
                ref_prob=metrics['ref_prob'],
                het_prob=metrics['het_prob'],
                alt_prob=metrics['alt_prob'],
                ref_count=metrics['ref_count'],
                alt_count=metrics['alt_count'],
                other_count=metrics['other_count']
            )
            outfile.write(out_line)

    def write_log(
        self, outfile
    ):
        # Create output
        log = (
            'Variants:\n'
            '  non biallelic: {non_biallelic}\n'
            '  no haplotype: {no_haplotype}\n'
            '  homozygous: {homozygous}\n'
            '  heterozygous: {heterozygous}\n'
            'Reads:\n'
            '  no variants: {no_variant}\n'
            '  biallelic heterozygotes variants: {biallelic_het_variant}\n'
            '  other variants: {other_variant}\n'
            'Read alelles:\n'
            '  ref: {ref}\n'
            '  alt: {alt}\n'
            '  other: {other}\n'
        ).format(**self.counter)
        outfile.write(log)

    def process_chromosome_variants(self, chromosome, outfile):
        # Read in variants for vcf
        self.vartree.read_vcf(chromosome)
        # Loop through variants and create default values
        variant_metrics = {}
        for interval in self.vartree.tree:
            variant = interval.data
            # Count and skip non-bial...
            if variant.is_biallelic():
                # Set default values for missing genotype and probabilities...
                if None in variant.gt or None in variant.pl:
                    self.counter['no_haplotype'] += 1
                    haplotype = 'NA'
                    ref_prob, het_prob, alt_prob = 'NA', 'NA', 'NA'
                # or extract and format genotype and probabilities
                else:
                    # Count variant type
                    if variant.is_heterozygous():
                        self.counter['heterozygous'] += 1
                    else:
                        self.counter['homozygous'] += 1
                    # Process haplotype and allele probabilites
                    haplotype = '|'.join(map(str, variant.gt))
                    probs = [10 ** ((-p) / 10) for p in variant.pl]
                    adj_probs = [p / sum(probs) for p in probs]
                    ref_prob, het_prob, alt_prob = [
                        '{:.2f}'.format(p) for p in adj_probs
                    ]
                # Add variant to dictionary
                variant_metrics[variant] = {
                    'haplotype': haplotype, 'ref_count': 0, 'alt_count': 0,
                    'other_count': 0, 'ref_prob': ref_prob,
                    'het_prob': het_prob, 'alt_prob': alt_prob
                }
            # or count and skip multiallelic variants
            else:
                self.counter['non_biallelic'] += 1
        # Get variant metrics
        if len(variant_metrics) > 0:
            for read in self.bam.fetch(contig=chromosome):
                # Extract read variants and their associayed metrics
                read_variants = self.vartree.get_read_variants(
                    read, partial=self.partial
                )
                # Get biallelic variants
                biallelic_het_variants = [
                    v for v in read_variants if v.is_biallelic_heterozygous()
                ]
                # Count and skip reads with no read variants
                if not read_variants:
                    self.counter['no_variant'] += 1
                # Count and skip reads with no biallelic heterozygous variants
                elif not biallelic_het_variants:
                    self.counter['other_variant'] += 1
                # Count and process reads with biallelic het variants
                else:
                    self.counter['biallelic_het_variant'] += 1
                    # Select variant and match read allele to variant alleles
                    selected_variant = random.choice(biallelic_het_variants)
                    try:
                        allele_index = selected_variant.alleles.index(
                            selected_variant.read_allele
                        )
                    except ValueError:
                        allele_index = None
                    # Count matches
                    if allele_index is None:
                        self.counter['other'] += 1
                        variant_metrics[selected_variant]['other_count'] += 1
                    elif allele_index == 0:
                        self.counter['ref'] += 1
                        variant_metrics[selected_variant]['ref_count'] += 1
                    elif allele_index == 1:
                        self.counter['alt'] += 1
                        variant_metrics[selected_variant]['alt_count'] += 1
            # Print metrics to file
            self.write_metrics(
                variant_metrics=variant_metrics, outfile=outfile
            )

    def process_all_variants(self, count_path, log_path):
        # Set counter to zero
        for key in self.counter:
            self.counter[key] = 0
        # Open count file
        with open(count_path, 'wt') as count_file:
            # Add commented input paths and parameters
            count_file.write(
                '#bam={}\n'.format(os.path.abspath(self.bam_path))
            )
            count_file.write(
                '#vcf={}\n'.format(os.path.abspath(self.vcf_path))
            )
            count_file.write('#sample={}\n'.format(self.sample))
            # Add commented header
            count_file.write(
                '#chrom\tposition\tid\tref\talt\thaplotype\tref_prob\t'
                'het_prob\talt_prob\tref_count\talt_count\tother_count\n'
            )
            # Get counts for each chromosome and add to file
            for chromosome in self.vartree.chromosomes:
                self.process_chromosome_variants(chromosome, count_file)
        # Write log to log file
        with open(log_path, 'wt') as log_file:
            self.write_log(log_file)


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
        "--bam", required=True, help=(
            "Filtered, coordinate sorted and indexed BAM file."
        )
    )
    parser.add_argument(
        "--vcf", required=True, help=(
            "Coordinate sorted, and tabix indexed, VCF file."
        )
    )
    parser.add_argument(
        "--sample", required=True, help=(
            "Sample for which to extract haplotype from VCF."
        )
    )
    parser.add_argument(
        "--outprefix", required=True, help=(
            "Prefix of output files."
        )
    )
    args = parser.parse_args()
    # Generate output files
    outfiles = {
        'initial': args.outprefix + '.variant_counts.txt',
        'compressed': args.outprefix + '.variant_counts.txt.gz',
        'log': args.outprefix + '.variant_counts_log.txt'
    }
    # Generate counts
    count_generator = GenerateCounts(
        bam=args.bam, vcf=args.vcf, sample=args.sample, partial=False
    )
    count_generator.process_all_variants(
        count_path=outfiles['initial'], log_path=outfiles['log']
    )
    count_generator.close()
    # Compress output file and remove uncompressed version
    pysam.tabix_compress(
        filename_in=outfiles['initial'], filename_out=outfiles['compressed'],
        force=True
    )
    os.remove(outfiles['initial'])
    # Index tabix file
    pysam.tabix_index(
        filename=outfiles['compressed'], seq_col=0, start_col=1, end_col=1,
        force=True, meta_char='#'
    )
