import argparse
import collections
import os
import pysam
import random
from mapping.vartree import VarTree


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
        self.vartree = VarTree(
            path=self.vcf_path, samples=[self.sample], check_phase=False
        )
        self.bam = pysam.AlignmentFile(self.bam_path)
        # Generate counter
        self.counter = {
            'non_biallelic': 0, 'no_haplotype': 0, 'homozygous': 0,
            'heterozygous': 0, 'no_variant': 0, 'no_bi_variant': 0,
            'bi_het_variant': 0, 'bi_nonhet_variant': 0,
            'ref': 0, 'alt': 0, 'other': 0
        }

    def close(
        self
    ):
        self.bam.close()

    def write_metrics(
        self, variant_metrics, outfile
    ):
        # Get counts for each variant
        for variant in variant_metrics.values():
            # Create output line
            out_line = (
                '{chromosome}\t{position}\t{id}\t{ref}\t{alts}\t{haplotype}\t'
                '{ref_prob}\t{het_prob}\t{alt_prob}\t{ref_as_count}\t'
                '{alt_as_count}\t{other_as_count}\t{ref_total_count}\t'
                '{alt_total_count}\t{other_total_count}\n'
            ).format(**variant)
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
            '  no biallelic variants: {no_bi_variant}\n'
            '  biallelic non-heterozygous variants: {bi_nonhet_variant}\n'
            '  biallelic heterozygotes variants: {bi_het_variant}\n'
            'Alelle specific reads:\n'
            '  ref: {ref}\n'
            '  alt: {alt}\n'
            '  other: {other}\n'
        ).format(**self.counter)
        outfile.write(log)

    def process_chromosome_variants(self, chromosome, outfile):
        # Read in variants for vcf
        self.vartree.read_vcf(chromosome)
        # Set random seed so that ouptut is reproducible
        random.seed(42)
        # Loop through variants and create default values
        variant_metrics = collections.OrderedDict()
        for variant in self.vartree.variants:
            # Count and skip non-bial...
            if variant.is_biallelic():
                # Set default values for missing genotype and probabilities...
                if None in variant.genotypes[self.sample]:
                    self.counter['no_haplotype'] += 1
                    haplotype = 'NA'
                    ref_prob, het_prob, alt_prob = 'NA', 'NA', 'NA'
                # or extract and format genotype and probabilities
                else:
                    # Count variant type
                    if variant.is_heterozygous(self.sample):
                        self.counter['heterozygous'] += 1
                    else:
                        self.counter['homozygous'] += 1
                    # Process haplotype and allele probabilites
                    haplotype = '|'.join(
                        map(str, variant.genotypes[self.sample])
                    )
                    ref_prob, het_prob, alt_prob = [
                        '{:.2f}'.format(p) for p in variant.probs[self.sample]
                    ]
                # Add variant to dictionary
                assert(variant.id not in variant_metrics)
                variant_metrics[variant.id] = {
                    'chromosome': chromosome, 'position': variant.start + 1,
                    'id': variant.id, 'ref': variant.alleles[0],
                    'alts': ','.join(variant.alleles[1:]),
                    'haplotype': haplotype, 'ref_as_count': 0,
                    'alt_as_count': 0, 'other_as_count': 0,
                    'ref_total_count': 0, 'alt_total_count': 0,
                    'other_total_count': 0, 'ref_prob': ref_prob,
                    'het_prob': het_prob, 'alt_prob': alt_prob
                }
            # or count and skip multiallelic variants
            else:
                self.counter['non_biallelic'] += 1
        # Get variant metrics
        if len(variant_metrics) > 0:
            for read in self.bam.fetch(contig=chromosome):
                # Count and skip reads without variants
                read_variants = self.vartree.get_read_variants(
                    read, partial=self.partial
                )
                if not read_variants:
                    self.counter['no_variant'] += 1
                    continue
                # Count and skip reads without biallelic variants
                bi_variants = [
                    rv for rv in read_variants if rv.is_biallelic()
                ]
                if not bi_variants:
                    self.counter['no_bi_variant'] += 1
                    continue
                # Add all allele counts for all biallelic variants
                for bi in bi_variants:
                    # Get read allele
                    try:
                        allele_index = bi.alleles.index(bi.read_allele)
                    except ValueError:
                        allele_index = None
                    # Add allele to counts
                    if allele_index is None:
                        variant_metrics[bi.id]['other_total_count'] += 1
                    elif allele_index == 0:
                        variant_metrics[bi.id]['ref_total_count'] += 1
                    elif allele_index == 1:
                        variant_metrics[bi.id]['alt_total_count'] += 1
                # Count and skip reads without biallelic heterozygous variants
                bihet_variants = [
                    bv for bv in bi_variants if bv.is_heterozygous(self.sample)
                ]
                if not bihet_variants:
                    self.counter['bi_nonhet_variant'] += 1
                    continue
                # Count biallelic het read and select biallelic het variant
                self.counter['bi_het_variant'] += 1
                bihet = random.choice(bihet_variants)
                # Select variant and match read allele to variant alleles
                try:
                    allele_index = bihet.alleles.index(bihet.read_allele)
                except ValueError:
                    allele_index = None
                # Count matches for selec
                if allele_index is None:
                    self.counter['other'] += 1
                    variant_metrics[bihet.id]['other_as_count'] += 1
                elif allele_index == 0:
                    self.counter['ref'] += 1
                    variant_metrics[bihet.id]['ref_as_count'] += 1
                elif allele_index == 1:
                    self.counter['alt'] += 1
                    variant_metrics[bihet.id]['alt_as_count'] += 1
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
                'het_prob\talt_prob\tref_as_count\talt_as_count\t'
                'other_as_count\tref_total_count\talt_total_count\t'
                'other_total_count\n'
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
        description="Extracts key metrics from the BAM file including "
        "allele counts for variants contained within each alignment. The "
        "input BAM file is expected to only contain properly aligned reads. "
        "The output text file contains 15 columns. The first 6 columns are "
        "chromosome, position, variant ID, reference allele, alternative "
        "allele and haplotype. Columns 7-9 are the probability that the "
        "genotype is homozygous reference allele, heterzoygous and homozygous "
        "alternative allele, respectively. Columns 10-12 contain the allele "
        "specific counts for the reference, alternative and other alleles. "
        "The allele specific counts are generated by randomly selecting a "
        "biallelic heterozygous variant, contained within each alignment, and "
        "determing the variant allele. Columns 13-15 contain the total "
        "reference, alternative and other alleles counts. These counts are "
        "generated by determining the variant allele for all variants "
        "contained in each alignment."
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
