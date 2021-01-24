import argparse
import collections
import intervaltree
import pysam
import random
import vartree


class bam_coverage(object):

    def __init__(self, seed=42):
        # Initialise trees
        self.ref = intervaltree.IntervalTree()
        self.alt = intervaltree.IntervalTree()
        self.other = intervaltree.IntervalTree()
        # Set random seed
        random.seed(seed)
        # Set hets
        self.hets = ('0/1', '1/0')
        # Create count tuple
        self.counts = collections.namedtuple(
            'counts', ['ref', 'alt', 'other']
        )

    def add_reads(self, reads, genotypes, alleles):
        # Get allele
        assert(len(genotypes) == len(alleles))
        if len(genotypes) > 0:
            het_alleles = [
                a for g, a in zip(genotypes, alleles) if g in self.hets
            ]
            if het_alleles:
                allele = random.choice(het_alleles)
            else:
                allele = random.choice(alleles)
        else:
            allele = 'ref'
        # Add reads to trees
        for read in reads:
            interval = intervaltree.Interval(
                read.reference_start, read.reference_end, read.query_name
            )
            if allele == 'ref':
                self.ref.add(interval)
            elif allele == 'alt':
                self.alt.add(interval)
            elif allele == 'other':
                self.other.add(interval)

    def get_allele_counts(self, start, end):
        # Loop through interval trees
        count_list = []
        for tree in (self.ref, self.alt, self.other):
            # Count unique reads overlapping inter
            reads = set()
            for interval in tree.overlap(start, end):
                reads.add(interval.data)
            read_counts = len(reads)
            count_list.append(read_counts)
        # Create tuple and return
        count_tuple = self.counts(*count_list)
        return(count_tuple)


class generate_counts(object):

    def __init__(
        self, bam, vcf, sample, window, paired_end, partial
    ):
        # Store input arguments
        self.bam_path = bam
        self.vcf_path = vcf
        self.sample = sample
        self.window = window
        self.paired_end = paired_end
        self.partial = partial
        # Create VarTree object
        self.vartree = vartree.VarTree(
            path=self.vcf_path, sample=self.sample
        )
        # Extract data from BAM file
        with pysam.AlignmentFile(self.bam_path) as bam:
            # Check pairing
            for read in bam.head(100):
                if self.paired_end:
                    assert(read.is_paired)
                else:
                    assert(not read.is_paired)
            # Get total counts
            if self.paired_end:
                self.total = bam.mapped // 2
            else:
                self.total = bam.mapped
            # Get chromosome lengths
            self.chrom_lengths = {}
            for chrom in bam.references:
                self.chrom_lengths[chrom] = bam.get_reference_length(chrom)

    def bam_generator(self, chromosome):
        # Set processing variables
        unpaired = {}
        # Loop thorugh chromosome reads
        with pysam.AlignmentFile(self.bam_path) as bam:
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
            # Check that all reads have been processed
            assert(len(unpaired) == 0)

    def get_genotype(self, genotype):
        ''' Function designed for multiple alternative alleles'''
        # Set all alternative alleles to 1 and return
        genotype = tuple(
            [None if g is None else min(g, 1) for g in genotype]
        )
        return(genotype)

    def get_het_prob(self, phred_probs, n_alt):
        ''' Function designed for multiple alternative alleles'''
        if None in phred_probs:
            prob_het = None
        else:
            het_indices = [int(i * ((i + 1) / 2)) for i in range(1, 1 + n_alt)]
            phred_prob_hets = [phred_probs[i] for i in het_indices]
            prob_het = sum([10 ** ((-x) / 10) for x in phred_prob_hets])
            prob_het = min(1, prob_het)
        return(prob_het)

    def get_read_variants(
        self, read_list
    ):
        # Get read variants for single end reads
        if len(read_list) == 1:
            read_variants = self.vartree.get_read_variants(
                read=read_list[0], partial=self.partial
            )
            variant_list = [read_variants]
        # Get read variants for paired end reads
        elif len(read_list) == 2:
            read1_variants, read2_variants, identical = (
                self.vartree.get_paired_read_variants(
                    read1=read_list[0], read2=read_list[1],
                    partial=self.partial
                )
            )
            assert(identical)
            variant_list = [read1_variants, read2_variants]
        # Loop trhough read sequences and variants
        read_variants = {}
        for read, variants in zip(read_list, variant_list):
            for position, variant in variants.items():
                # Match read sequence to ref, alt, other
                read_sequence = read.query_sequence[variant.start:variant.end]
                if read_sequence == variant.ref:
                    allele = 'ref'
                elif read_sequence in variant.alts:
                    allele = 'alt'
                else:
                    allele = 'other'
                # Get variant id
                variant_id = (position, variant.ref, variant.alts)
                read_variants[variant_id] = allele
        return(read_variants)

    def add_window_counts(self, metrics, coverage, chromosome):
        # Get chromosome length
        chromosome_length = self.chrom_lengths[chromosome]
        # Loop through variants
        for variant_id in metrics.keys():
            # Extract variant data
            start, ref = variant_id[0:2]
            end = start + len(ref)
            # Get window counts
            window_start = max(start - self.window, 0)
            window_end = min(end + self.window, chromosome_length)
            window_counts = coverage.get_allele_counts(
                window_start, window_end
            )
            # Add counts to metrics
            window_metrics = {
                'window_start': window_start,
                'window_end': window_end,
                'ref_window': window_counts.ref,
                'alt_window': window_counts.alt,
                'other_window': window_counts.other
            }
            metrics[variant_id].update(window_metrics)
        # Return updated metrics
        return(metrics)

    def write_metrics(self, chromosome, metrics, coverage, outfile):
        # Extract variants and sort by positions
        variant_ids = list(metrics.keys())
        variant_ids = sorted(variant_ids, key=lambda x: x[0])
        # Get counts for each variant
        for variant_id in variant_ids:
            variant_metrics = metrics[variant_id]
            # Get variant id metrics
            start, ref, alts = variant_id
            end = start + len(ref)
            alts = ','.join(alts)
            # Process genotype and probabilities
            genotype = '/'.join(map(str, variant_metrics['genotype']))
            genotype = genotype.replace('None', '.')
            if variant_metrics['hetprob'] is None:
                hetprob = '.'
            else:
                hetprob = '{:.2f}'.format(variant_metrics['hetprob'])
            # Create output line
            out_line = (
                '{chromosome}\t{start}\t{end}\t{window_start}\t{window_end}\t'
                '{ref}\t{alts}\t{sample}\t{genotype}\t{hetprob}\t'
                '{ref_count}\t{alt_count}\t{other_count}\t{ref_count_window}\t'
                '{alt_count_window}\t{other_count_window}\t{total_count}\n'
            ).format(
                chromosome=chromosome,
                start=start + 1,
                end=end,
                window_start=variant_metrics['window_start'] + 1,
                window_end=variant_metrics['window_end'],
                ref=ref,
                alts=alts,
                sample=self.sample,
                genotype=genotype,
                hetprob=hetprob,
                ref_count=variant_metrics['ref'],
                alt_count=variant_metrics['alt'],
                other_count=variant_metrics['other'],
                ref_count_window=variant_metrics['ref_window'],
                alt_count_window=variant_metrics['alt_window'],
                other_count_window=variant_metrics['other_window'],
                total_count=self.total
            )
            outfile.write(out_line)

    def process_chromosome_variants(self, chromosome, outfile):
        # Create coverage object and read vcf
        coverage = bam_coverage()
        self.vartree.read_vcf(chromosome)
        # Populate metrics dictionary with default values
        metrics = {}
        for interval in self.vartree.tree:
            # Process genotype and probabilites
            genotype = self.get_genotype(interval.data.genotype)
            hetprob = self.get_het_prob(
                interval.data.probs, len(interval.data.alts)
            )
            # Generate key for variant
            variant_id = (
                interval.begin, interval.data.ref, interval.data.alts
            )
            # Add variant to dictionary
            metrics[variant_id] = {
                'ref': 0, 'alt': 0, 'other': 0, 'genotype': genotype,
                'hetprob': hetprob
            }
        # Get variant metrics
        if len(metrics) > 0:
            for read_list in self.bam_generator(chromosome):
                # Extract read variants and their associayed metrics
                read_variants = self.get_read_variants(read_list=read_list)
                variant_ids = list(read_variants.keys())
                alleles = [read_variants[i] for i in variant_ids]
                genotypes = [metrics[i]['genotype'] for i in variant_ids]
                # Add allele counts
                for variant_id, allele in zip(variant_ids, alleles):
                    metrics[variant_id][allele] += 1
                # Add coverage
                coverage.add_reads(
                    reads=read_list, genotypes=genotypes, alleles=alleles
                )
            # Add window counts to metrics
            metrics = self.add_window_counts(
                metrics, coverage, chromosome
            )
            # Print metrics to file
            self.write_metrics(
                chromosome=chromosome, metrics=metrics, coverage=coverage,
                outfile=outfile
            )

    def process_all_variants(self, outpath):
        # Open outfile
        with open(outpath, 'wt') as outfile:
            for chromosome in self.vartree.chromosomes:
                self.process_chromosome_variants(chromosome, outfile)


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
            "Coordinate-sorted input BAM file."
        )
    )
    parser.add_argument(
        "--vcf", required=True, help=(
            "Coordinate sorted, and tabix indexed, VCF file."
        )
    )
    parser.add_argument(
        "--outprefix", required=True, help=(
            "Prefix of output files."
        )
    )
    parser.add_argument(
        "--sample", required=True, help=(
            "Sample for which to extract genotype from VCF (default=None)"
        )
    )
    parser.add_argument(
        "--window", required=True, type=int, help=(
            "Size of window around variants in which to calculate coverage"
        )
    )
    parser.add_argument(
        "--paired_end", action='store_true', help=(
            "BAM file contains paired end alignments (default=False)"
        )
    )
    args = parser.parse_args()
    # Generate output files
    outfiles = {
        'initial': args.outprefix + '.variant_counts.txt',
        'compressed': args.outprefix + '.variant_counts.txt.gz'
    }
    # Generate counts
    count_generator = generate_counts(
        bam=args.bam, vcf=args.vcf, sample=args.sample,
        paired_end=args.paired_end, window=200, partial=False
    )
    count_generator.process_all_variants(outfiles['initial'])
    # Compress and index file
    pysam.tabix_compress(
        filename_in=outfiles['initial'], filename_out=outfiles['compressed'],
        force=True
    )
    pysam.tabix_index(
        filename=outfiles['compressed'], seq_col=0, start_col=1, end_col=1,
        force=True
    )
