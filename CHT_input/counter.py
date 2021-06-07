import collections
import intervaltree
import itertools
import math
import pysam
import sys


def add_logs(
    loga, logb
):
    sum_logs = max(loga, logb) + math.log(1 + math.exp(-abs(loga - logb)))
    return(sum_logs)


def calculate_posterior_hetp(
    prior, ref, alt, error=0.01
):
    # Adjust prior
    prior = min(0.99, prior)
    # Calculate genotype likelihoods
    badlike = add_logs(
        (math.log(error) * ref) + (math.log(1 - error) * alt),
        (math.log(1 - error) * ref) + (math.log(error) * alt)
    )
    goodlike = (math.log(0.5) * ref) + (math.log(0.5) * alt)
    # avoid overflow (very close to 1.0)
    if goodlike - badlike > 40:
        hetp_post = 1.0
    # Or calculate posterior
    else:
        hetp_post = (
            prior * math.exp(goodlike - badlike) /
            (prior * math.exp(goodlike - badlike) + (1.0 - prior))
        )
    # Return posterior
    return(hetp_post)


VariantTuple = collections.namedtuple(
    '_Variant', [
        'chrom', 'start', 'end', 'id', 'ref', 'alt', 'haplotype',
        'genotype', 'het_prob', 'ref_as_count', 'alt_as_count',
        'other_as_count'
    ]
)


class IndividualVariant(VariantTuple):

    @classmethod
    def from_line(cls, line, adjhetprob=None):
        # Split line data and extract values
        line_data = line.strip().split('\t')
        assert(len(line_data) == 15)
        # Extract variant description
        chrom = line_data[0]
        position = int(line_data[1])
        id, ref, alt = line_data[2:5]
        start = position - 1
        end = start + len(ref)
        # Extract genotype information
        if 'NA' in line_data[5:9]:
            assert(line_data[5:9] == ['NA', 'NA', 'NA', 'NA'])
            haplotype = None
            het_prob = None
            genotype = None
        else:
            haplotype = line_data[5]
            het_prob, alt_prob = map(float, line_data[7:9])
            genotype = het_prob + (2 * alt_prob)
        # Get counts
        ref_as_count, alt_as_count, other_as_count = map(
            int, line_data[9:12]
        )
        ref_total_count, alt_total_count, other_total_count = map(
            int, line_data[12:15]
        )
        # Calculate adjusted heterozygous possibility
        if het_prob is not None and adjhetprob is not None:
            if adjhetprob == 'as':
                het_prob = calculate_posterior_hetp(
                    het_prob, ref=ref_as_count, alt=alt_as_count
                )
            elif adjhetprob == 'total':
                het_prob = calculate_posterior_hetp(
                    het_prob, ref=ref_total_count, alt=alt_total_count
                )
            else:
                raise ValueError('adjhetprob variable not recognised')
        # Generate hash
        new_variant = cls(
            chrom=chrom, start=start, end=end, id=id, ref=ref, alt=alt,
            haplotype=haplotype, genotype=genotype, het_prob=het_prob,
            ref_as_count=ref_as_count, alt_as_count=alt_as_count,
            other_as_count=other_as_count
        )
        return(new_variant)

    def __hash__(
        self
    ):
        return(self.hash)

    def __eq__(
        self, other
    ):
        if self.__class__ == other.__class__ and self.hash == other.hash:
            return(True)
        return(False)

    def __ne__(
        self, other
    ):
        return(not self.__eq__(other))


class CountTree(object):

    def __init__(self, path, adjhetprob):
        # Store initial arguments
        self.path = path
        self.adjhetprob = adjhetprob
        # Get vcf chromosomes and check samples
        with pysam.TabixFile(self.path) as count_file:
            self.chromosomes = count_file.contigs
        # Create empty slots for data
        self.previous_chromosomes = []
        self.current_chromosome = None
        self.variants = None
        self.regions = None
        self.tree = None
        # Tuples containing possible haplotypes
        self.haplotypes = set(['0|0', '0|1', '1|0', '1|1'])
        self.heterozygotes = set(['0|1', '1|0'])

    def read_counts(self, chromosome):
        '''Retrieves variant counts from file and creates an intervaltree
        containing variant locations for rapid retreival by location

        Parameters
        ----------
        chromosome:
            Name of chromosome for which to get variants

        Raises
        ------
        ValueError:
            If chromosome has been requested previously. This ensures that
            chromosomes are processed sequentially which reduces run time.
        '''
        # Check chromosome is new
        if chromosome in self.previous_chromosomes:
            raise ValueError('input must be sorted by chromosome')
        self.previous_chromosomes.append(chromosome)
        self.current_chromosome = chromosome
        # Create empty iterator for missing chromosomes or...
        if chromosome not in self.chromosomes:
            warning = "WARNING: {} not in VCF header\n".format(chromosome)
            sys.stderr.write(warning)
            chrom_iter = []
        # ...create iterator without sample data
        else:
            count_file = pysam.TabixFile(self.path)
            chrom_iter = count_file.fetch(reference=chromosome)
        # Set variables to process data
        self.variants = []
        intervals = []
        self.region_variants = {}
        self.test_variants = {}
        # Loop thorugh chromosome variants in vcf file
        for index, line in enumerate(chrom_iter):
            # Create variant and store
            variant = IndividualVariant.from_line(line, self.adjhetprob)
            self.variants.append(variant)
            # Create intervaltree interval and add to list
            interval = intervaltree.Interval(
                variant.start, variant.end, index
            )
            intervals.append(interval)
        # Create intervaltree IntervalTree from list of intervals
        self.tree = intervaltree.IntervalTree(intervals)

    def get_variants(self, starts, ends):
        '''Retrieves variants overlapping the specified interval

        Parameters
        ----------
        starts:
            an iterator of 0-based interval starts
        ends:
            an iterator of 0-based interval ends

        Returns
        -------
        variants: list
            A list of unique IndividualVariant object containing variant data
        '''
        # Get index of all variant overlapping data
        indices = set()
        for start, end in itertools.zip_longest(starts, ends):
            for interval in self.tree.envelop(start, end):
                indices.add(interval.data)
        # Extract variant for each index
        variants = [self.variants[i] for i in sorted(list(indices))]
        return(variants)

    def generate_test_string(self, start, end, ref, alt):
        '''Retreives specified test variant and generates a string suitable
        for insertion into a CHT input file.

        Parameters
        ----------
        start int:
            start position of test variant
        end int:
            end position of test variant
        ref str:
            refrenence allele of test variant
        alt str:
            alternative allele of test variant

        Returns
        -------
        region_str: str
            A string containing metrics for the test variant for the combined
            halotype test.

        Raises:
        -------
        AssertionError:
            If a single variant matching reference and alternative alleles
            is not found.
        '''
        # Get variant
        test_variants = self.get_variants((start,), (end,))
        assert(len(test_variants) == 1)
        test_variant = test_variants[0]
        assert(test_variant.ref == ref)
        assert(test_variant.alt == alt)
        # Generate string
        test_list = [
            self.current_chromosome,
            str(test_variant.start + 1),
            test_variant.id,
            test_variant.ref,
            test_variant.alt,
            '{:.2f}'.format(test_variant.genotype),
            test_variant.haplotype
        ]
        test_str = ' '.join(test_list)
        return(test_str)

    def get_test_string(self, start, end, ref, alt):
        '''Will retreive test variant string if test variant has been
        observed previously or will generate and store test variant string
        if test variant is novel.

        Parameters
        ----------
        start int:
            start position of test variant
        end int:
            end position of test variant
        ref str:
            refrenence allele of test variant
        alt str:
            alternative allele of test variant

        Returns
        -------
        region_str: str
            A string containing metrics for the test variant for the combined
            halotype test.

        Raises:
        -------
        AssertionError:
            If a single variant matching reference and alternative alleles
            is not found.
        '''
        test_tuple = (start, end, ref, alt)
        if test_tuple in self.test_variants:
            test_str = self.test_variants[test_tuple]
        else:
            test_str = self.generate_test_string(
                start=start, end=end, ref=ref, alt=alt
            )
            self.test_variants[test_tuple] = test_str
        return(test_str)

    def generate_region_string(self, test_haplotype, starts, ends):
        '''Will retreive variants overlapping the specified interval and format
        them into a string suitable for insertion into a CHT input file

        Parameters
        ----------
        test_haplotype:
            haplotype of the test snp (one of 0|0, 0|1, 1|0 or 1|1)
        starts:
            an iterator of 0-based interval starts
        ends:
            an iterator of 0-based interval ends

        Returns
        -------
        region_str: str
            A list of unique IndividualVariant object containing variant data
        '''
        # Check arguments
        assert(test_haplotype in self.haplotypes)
        # Get region variants
        region_variants = self.get_variants(starts, ends)
        # Process heterozygotic test variants...
        if test_haplotype in self.heterozygotes:
            # Set zero counts
            ref_hap_counts = []
            alt_hap_counts = []
            other_hap_counts = []
            # Loop through region variants and check haplotype
            for variant in region_variants:
                assert(variant.haplotype in self.haplotypes)
                # Extract counts for heterozygotic variants...
                if variant.haplotype in self.heterozygotes:
                    if variant.haplotype == test_haplotype:
                        ref_hap_counts.append(variant.ref_as_count)
                        alt_hap_counts.append(variant.alt_as_count)
                        other_hap_counts.append(variant.other_as_count)
                    # Get counts for discordant heterozygous haplotypes
                    else:
                        ref_hap_counts.append(variant.alt_as_count)
                        alt_hap_counts.append(variant.ref_as_count)
                        other_hap_counts.append(variant.other_as_count)
                # or set counts to zero
                else:
                    ref_hap_counts.append(0)
                    alt_hap_counts.append(0)
                    other_hap_counts.append(0)
        # or set counts to zero for homozygotes
        else:
            ref_hap_counts = [0] * len(region_variants)
            alt_hap_counts = [0] * len(region_variants)
            other_hap_counts = [0] * len(region_variants)
        # Create string and return
        region_positions = [v.start + 1 for v in region_variants]
        region_hetprobs = [v.het_prob for v in region_variants]
        region_linkage = ['1.00' for v in region_variants]
        # Merge strings
        region_list = [
            ';'.join([str(s + 1) for s in starts]),
            ';'.join([str(e) for e in ends]),
            ';'.join([str(p) for p in region_positions]),
            ';'.join(['{:.2f}'.format(h) for h in region_hetprobs]),
            ';'.join([k for k in region_linkage]),
            ';'.join([str(r) for r in ref_hap_counts]),
            ';'.join([str(a) for a in alt_hap_counts]),
            ';'.join([str(o) for o in other_hap_counts])
        ]
        region_str = ' '.join(region_list)
        return(region_str)

    def get_region_string(self, test_haplotype, starts, ends):
        '''Will retreive region string if test haplotype and region have been
        observed previously or will generate and store region string if test
        haplotype and region are novel.

        Parameters
        ----------
        test_haplotype:
            haplotype of the test snp (one of 0|0, 0|1, 1|0 or 1|1)
        starts:
            an iterator of 0-based interval starts
        ends:
            an iterator of 0-based interval ends

        Returns
        -------
        region_str: str
            A list of unique IndividualVariant object containing variant data
        '''
        # Generate region tuple
        region_tuple = (test_haplotype, tuple(starts), tuple(ends))
        # Extract region string if tuple has been observed previously...
        if region_tuple in self.region_variants:
            region_str = self.region_variants[region_tuple]
        # or generate and store region string if tuple is novel
        else:
            region_str = self.generate_region_string(
                test_haplotype=test_haplotype, starts=starts, ends=ends
            )
            self.region_variants[region_tuple] = region_str
        # Return str
        return(region_str)


class BamCounts(object):

    def __init__(self, path):
        # Store paths and open BAM files
        self.path = path
        self.bam_file = pysam.AlignmentFile(self.path)
        # Get total counts
        self.total = self.bam_file.mapped
        # Get chromosome lengths
        self.chrom_lengths = {}
        for chrom in self.bam_file.references:
            chrom_len = self.bam_file.get_reference_length(chrom)
            self.chrom_lengths[chrom] = chrom_len
        # Generate variable to store counts
        self.region_counts = {}

    def generate_count_string(self, chrom, starts, ends):
        counts = 0
        for start, end in zip(starts, ends):
            counts += self.bam_file.count(chrom, start, end)
        count_list = [
            str(counts),
            str(self.total)
        ]
        count_str = ' '.join(count_list)
        return(count_str)

    def get_count_string(self, chrom, starts, ends):
        # Generate region tuple
        region_tuple = (chrom, tuple(starts), tuple(ends))
        if region_tuple in self.region_counts:
            region_str = self.region_counts[region_tuple]
        else:
            region_str = self.generate_count_string(
                chrom=chrom, starts=starts, ends=ends
            )
            self.region_counts[region_tuple] = region_str
        return(region_str)

    def close(self):
        self.bam_file.close()


x = CountTree(
    '/Users/rabinowi/wasp/test_data/alignments/paired_end.variant_counts.txt.gz',
    adjhetprob='total'
)
x.read_counts('chr4')
region_str = x.get_region_string('0|1', [0], [1000])
print(region_str)
region_str = x.get_region_string('0|1', [0], [1000])
print(region_str)
test_str = x.get_test_string(337, 400, 'T', 'A')
print(test_str)
test_str = x.get_test_string(337, 400, 'T', 'A')
print(test_str)
y = BamCounts(
    '/Users/rabinowi/wasp/test_data/alignments/paired_end.filtered.rmdup.bam'
)
bam_str = y.get_count_string('chr2L', (0,), (10000,))
print(bam_str)
