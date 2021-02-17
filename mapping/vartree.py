import bisect
import collections
import intervaltree
import pysam
import sys


class VariantFunctions(object):

    def __hash__(
        self
    ):
        return(hash(self.id))

    def __eq__(
        self, other
    ):
        if (
            self.__class__ == other.__class__ and
            self.__hash__() == other.__hash__()
        ):
            return(True)
        return(False)

    def __ne__(
        self, other
    ):
        return(not self.__eq__(other))

    def is_snv(self):
        for allele in self.alleles:
            if len(allele) != 1:
                return(False)
        return(True)

    def is_biallelic(self):
        if len(self.alleles) == 2:
            return(True)
        return(False)

    def is_heterozygous(self):
        genotypes = set(self.genotype)
        if None not in genotypes and len(genotypes) == 2:
            return(True)
        return(False)

    def is_biallelic_heterozygous(self):
        if self.is_biallelic() and self.is_heterozygous():
            return(True)
        else:
            return(False)


class IntervalVariant(
    VariantFunctions, collections.namedtuple(
        'IntervalVariant', ['alleles', 'id', 'genotype', 'probs']
    )
):
    pass


class OverlappingVariant(
    VariantFunctions, collections.namedtuple(
        'OverlappingVariant', [
            'chrom', 'start', 'end', 'alleles', 'id', 'genotype', 'probs',
            'read_start', 'read_end', 'read_allele'
        ]
    )
):
    pass


class VarTree(object):

    def __init__(self, path, sample=None):
        # Store initial arguments
        self.path = path
        self.sample = sample
        # Get vcf chromosomes and check samples
        with pysam.VariantFile(self.path) as vcf:
            self.chromosomes = list(vcf.header.contigs)
            if self.sample:
                assert(self.sample in list(vcf.header.samples))
        # Set current chromosome
        self.current_chromosome = None
        # Create allele regx
        self.ref_set = set(['A', 'C', 'G', 'T'])
        self.alt_set = set(['A', 'C', 'G', 'T', '*'])
        # Create slot for tree
        self.tree = None

    def read_vcf(self, chromosome):
        """read in SNPs and indels from text input file"""
        # Set current chromosome
        self.current_chromosome = chromosome
        # Create empty iterator for missing chromosomes or...
        if chromosome not in self.chromosomes:
            sys.stderr.write(
                "WARNING: chromosome {} not in VCF header\n".format(
                    chromosome
                )
            )
            vcf = None
            chrom_iter = []
        # Create iterator without sample data or...
        elif self.sample is None:
            vcf = pysam.VariantFile(self.path, drop_samples=True)
            chrom_iter = vcf.fetch(contig=chromosome)
        # Create iterator for samples
        else:
            vcf = pysam.VariantFile(self.path, drop_samples=False)
            vcf.subset_samples([self.sample])
            chrom_iter = vcf.fetch(contig=chromosome)
        # Read file and create list of intervaltree Intervals
        interval_list = []
        for entry in chrom_iter:
            # Extract alleles and check
            alleles = entry.alleles
            for ref_base in alleles[0]:
                assert(ref_base in self.ref_set)
            for alt in alleles[1:]:
                for alt_base in alt:
                    assert(alt_base in self.alt_set)
            # Get genotype and associated probabilities
            if self.sample:
                # Get sample data and genotype
                sample_data = entry.samples[self.sample]
                genotype = sample_data['GT']
                # Set probabilities to None if genotype unknown or...
                if None in genotype:
                    probs = None
                # Calculate probabilities
                else:
                    if 'GL' in sample_data:
                        raw_prob = [10 ** p for p in sample_data['GL']]
                    elif 'PL' in sample_data:
                        raw_prob = [10 ** (-p / 10) for p in sample_data['PL']]
                    else:
                        raise KeyError('absent genotype probability')
                    # Normalise probability
                    probs = tuple([p / sum(raw_prob) for p in raw_prob])
            else:
                genotype = None
                probs = None
            # Get variant id
            if entry.id is None:
                variant_id = '{}_{}_{}'.format(
                    chromosome, entry.start + 1, '_'.join(alleles)
                )
            else:
                variant_id = entry.id
            # Create intervaltree interval and add to list
            interval_variant = IntervalVariant(
                alleles=alleles, id=variant_id, genotype=genotype,
                probs=probs
            )
            interval = intervaltree.Interval(
                entry.start, entry.stop, interval_variant
            )
            interval_list.append(interval)
        # Create intervaltree IntervalTree from list of intervals
        self.tree = intervaltree.IntervalTree(interval_list)

    def get_overlapping_variants(self, start, end):
        '''Retrieves variants overlapping the specified interval

        Parameters
        ----------
        start: int
            0-based interval start
        end: int
            0-based interval end

        Returns
        -------
        variants: list
            A list of IndividualVariant object containing variant data
        '''
        # include intervals partially contained within interval or...
        intervals = self.tree.overlap(start, end)
        variants = [
            OverlappingVariant(
                chrom=self.current_chromosome, start=i.begin, end=i.end,
                alleles=i.data.alleles, id=i.data.id, genotype=i.data.genotype,
                probs=i.data.probs, read_start=None, read_end=None,
                read_allele=None
            ) for i in intervals
        ]
        return(variants)

    def get_read_variants(
        self, read, partial
    ):
        '''Retrieves variants overlapping the supplied read

        Parameters
        ----------
        read: pysam.AlignedSegment
            An aligned read
        partial: bool
            Return variants partially overlapping read end and start

        Returns
        -------
        read_variants: list

        '''
        # Check chromosome
        if not self.current_chromosome == read.reference_name:
            raise ValueError('mismatched chromosomes')
        # Get all variants
        overlapping_variants = set()
        for start, end in read.get_blocks():
            for variant in self.get_overlapping_variants(start=start, end=end):
                overlapping_variants.add(variant)
        # Convert genomic variants into a list and sort
        overlapping_variants = list(overlapping_variants)
        overlapping_variants = sorted(
            overlapping_variants, key=lambda x: x.start
        )
        # Filter partially overlapping start and end
        if not partial:
            overlapping_variants = [
                ov for ov in overlapping_variants if
                ov.start >= read.reference_start and
                ov.end <= read.reference_end
            ]
        # Create read variants from genomic variants
        read_variants = []
        if overlapping_variants:
            # Extract read data
            positions = read.get_reference_positions(full_length=True)
            read_align_start = read.query_alignment_start
            read_align_end = read.query_alignment_end
            # Check all bases have assigned positions
            if None in positions[read_align_start:read_align_end]:
                # Generate cigar covering aligned segment of read
                read_cigar = ''.join([
                    str(operation) * length for
                    operation, length in read.cigartuples if
                    operation in (0, 1, 4, 6, 7, 8)
                ])
                assert(len(read_cigar) == len(positions))
                # Assign positions of insertion to previous bases
                for i in range(read_align_start, read_align_end):
                    if positions[i] is None:
                        assert(read_cigar[i] == '1')
                        positions[i] = positions[i - 1]
                # Raise error for aligned bases with no positions
                assert(None not in positions[read_align_start:read_align_end])
            # Loop through putative variant end get read variants
            for variant in overlapping_variants:
                # Add indices of variant within read
                read_start = bisect.bisect_left(
                    a=positions, x=variant.start, lo=read_align_start,
                    hi=read_align_end
                )
                read_end = bisect.bisect_left(
                    a=positions, x=variant.end, lo=read_align_start,
                    hi=read_align_end
                )
                read_allele = read.query_sequence[read_start:read_end]
                # Create read variant and store
                read_variant = variant._replace(
                    read_start=read_start, read_end=read_end,
                    read_allele=read_allele
                )
                read_variants.append(read_variant)
        # Return errors and passed variants
        return(read_variants)
