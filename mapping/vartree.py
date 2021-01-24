import collections
import bisect
import intervaltree
import pysam
import sys


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
        # Create variant named tuple
        self.variant = collections.namedtuple(
            'variant', ['start', 'end', 'ref', 'alts', 'genotype', 'probs']
        )
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
            # Extract ref and alt alleles and check
            ref, alts = entry.ref, entry.alts
            for base in ref:
                assert(base in self.ref_set)
            for alt in alts:
                for base in alt:
                    assert(base in self.alt_set)
            # Get genotype
            if self.sample:
                genotype = entry.samples[self.sample]['GT']
                probs = entry.samples[self.sample]['PL']
            else:
                genotype = None
                probs = None
            # Create intervaltree interval and add to list
            interval = intervaltree.Interval(
                entry.start, entry.stop, self.variant(
                    start=entry.start, end=entry.stop, ref=ref, alts=alts,
                    genotype=genotype, probs=probs
                )
            )
            interval_list.append(interval)
        # Create intervaltree IntervalTree from list of intervals
        self.tree = intervaltree.IntervalTree(interval_list)

    def get_overlapping_variants(self, start, end, partial):
        '''Retrieves variants overlapping the specified interval

        Parameters
        ----------
        start: int
            0-based interval start
        end: int
            0-based interval end
        partial: bool
            Return variants partially overlapping interval

        Returns
        -------
        variants: list
            A list of collections.namedtuple cotaining the genomic
            start and end positions of the variant as well as the
            reference and alternative alleles. List is sorted by
            variant start position.
        '''
        # Check partial argument and either...
        if not isinstance(partial, bool):
            raise TypeError('partial argument must be boolean')
        # include intervals partially contained within interval or...
        if partial:
            intervals = self.tree.overlap(start, end)
        # exclude intervals partially contained within interval
        else:
            intervals = self.tree.envelop(start, end)
        # Extract variants, sort and return
        variants = [x.data for x in intervals]
        variants = sorted(variants, key=lambda x: x.start)
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
            Return variants partially overlapping read

        Returns
        -------
        read_variants: collections.OrderedDict
            A dictionary of variants. The key is the 0-based position of
            the variant within the genome and the value is a named tuple
            cotaining the start and end positions of the variant within
            the read as well as the reference and alternative alleles.
            The dictionary is ordered by genomic position of the variants.
        '''
        # Set output variables
        read_variants = collections.OrderedDict()
        # Get variants
        variants = self.get_overlapping_variants(
            start=read.reference_start, end=read.reference_end, partial=partial
        )
        # Process supplied variants
        if len(variants) > 0:
            # Extract read data
            positions = read.get_reference_positions(full_length=True)
            read_start = read.query_alignment_start
            read_end = read.query_alignment_end
            # Check all bases have assigned positions
            if None in positions[read_start:read_end]:
                # Generate cigar covering aligned segment of read
                read_cigar = ''.join([
                    str(operation) * length for
                    operation, length in read.cigartuples if
                    operation in (0, 1, 4, 6, 7, 8)
                ])
                assert(len(read_cigar) == len(positions))
                # Assign positions of insertion to previous bases
                for i in range(read_start, read_end):
                    if positions[i] is None:
                        assert(read_cigar[i] == '1')
                        positions[i] = positions[i - 1]
                assert(None not in positions[read_start:read_end])
            # Loop through variants
            for variant in variants:
                # Get indices within reads
                start_index = bisect.bisect_left(
                    positions, variant.start, lo=read_start, hi=read_end
                )
                end_index = bisect.bisect_left(
                    positions, variant.end, lo=read_start, hi=read_end
                )
                # Store varaints mapped to reads
                mapped_variant = variant._replace(
                    start=start_index, end=end_index
                )
                read_variants[variant.start] = mapped_variant
        # Return errors and passed variants
        return(read_variants)

    def get_paired_read_variants(
        self, read1, read2, partial
    ):
        '''Retrieves variants overlapping the supplied read

        Parameters
        ----------
        read1: pysam.AlignedSegment
            First read of an aligned pair
        read2: pysam.AlignedSegment
            Second read of an aligned pair
        partial: bool
            Return variants partially overlapping reads

        Returns
        -------
        read1_variants: collections.OrderedDict
            A dictionary of variants for read1. The key is the 0-based
            position of the variant within the genome and the value is a
            named tuple cotaining the start and end positions of the variant
            within the read as well as the reference and alternative alleles.
            The dictionary is ordered by genomic position of the variants.
        read2_variants: collections.OrderedDict
            Same as read1_variants but for read2
        identical_sequence: bool
            A boolean indicating if the read sequences are identical at
            variants common to read1 and read2 match.
        '''
        # Get variants for the reads
        read1_variants = self.get_read_variants(read=read1, partial=partial)
        read2_variants = self.get_read_variants(read=read2, partial=partial)
        # Check common variants are identical
        common_variants = read1_variants.keys() & read2_variants.keys()
        identical_variants = True
        for variant_position in common_variants:
            # Get read sequences for each read
            read1_variant = read1_variants[variant_position]
            read2_variant = read2_variants[variant_position]
            read1_seq = read1.query_sequence[
                read1_variant.start:read1_variant.end
            ]
            read2_seq = read2.query_sequence[
                read2_variant.start:read2_variant.end
            ]
            # Deterimine if the sequences are identical
            if read1_seq != read2_seq:
                identical_variants = False
                break
        # Return variant
        return(read1_variants, read2_variants, identical_variants)
