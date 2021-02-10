import bisect
import collections
import copy
import intervaltree
import pysam
import sys


VariantTuple = collections.namedtuple(
    'VariantTuple', [
        'chrom', 'start', 'end', 'alleles', 'id', 'gt', 'pl',
        'read_start', 'read_end', 'read_allele', 'hash'
    ]
)


class IndividualVariant(VariantTuple):

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
        if None not in self.gt and len(set(self.gt)) == 2:
            return(True)
        return(False)

    def is_biallelic_heterozygous(self):
        if self.is_biallelic() and self.is_heterozygous():
            return(True)
        else:
            return(False)


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
            # Get genotype
            if self.sample:
                gt = entry.samples[self.sample]['GT']
                pl = entry.samples[self.sample]['PL']
            else:
                gt = (None,)
                pl = (None,)
            # Get variant id
            if entry.id is None:
                variant_id = '{}_{}_{}'.format(
                    chromosome, entry.start + 1, '_'.join(alleles)
                )
            else:
                variant_id = entry.id
            # Create intervaltree interval and add to list
            new_variant = IndividualVariant(
                chrom=chromosome, start=entry.start, end=entry.stop,
                alleles=alleles, id=variant_id, gt=gt, pl=pl, read_start=None,
                read_end=None, read_allele=None,
                hash=hash((chromosome, entry.start, alleles))
            )
            new_interval = intervaltree.Interval(
                entry.start, entry.stop, new_variant
            )
            interval_list.append(new_interval)
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
        variants = [x.data for x in intervals]
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
        # Get all variants
        genomic_variants = set()
        for start, end in read.get_blocks():
            for variant in self.get_overlapping_variants(start=start, end=end):
                genomic_variants.add(variant)
        # Convert genomic variants into a list and sort
        genomic_variants = list(genomic_variants)
        genomic_variants = sorted(genomic_variants, key=lambda x: x.start)
        # Filter partially overlapping start and end
        if not partial:
            genomic_variants = [
                gv for gv in genomic_variants if
                gv.start >= read.reference_start and
                gv.end <= read.reference_end
            ]
        # Create read variants from genomic variants
        read_variants = []
        if genomic_variants:
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
                assert(None not in positions[read_align_start:read_align_end])
            # Loop through putative variant end get read variants
            for variant in genomic_variants:
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
