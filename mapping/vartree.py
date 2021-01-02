import collections
import bisect
import intervaltree
import pysam
import re
import sys


class VarTable(object):

    def __init__(self, path):
        # Set path to vcf file
        self.vcf = path
        # Create variant named tuple
        self.variant = collections.namedtuple(
            'variant', ['start', 'end', 'ref', 'alt']
        )
        # Create allele regx
        ref_characters = ('A', 'C', 'G', 'T')
        ref_regx_string = '^(' + '|'.join(ref_characters) + ')+$'
        self.ref_regx = re.compile(ref_regx_string)
        alt_characters = ('A', 'C', 'G', 'T', '\\*', ',')
        alt_regx_string = '^(' + '|'.join(alt_characters) + ')*$'
        self.alt_regx = re.compile(alt_regx_string)
        # Clear allele data
        self.clear()

    def clear(self):
        self.chromosome = None
        self.alleles = None

    def read_vcf(self, chromosome):
        """read in SNPs and indels from text input file"""
        # Clear any prexisting alleles
        self.clear()
        # Extract intervaltree intervals for chromomosome variants in VCF
        interval_list = []
        with pysam.TabixFile(self.vcf) as vcf:
            # Create iterator for chromosome variants
            try:
                vcf_iterator = vcf.fetch(chromosome)
            except ValueError:
                sys.stderr.write(
                    "WARNING: cannot find chromosome {} in VCF\n".format(
                        chromosome
                    )
                )
                vcf_iterator = []
            # Read file and create list of intervaltree Intervals
            for line in vcf_iterator:
                # Extract variant data from line
                line_data = line.strip().split('\t')
                start = int(line_data[1])
                ref, alt = line_data[3:5]
                # Check position and alleles
                assert(start > 0)
                assert(self.ref_regx.match(ref) is not None)
                assert(self.alt_regx.match(alt) is not None)
                # Process position and allele
                start = int(start) - 1
                alt = alt.replace("*", "")
                end = start + len(ref)
                # Create intervaltree interval and add to list
                interval = intervaltree.Interval(
                    start, end, self.variant(start, end, ref, alt)
                )
                interval_list.append(interval)
        # Create intervaltree IntervalTree from list of intervals
        self.chromosome = chromosome
        self.alleles = intervaltree.IntervalTree(interval_list)

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
            intervals = self.alleles.overlap(start, end)
        # exclude intervals partially contained within interval
        else:
            intervals = self.alleles.envelop(start, end)
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
