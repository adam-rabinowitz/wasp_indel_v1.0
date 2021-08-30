import bisect
import collections
import intervaltree
import itertools
import pysam
import sys


class Variant(
    collections.namedtuple(
        'Variant', [
            'chrom', 'start', 'end', 'alleles', 'id', 'genotypes', 'probs',
            'read_start', 'read_end', 'read_allele'
        ]
    )
):

    def __hash__(self):
        return(hash(self.id))

    def __eq__(self, other):
        if (
            self.__class__ == other.__class__ and
            self.__hash__() == other.__hash__()
        ):
            return(True)
        return(False)

    def __ne__(self, other):
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

    def is_heterozygous(self, sample):
        genotypes = set(self.genotypes[sample])
        if None not in genotypes and len(genotypes) == 2:
            return(True)
        return(False)

    def is_biallelic_snv(self):
        if self.is_biallelic() and self.is_snv():
            return(True)
        return(False)

    def is_biallelic_heterozygous(self, sample):
        if self.is_biallelic() and self.is_heterozygous(sample):
            return(True)
        return(False)


class VarTree(object):

    def __init__(self, path, samples=None, check_phase=True):
        # Store initial arguments
        self.path = path
        self.samples = samples if samples else []
        self.check_phase = check_phase
        # Get vcf chromosomes and check samples
        with pysam.VariantFile(self.path) as vcf:
            self.chromosomes = list(vcf.header.contigs)
            for sample in self.samples:
                assert(sample in list(vcf.header.samples))
        # Create allele regx
        self.ref_set = set(['A', 'C', 'G', 'T'])
        self.alt_set = set(['A', 'C', 'G', 'T', '*'])
        # Create empty slots for variant data
        self.current_chromosome = None
        self.variant_tree = None
        self.variants = None

    def read_vcf(self, chromosome):
        """read in SNPs and indels from text input file"""
        # Store chromosome and create variables to process data
        self.current_chromosome = chromosome
        self.variants = []
        intervals = []
        variant_ids = set()
        entry_index = 0
        # Create empty iterator for missing chromosomes or...
        if chromosome not in self.chromosomes:
            warning = "WARNING: {} not in VCF header\n".format(chromosome)
            sys.stderr.write(warning)
            vcf = None
            chrom_iter = []
        # ...create iterator with sample data or...
        elif self.samples:
            vcf = pysam.VariantFile(self.path, drop_samples=False)
            vcf.subset_samples(self.samples)
            chrom_iter = vcf.fetch(contig=chromosome)
        # ...create iterator without sample data
        else:
            vcf = pysam.VariantFile(self.path, drop_samples=True)
            chrom_iter = vcf.fetch(contig=chromosome)
        # Loop thorugh chromosome variants in vcf file
        for entry in chrom_iter:
            # Extract alleles and check bases
            alleles = entry.alleles
            for ref_base in alleles[0]:
                assert(ref_base in self.ref_set)
            for alt in alleles[1:]:
                for alt_base in alt:
                    assert(alt_base in self.alt_set)
            # Check for non reference bases
            if self.samples:
                # Get allele indices for all genotypes
                gt_values = [
                    allele_index for sample in self.samples
                    for allele_index in entry.samples[sample]['GT']
                ]
                # Skip variants where all alleles are reference
                if (max(gt_values) == 0):
                    continue
            # Get genotype probabilities for samples
            genotypes = {}
            probs = {}
            # Populate genotype probabilities for named samples
            for sample in self.samples:
                # Get sample data, check phase and get genotype
                sample_data = entry.samples[sample]
                if self.check_phase and not sample_data.phased:
                    raise ValueError('unphased variant')
                sample_genotype = sample_data['GT']
                # Set probabilities to None if genotype unknown or...
                if None in sample_genotype:
                    sample_probs = (None,)
                # Calculate probabilities
                else:
                    if 'GL' in sample_data:
                        sample_gl = sample_data['GL']
                        raw_probs = [10 ** p for p in sample_gl]
                    elif 'PL' in sample_data:
                        sample_pl = sample_data['PL']
                        raw_probs = [10 ** (-p / 10) for p in sample_pl]
                    else:
                        raise KeyError('absent genotype probability')
                    # Normalise probability
                    sample_probs = tuple(
                        [p / sum(raw_probs) for p in raw_probs]
                    )
                # Store genotype and probs
                genotypes[sample] = sample_genotype
                probs[sample] = sample_probs
            # Get variant id
            if entry.id is None:
                variant_id = '{}_{}_{}'.format(
                    chromosome, entry.start + 1, '_'.join(alleles)
                )
            else:
                variant_id = entry.id
            # Check id is unique
            assert(variant_id not in variant_ids)
            variant_ids.add(variant_id)
            # Create variant and add to list
            variant = Variant(
                chrom=chromosome, start=entry.start, end=entry.stop,
                alleles=alleles, id=variant_id, genotypes=genotypes,
                probs=probs, read_start=None, read_end=None,
                read_allele=None
            )
            self.variants.append(variant)
            # Create intervaltree interval and add to list
            interval = intervaltree.Interval(
                entry.start, entry.stop, entry_index
            )
            intervals.append(interval)
            # Increase index
            entry_index += 1
        # Create intervaltree IntervalTree from list of intervals
        assert(len(self.variants) == entry_index)
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
            for interval in self.tree.overlap(start, end):
                indices.add(interval.data)
        # Extract variant for each index
        variants = [self.variants[i] for i in sorted(list(indices))]
        return(variants)

    def get_read_positions(self, read):
        # Get initial positions
        read_align_start = read.query_alignment_start
        read_align_end = read.query_alignment_end
        positions = read.get_reference_positions(full_length=True)
        # Attempt to add positions where None is present
        if None in positions[read_align_start:read_align_end]:
            # Generate cigar covering aligned segment of read
            cigar = ''.join([
                str(operation) * length for
                operation, length in read.cigartuples if
                operation in (0, 1, 4, 7, 8)
            ])
            assert(len(cigar) == len(positions))
            # Assign positions of insertion to previous bases
            for i in range(read_align_start, read_align_end):
                if positions[i] is None:
                    assert(cigar[i] == '1')
                    positions[i] = positions[i - 1]
        # Check all positions have now been assigned
        if None in positions[read_align_start:read_align_end]:
            message = 'Abnormal Alignment: {} {}'.format(
                read.query_name, read.cigarstring
            )
            raise ValueError(message)
        # Return alignment positions
        return(read_align_start, read_align_end, positions)

    def get_read_variants(self, read, partial):
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
        overlapping_variants = self.get_variants(
            *zip(*read.get_blocks())
        )
        # Filter partially overlapping start and end
        if not partial:
            overlapping_variants = [
                ov for ov in overlapping_variants if
                ov.start >= read.reference_start and
                ov.end <= read.reference_end
            ]
        # Add read data to variants
        read_variants = []
        if overlapping_variants:
            # Get positions of bases
            read_align_start, read_align_end, positions = (
                self.get_read_positions(read)
            )
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
        # Return read variants
        return(read_variants)
