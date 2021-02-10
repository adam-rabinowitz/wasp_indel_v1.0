MetricsTuple = collections.namedtuple(
    '_Variant', [
        'chrom', 'position', 'id', 'ref', 'alt', 'haplotype', 'ref_prob',
        'het_prob', 'alt_prob', 'ref_count', 'alt_count', 'other_count',
        'hash'
    ]
)


class VariantMetrics(MetricTuple):

    @classmethod
    def from_line(cls, line):
        # Split line data and extract values
        line_data = line.strip().split('\t')
        assert(len(line_data) == 12)
        # Extract test variant
        chrom = line_data[0]
        position = int(line_data[1])
        id, ref, alt = line_data[2:5]
        # Extract genotype information
        if 'NA' in line_data[5:9]:
            assert(line_data[5:9] == ['NA', 'NA', 'NA', 'NA'])
            haplotype = None
            ref_prob, het_prob, alt_prob = None, None, None
        else:
            haplotype = line_data[5]
            ref_prob, het_prob, alt_prob = map(float, line_data[6:9])
        # Get counts
        ref_count, alt_count, other_count = map(int, line_data[9:12])
        # Generate hash
        new_variant = cls(
            chrom=chrom, position=position, id=id, ref=ref, alt=alt,
            haplotype=haplotype, ref_prob=ref_prob, het_prob=het_prob,
            alt_prob=alt_prob, ref_count=ref_count, alt_count=alt_count,
            other_count=other_count, hash=hash((chrom, position, ref, alt))
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


class BamCounts(object):

    def __init__(self, path):
        # Store paths and open BAM files
        self.path = path
        self.bam = pysam.AlignmentFile(self.path)
        # Get total counts
        self.total = bam.mapped
        # Get chromosome lengths
        self.chrom_lengths = {}
        for chrom in self.bams[0].references:
            chrom_len = self.bams[0].get_reference_length(chrom)
            self.chrom_lengths[chrom] = chrom_len

    def get_counts(self, chrom, start, end):
        counts = bam.count(chrom, start, end)
        return(counts)

    def close(self):
        for bam in self.bams:
            bam.close()


class VariantCounts(object):

    def __init__(self, paths):
        # Set paths and open files
        self.paths = paths
        self.counts = [pysam.TabixFile(p) for p in self.paths]

    def generator(self, chrom=None, start=None, end=None):
        # Generate iterators and loop through them
        iterators = [
            counts.fetch(
                reference=chrom, start=start, end=end,
                multiple_iterators=True
            ) for
            counts in self.counts
        ]
        while True:
            # Get next entry in count files
            try:
                lines = [iterator.next() for iterator in iterators]
            except StopIteration:
                break
            # Get variants
            sample_variants = [
                IndividualVariant.from_line(line) for line in lines
            ]
            # Check consistency of entries
            for sv in sample_variants[1:]:
                assert(sv == sample_variants[0])
            # Skip variants partially overlapping target regions
            variant_start = sample_variants[0].position - 1
            variant_end = variant_start + len(sample_variants[0].ref)
            if start is not None:
                if variant_start < start:
                    continue
            if end is not None:
                if variant_end > end:
                    continue
            # Return variants
            yield(sample_variants)

    def close(self):
        for counts in self.counts:
            counts.close()