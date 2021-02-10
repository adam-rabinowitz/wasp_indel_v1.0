import collections
import pysam


VariantTuple = collections.namedtuple(
    '_Variant', [
        'chrom', 'start', 'end', 'id', 'ref', 'alt', 'haplotype', 'ref_prob',
        'het_prob', 'alt_prob', 'ref_count', 'alt_count', 'other_count',
        'hash'
    ]
)


class IndividualVariant(VariantTuple):

    @classmethod
    def from_line(cls, line):
        # Split line data and extract values
        line_data = line.strip().split('\t')
        assert(len(line_data) == 12)
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
            ref_prob, het_prob, alt_prob = None, None, None
        else:
            haplotype = line_data[5]
            ref_prob, het_prob, alt_prob = map(float, line_data[6:9])
        # Get counts
        ref_count, alt_count, other_count = map(int, line_data[9:12])
        # Generate hash
        new_variant = cls(
            chrom=chrom, start=start, end=end, id=id, ref=ref, alt=alt,
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


class VariantCounts(object):

    def __init__(self, path):
        # Set paths and open files
        self.path = path
        self.count_file = pysam.TabixFile(self.path)

    def get_target_variant(self, chrom, position, ref, alt):
        # Get positions of test variant
        variant_start = position - 1
        variant_end = variant_start + len(ref)
        # Get putative variants
        variants = [
            IndividualVariant.from_line(line) for
            line in self.count_file.fetch(
                reference=chrom, start=variant_start, end=variant_end,
                multiple_iterators=True
            )
        ]
        # Filter variants
        variants = [v for v in variants if v.ref == ref and v.alt == alt]
        assert(len(variants) == 1)
        return(variants[0])

    def get_region_variants(self, chrom=None, start=None, end=None):
        # Sequentially generate variants in region
        for line in self.count_file.fetch(
            reference=chrom, start=start, end=end,
            multiple_iterators=True
        ):
            variant = IndividualVariant.from_line(line)
            # Check variant is contained within supplied region
            if start and variant.start < start:
                continue
            if end and variant.end > end:
                continue
            # Yield variant
            yield(variant)

    def close(self):
        self.count_file.close()


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

    def get_counts(self, chrom, start, end):
        counts = self.bam_file.count(chrom, start, end)
        return(counts)

    def close(self):
        self.bam_file.close()
