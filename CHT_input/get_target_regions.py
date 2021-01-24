import collections
import pysam


class CountGenerator(object):

    def __init__(self, paths):
        # Set variables
        self.paths = paths
        # Open files
        self.counts = [
            pysam.TabixFile(p, parser=pysam.asTuple()) for p in self.paths
        ]
        # Create named tuple
        self.count_tup = collections.namedtuple(
            'count_tup', [
                'chrom', 'start', 'end', 'window_start', 'window_end', 'ref',
                'alt', 'sample', 'genotype', 'ref_count', 'alt_count',
                'other_count', 'ref_count_window', 'alt_count_window',
                'other_count_window', 'total_count'
            ]
        )

    def count_generator(self, chrom=None, start=None, end=None):
        # Generate iterators and loop through them
        iterators = [
            counts.fetch(reference=chrom, start=start, end=end) for
            counts in self.counts
        ]
        while True:
            # Get next entry in count files
            try:
                entries = [list(iterator.next()) for iterator in iterators]
            except StopIteration:
                break
            # Check consistency of entries
            variant = entries[0][0:7]
            for entry in entries[1:]:
                assert(variant == entry[0:7])
            # Format data
            count_tup = self.count_tup(
                chrom=variant[0],
                start=int(variant[1]),
                end=int(variant[2]),
                window_start=int(variant[3]),
                window_end=int(variant[4]),
                ref=variant[5],
                alt=variant[6],
                sample=[x[7] for x in entries],
                genotype=[x[8] for x in entries],
                ref_count=[int(x[9]) for x in entries],
                alt_count=[int(x[10]) for x in entries],
                other_count=[int(x[11]) for x in entries],
                ref_count_window=[int(x[12]) for x in entries],
                alt_count_window=[int(x[13]) for x in entries],
                other_count_window=[int(x[14]) for x in entries],
                total_count=[int(x[15]) for x in entries]
            )
            yield(count_tup)

    def close(self):
        for counts in self.counts:
            counts.close()


class FilterVariants(object):

    def __init__(self, inpaths, outpath):
        # Store paths and create generator
        self.paths = paths
        self.generator = CountGenerator(self.paths)

    def filter(self, min_het, min_minor_count, min_total_count):
        for variant in self.generator.count_generator():
            # Skip variants with limited heterozygous genotypes
            het_count = (
                variant.genotype.count('0/1') +
                variant.genotype.count('1/0')
            )
            if het_count < min_het:
                continue
            # Skip variants with low minor allele count
            ref_count = sum(variant.ref_count)
            alt_count = sum(variant.alt_count)
            minor_count = min(ref_count, alt_count)
            if minor_count < min_minor_count:
                continue
            # Skip variant with limited total counts in window
            total_count = (
                sum(variant.ref_count_window) +
                sum(variant.alt_count_window) +
                sum(variant.other_count_window)
            )
            if total_count < min_total_count:
                continue 

 

    
#     def add_counts(self, chrom, variants):
#         # Create iterators and loop through them
#         iterators = [counts.fetch(chrom) for counts in self.counts]
#         while True:
#             # Get next entry in count files
#             try:
#                 entries = [list(iterator.next()) for iterator in iterators]
#             except StopIteration:
#                 break
#             # Get variant and check consistency
#             variant = entries[0][0:4]
#             for entry in entries[1:]:
#                 assert(variant == entry[0:4])
#             # Format variant
#             variant[1] = int(variant[1])
#             variant = tuple(variant)
#             # Get counts
#             variant_counts = [map(int, entry[4:7]) for entry in entries]
#             ref_count, alt_count, other_count = zip(*variant_counts)
#             # Store counts
#             variants[variant]['ref_count'] = ref_count
#             variants[variant]['alt_count'] = alt_count
#             variants[variant]['other_count'] = other_count
#         return(variants)

#     def close(self):
#         for count in self.counts:
#             count.close()

# samples = [
#     'DGRP-28', 'DGRP-307', 'DGRP-399', 'DGRP-57', 'DGRP-639',
#     'DGRP-712', 'DGRP-714', 'DGRP-852', 'vgn'
# ]
# genotype = variant_genotype(
#     '/Users/rabinowi/wasp/test_data/variants/F1_haplotype_joint_call_GATK_stringent.vcf.gz',
#     samples=samples
# )
# counts = variant_counts([
#     '/Users/rabinowi/wasp/test_data/alignments/single_end.variant_counts.txt.gz',
#     '/Users/rabinowi/wasp/test_data/alignments/single_end.variant_counts.txt.gz'
# ])
# chrom = 'chr2L'
# chrom_variants = genotype.get_variant_genotypes('chr2L')
# chrom_variants = counts.add_counts('chr2L', chrom_variants)
# for variant in chrom_variants:
#     if len(chrom_variants[variant]) != 4:
#         raise ValueError('no counts for variant {} {} {} {}'.format(variant))
# count = 0
# for k, v in chrom_variants.items():
#     print(k, v)
#     count += 1
#     if count == 10:
#         break
