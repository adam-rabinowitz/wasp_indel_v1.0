import argparse
import counts
import gzip


class MultipleBamCounts(object):

    def __init__(self, paths):
        # Store paths and create objects
        self.paths = paths
        self.bam_files = [counts.BamCounts(p) for p in self.paths]
        # Get chromosome lengths and check consistency
        self.chrom_lengths = self.bam_files[0].chrom_lengths
        for bam in self.bam_files[1:]:
            assert(self.chrom_lengths == bam.chrom_lengths)

    def get_counts(self, chrom, start, end):
        counts = [bam.get_counts(chrom, start, end) for bam in self.bam_files]
        count_sum = sum(counts)
        return(count_sum)

    def close(self):
        for bam in self.bam_files:
            bam.close()


class MultipleVariantCounts(object):

    def __init__(self, paths):
        # Store paths and create objects:
        self.paths = paths
        self.count_files = [counts.VariantCounts(p) for p in self.paths]
        # Set heterozygous haplotypes
        self.heterozygous = ('0|1', '1|0')

    def get_all_variants(self):
        # Create generators
        variant_generators = [
            count_file.get_region_variants() for
            count_file in self.count_files
        ]
        # Loop through generators and get variants
        while True:
            # Get next variants
            try:
                variants = [
                    generator.__next__() for
                    generator in variant_generators
                ]
            except StopIteration:
                break
            # Check identity and yield
            for variant in variants[1:]:
                assert(variant == variants[0])
            yield(variants)

    def get_as_count(self, chrom, start, end):
        # Set output allele specific counts
        as_count = 0
        # Loop through region in all count files
        for count_file in self.count_files:
            for variant in count_file.get_region_variants(
                chrom=chrom, start=start, end=end
            ):
                # Add ref and alt counts if variant is heterozygous
                if variant.haplotype in self.heterozygous:
                    as_count += variant.ref_count
                    as_count += variant.alt_count
        # Return allele specific counts
        return(as_count)

    def close(self):
        for count_file in self.count_files:
            count_file.close()


class FilterVariants(object):

    def __init__(
        self, bam_paths, variant_paths
    ):
        # Check arguments
        assert(len(bam_paths) == len(variant_paths))
        # Store arguments
        self.bam_paths = bam_paths
        self.variant_paths = variant_paths
        # Create classes
        self.bam_counts = MultipleBamCounts(self.bam_paths)
        self.variant_counts = MultipleVariantCounts(self.variant_paths)
        # Create counter
        self.counter = {
            'low_het': 0, 'low_minor': 0, 'low_as_reads': 0,
            'low_total_reads': 0, 'passed': 0
        }

    def filter_windows(
        self, min_het, min_minor, min_as_reads, min_total_reads, window,
        outpath
    ):
        # Open outfile
        if outpath.endswith('.gz'):
            outfile = gzip.open(outpath, 'wt')
        else:
            outfile = open(outpath, 'wt')
        # Loop through variants
        for test_variants in self.variant_counts.get_all_variants():
            # Get count of heterozygous variants
            test_haplotypes = [tv.haplotype for tv in test_variants]
            het_count = (
                test_haplotypes.count('0|1') + test_haplotypes.count('1|0')
            )
            # Count and skip variants with low number of heterozygotes
            if het_count < min_het:
                self.counter['low_het'] += 1
                continue
            # Calculate minor allele counts
            ref_count = sum([tv.ref_count for tv in test_variants])
            alt_count = sum([tv.alt_count for tv in test_variants])
            minor_count = min(ref_count, alt_count)
            # Count and skip variants with low number of minor alleles
            if minor_count < min_minor:
                self.counter['low_minor'] += 1
                continue
            # Define region
            chrom = test_variants[0].chrom
            chrom_length = self.bam_counts.chrom_lengths[chrom]
            region_start = max(test_variants[0].start - window, 0)
            region_end = min(test_variants[0].end + window, chrom_length)
            # Get total read counts within window
            total_reads = self.bam_counts.get_counts(
                chrom=chrom, start=region_start, end=region_end
            )
            # Count and skip regions with low number of total reads
            if total_reads < min_total_reads:
                self.counter['low_total_reads'] += 1
                continue
            # Get allele specific counts within region
            as_reads = self.variant_counts.get_as_count(
                chrom=chrom, start=region_start, end=region_end
            )
            # Count and skip regions with low number of as reads
            if as_reads < min_as_reads:
                self.counter['low_as_reads'] += 1
                continue
            # Write data to file
            self.counter['passed'] += 1
            outline = (
                '{chrom}\t{pos}\t{ref}\t{alt}\t{start}\t{end}\n'
            ).format(
                chrom=chrom, pos=test_variants[0].start + 1,
                ref=test_variants[0].ref, alt=test_variants[0].alt,
                start=region_start + 1, end=region_end
            )
            outfile.write(outline)
        print(self.counter)

    def close(self):
        self.bam_counts.close()
        self.variant_counts.close()


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="This program generates target regions for subsequent "
        "analysis using the CHT test. The input count files should be "
        "generated using the 'get_counts.py' script and each count file "
        "should contain data for an identical set of variants."
    )
    parser.add_argument(
        "--variants", required=True, nargs='+', help=(
            "Input variant count files."
        )
    )
    parser.add_argument(
        "--bams", required=True, nargs='+', help=(
            "Input variant BAM files."
        )
    )
    parser.add_argument(
        "--outfile", required=True, help=(
            "Output region file."
        )
    )
    parser.add_argument(
        "--min_het", required=True, type=int, help=(
            "Minimum number of individuals with heterozygous genotypes for "
            "the test variant."
        )
    )
    parser.add_argument(
        "--min_minor", required=True, type=int, help=(
            "Minimum total number of minor alleles counts for the test "
            "variant across all individuals."
        )
    )
    parser.add_argument(
        "--min_as_reads", required=True, type=int, help=(
            "Minimum total number of allele-specific reads in window "
            "around test variant across all individuals."
        )
    )
    parser.add_argument(
        "--min_total_reads", required=True, type=int, help=(
            "Minimum total number of reads in window around test variant "
            "across all individuals."
        )
    )
    parser.add_argument(
        "--window", required=True, type=int, help=(
            "Distance upstream and downstream of the test variant in which "
            "to exract allele specific and total read counts."
        )
    )
    args = parser.parse_args()
    # Filter variants
    variant_filter = FilterVariants(
        bam_paths=args.bams, variant_paths=args.variants
    )
    variant_filter.filter_windows(
        min_het=args.min_het, min_minor=args.min_minor,
        min_as_reads=args.min_as_reads, min_total_reads=args.min_total_reads,
        window=args.window, outpath=args.outfile
    )
    variant_filter.close()
