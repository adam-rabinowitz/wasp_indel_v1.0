import argparse
import counts
import gzip
import itertools


class FilterCounts(object):

    def __init__(
        self
    ):
        # Set initial count values
        self.variants = 0
        self.low_het = 0
        self.low_minor = 0
        self.variants_passed = 0
        self.regions = 0
        self.low_as_reads = 0
        self.low_total_reads = 0
        self.regions_passed = 0

    def variant_counts(
        self
    ):
        variant_counts = (
            'Variants:\n'
            '  considered: {variants}\n'
            '  low heterozygous alleles: {low_het}\n'
            '  low minor allele count: {low_minor}\n'
            '  passed: {variants_passed}'
        ).format(
            variants=self.variants, low_het=self.low_het,
            low_minor=self.low_minor, variants_passed=self.variants_passed
        )
        return(variant_counts)

    def region_counts(
        self
    ):
        region_counts = (
            'Regions:\n'
            '  considered: {regions}\n'
            '  low allele specific reads: {low_as_reads}\n'
            '  low total reads: {low_total_reads}\n'
            '  passed: {regions_passed}'
        ).format(
            regions=self.regions, low_as_reads=self.low_as_reads,
            low_total_reads=self.low_total_reads,
            regions_passed=self.regions_passed
        )
        return(region_counts)


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
        counts = [
            bam.get_counts(chrom, [start], [end]) for bam in self.bam_files
        ]
        count_sum = sum(counts)
        return(count_sum)

    def close(self):
        for bam in self.bam_files:
            bam.close()


class MultipleVariantCounts(object):

    def __init__(
        self, paths
    ):
        # Store paths and create objects:
        self.paths = paths
        self.count_trees = [counts.CountTree(p, None) for p in self.paths]
        # Extract chromosomes
        chromosomes = [ct.chromosomes for ct in self.count_trees]
        for chrom in chromosomes[1:]:
            if chrom != chromosomes[0]:
                raise ValueError('inconsistent chromosomes')
        self.chromosomes = chromosomes[0]
        # Set heterozygous haplotypes
        self.heterozygous = ('0|1', '1|0')

    def read_chrom(
        self, chrom
    ):
        # Read variants for chromosome
        [ct.read_counts(chrom) for ct in self.count_trees]

    def variant_generator(
        self
    ):
        # Loop through chromosome
        variant_iter = itertools.zip_longest(
            *[ct.variants for ct in self.count_trees]
        )
        for variants in variant_iter:
            # Check variants
            for v in variants[1:]:
                if v != variants[0]:
                    raise StopIteration('differing variants')
            # Yield variants
            yield variants

    def get_variants(
        self, start=None, end=None
    ):
        # Create generators
        variant_iter = itertools.zip_longest(
            *[ct.get_variants([start], [end]) for ct in self.count_trees]
        )
        # Loop through generators and get variants
        for variants in variant_iter:
            # Check variants
            for v in variants[1:]:
                if v != variants[0]:
                    raise StopIteration('differing variants')
            # Yield variants
            yield variants

    def get_as_count(self, start, end):
        # Set output allele specific counts
        as_count = 0
        # Loop through all region variants in all count files
        for variants in self.get_variants(start=start, end=end):
            for variant in variants:
                # Add ref and alt counts if variant is heterozygous
                if variant.haplotype in self.heterozygous:
                    as_count += variant.ref_as_count
                    as_count += variant.alt_as_count
        # Return allele specific counts
        return(as_count)


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

    def close(
        self
    ):
        self.bam_counts.close()

    def region_generator(
        self, path
    ):
        # Open region file
        if path.endswith('.gz'):
            infile = gzip.open(path, 'rt')
        else:
            infile = open(path, 'rt')
        # Loop through input file skipping empty lines
        for line in infile:
            line = line.strip()
            if not line:
                continue
            # Extract data from line
            line_list = line.split('\t')
            chrom = line_list[0]
            test_start, test_end, region_start, region_end = map(
                int, line_list[1:]
            )
            # Yield data
            yield(chrom, test_start, test_end, region_start, region_end)
        # Close file
        infile.close()

    def count_hets(
        self, variants
    ):
        haplotypes = [v.haplotype for v in variants]
        het_count = haplotypes.count('0|1') + haplotypes.count('1|0')
        return(het_count)

    def minor_as_count(
        self, variants
    ):
        ref_as_count = sum([v.ref_as_count for v in variants])
        alt_as_count = sum([v.alt_as_count for v in variants])
        minor_as_count = min(ref_as_count, alt_as_count)
        return(minor_as_count)

    def filter_windows(
        self, min_het, min_minor, min_as_reads, min_total_reads, window,
        outpath
    ):
        # Open outfile
        if outpath.endswith('.gz'):
            outfile = gzip.open(outpath, 'wt')
        else:
            outfile = open(outpath, 'wt')
        # Create counter
        filter_counts = FilterCounts()
        # Loop through chromosome
        for chrom in self.variant_counts.chromosomes:
            self.variant_counts.read_chrom(chrom)
            # Loop through variants
            for test_variants in self.variant_counts.variant_generator():
                filter_counts.variants += 1
                # Skip test variants if heterozygous count is too low
                if self.count_hets(test_variants) < min_het:
                    filter_counts.low_het += 1
                    continue
                # Skip test variants if minor allele count is too low
                if self.minor_as_count(test_variants) < min_minor:
                    filter_counts.low_minor += 1
                    continue
                # Count passed variants and considered regions
                filter_counts.variants_passed += 1
                filter_counts.regions += 1
                # Define region
                chrom_length = self.bam_counts.chrom_lengths[chrom]
                region_start = max(test_variants[0].start - window, 0)
                region_end = min(test_variants[0].end + window, chrom_length)
                # Skip region if total reads is too low
                total_reads = self.bam_counts.get_counts(
                    chrom=chrom, start=region_start, end=region_end
                )
                if total_reads < min_total_reads:
                    filter_counts.low_total_reads += 1
                    continue
                # Skip region if allele specific reads are too low
                as_reads = self.variant_counts.get_as_count(
                    start=region_start, end=region_end
                )
                if as_reads < min_as_reads:
                    filter_counts.low_as_reads += 1
                    continue
                # Write acceptable test variants and regions to file
                filter_counts.regions_passed += 1
                outline = (
                    '{chrom}\t{pos}\t{ref}\t{alt}\t{start}\t{end}\n'
                ).format(
                    chrom=chrom, pos=test_variants[0].start + 1,
                    ref=test_variants[0].ref, alt=test_variants[0].alt,
                    start=region_start + 1, end=region_end
                )
                outfile.write(outline)
        # Print filter metrics
        print(filter_counts.variant_counts())
        print(filter_counts.region_counts())

    def filter_regions(
        self, min_het, min_minor, min_as_reads, min_total_reads,
        regions, outpath
    ):
        # Open outfile
        if outpath.endswith('.gz'):
            outfile = gzip.open(outpath, 'wt')
        else:
            outfile = open(outpath, 'wt')
        # Create counter
        filter_counts = FilterCounts()
        # Set current chromosome
        current_chrom = None
        previous_chroms = []
        # Loop thorugh region file
        reg_gen = self.region_generator(regions)
        for chrom, test_start, test_end, region_start, region_end in reg_gen:
            filter_counts.regions += 1
            # Check chromosome
            if chrom != current_chrom:
                # Raise error if chromosome has been seen before
                if chrom in previous_chroms:
                    raise ValueError('region file must chromosome sorted')
                # Store new chromosome and read variants
                previous_chroms.append(chrom)
                current_chrom = chrom
                self.variant_counts.read_chrom(chrom)
            # Adjust stars to account for zero indexing
            test_start -= 1
            region_start -= 1
            # Skip region if total reads is too low
            total_reads = self.bam_counts.get_counts(
                chrom=chrom, start=region_start, end=region_end
            )
            if total_reads < min_total_reads:
                filter_counts.low_total_reads += 1
                continue
            # Skip region if allele specific reads are too low
            as_reads = self.variant_counts.get_as_count(
                start=region_start, end=region_end
            )
            if as_reads < min_as_reads:
                filter_counts.low_as_reads += 1
                continue
            # Loop through variants
            filter_counts.regions_passed += 1
            for test_variants in self.variant_counts.get_variants(
                start=test_start, end=test_end
            ):
                filter_counts.variants += 1
                # Skip test variants if heterozygous count is too low
                if self.count_hets(test_variants) < min_het:
                    filter_counts.low_het += 1
                    continue
                # Skip test variants if minor allele count is too low
                if self.minor_as_count(test_variants) < min_minor:
                    filter_counts.low_minor += 1
                    continue
                # Write acceptable test variants and regions to file
                filter_counts.variants_passed += 1
                outline = (
                    '{chrom}\t{pos}\t{ref}\t{alt}\t{start}\t{end}\n'
                ).format(
                    chrom=chrom, pos=test_variants[0].start + 1,
                    ref=test_variants[0].ref, alt=test_variants[0].alt,
                    start=region_start + 1, end=region_end
                )
                outfile.write(outline)
        # Print filter metrics
        print(filter_counts.region_counts())
        print(filter_counts.variant_counts())


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="This program generates target regions for subsequent "
        "analysis using the CHT test. The input variant count files should "
        "be generated using the 'get_counts.py' script and each file "
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
            "Minimum total number of allele-specific reads in region "
            "associated to test variant across all individuals."
        )
    )
    parser.add_argument(
        "--min_total_reads", required=True, type=int, help=(
            "Minimum total number of reads in region associated to test "
            "variant across all individuals."
        )
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--window", type=int, help=(
            "Distance upstream and downstream of the test variant in which "
            "to exract allele specific and total read counts."
        )
    )
    group.add_argument(
        "--regions", help=(
            "File containing region coordinates. File should contain five "
            "columns: chrom, start and end of test variant region and "
            "start and end of region in which to perform test."
        )
    )
    args = parser.parse_args()
    # Create variant filter object
    variant_filter = FilterVariants(
        bam_paths=args.bams, variant_paths=args.variants
    )
    # Get test variants and associate windows if window set else...
    if args.window:
        variant_filter.filter_windows(
            min_het=args.min_het, min_minor=args.min_minor,
            min_as_reads=args.min_as_reads,
            min_total_reads=args.min_total_reads,
            window=args.window, outpath=args.outfile
        )
    # or get test variants within regions
    elif args.regions:
        variant_filter.filter_regions(
            min_het=args.min_het, min_minor=args.min_minor,
            min_as_reads=args.min_as_reads,
            min_total_reads=args.min_total_reads,
            regions=args.regions, outpath=args.outfile
        )
    variant_filter.close()
