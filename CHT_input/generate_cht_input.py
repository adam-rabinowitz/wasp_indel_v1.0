import argparse
import collections
import itertools
import math
import pysam


VariantTuple = collections.namedtuple(
    '_Variant', [
        'chrom', 'position', 'id', 'ref', 'alt', 'haplotype', 'ref_prob',
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

    def __init__(self, paths):
        # Store paths and open BAM files
        self.paths = paths
        self.bams = [pysam.AlignmentFile(p) for p in self.paths]
        # Get total counts
        self.total_counts = [bam.mapped for bam in self.bams]
        # Get chromosome lengths
        self.chrom_lengths = {}
        for chrom in self.bams[0].references:
            chrom_len = self.bams[0].get_reference_length(chrom)
            self.chrom_lengths[chrom] = chrom_len

    def get_counts(self, chrom, start, end):
        bam_counts = [
            bam.count(chrom, start, end) for bam in self.bams
        ]
        return(bam_counts)

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


class WriteOutput(object):

    def __init__(self, paths, adjust_hetp):
        # Store input arguments
        self.paths = paths
        self.adjust_hetp = adjust_hetp
        # Open files and write header to each
        self.header = (
            "CHROM TEST.SNP.POS TEST.SNP.ID TEST.SNP.REF.ALLELE "
            "TEST.SNP.ALT.ALLELE TEST.SNP.GENOTYPE TEST.SNP.HAPLOTYPE "
            "REGION.START REGION.END REGION.SNP.POS REGION.SNP.HET.PROB "
            "REGION.SNP.LINKAGE.PROB REGION.SNP.REF.HAP.COUNT "
            "REGION.SNP.ALT.HAP.COUNT REGION.SNP.OTHER.HAP.COUNT "
            "REGION.READ.COUNT GENOMEWIDE.READ.COUNT\n"
        )
        self.files = [open(p, 'wt') for p in self.paths]
        for f in self.files:
            f.write(self.header)
        # Tuple containing expected haplotypes
        self.haplotypes = set(['0|0', '0|1', '1|0', '1|1'])
        self.heterozygotes = set(['0|1', '1|0'])

    def addlogs(
        self, loga, logb
    ):
        sum_logs = max(loga, logb) + math.log(1 + math.exp(-abs(loga - logb)))
        return(sum_logs)

    def calculate_posterior_hetp(
        self, prior, ref, alt, error=0.01
    ):
        # Adjust prior
        prior = min(0.99, prior)
        # Calculate genotype likelihoods
        badlike = self.addlogs(
            (math.log(error) * ref) + (math.log(1 - error) * alt),
            (math.log(1 - error) * ref) + (math.log(error) * alt)
        )
        goodlike = (math.log(0.5) * ref) + (math.log(0.5) * alt)
        # avoid overflow (very close to 1.0)
        if goodlike - badlike > 40:
            hetp_post = 1.0
        # Or calculate posterior
        else:
            hetp_post = (
                prior * math.exp(goodlike - badlike) /
                (prior * math.exp(goodlike - badlike) + (1.0 - prior))
            )
        # Return posterior
        return(hetp_post)

    def get_hetp(
        self, region_variants
    ):
        # Get probabilities
        if self.adjust_hetp:
            # Adjust probabilites
            hetp = [
                self.calculate_posterior_hetp(
                    prior=v.het_prob, ref=v.ref_count, alt=v.alt_count
                ) for v in region_variants
            ]
        else:
            hetp = [v.het_prob for v in region_variants]
        # Convert output to string and return
        hetp_str = ';'.join('{:.2f}'.format(p) for p in hetp)
        return(hetp_str)

    def get_haplotype_counts(
        self, test_variant, region_variants
    ):
        # Check arguments
        assert(test_variant.haplotype in self.haplotypes)
        assert(len(region_variants) > 0)
        # Process heterozygotic test variants...
        if test_variant in self.heterozygotes:
            # Set zero counts
            ref_hap_counts = []
            alt_hap_counts = []
            other_hap_counts = []
            # Loop through region variants and check haplotype
            for variant in region_variants:
                assert(variant.haplotype in self.haplotypes)
                # Extract counts for heterozygotic variants...
                if variant.haplotype in self.heterozygotes:
                    if variant.haplotype == test_variant.haplotype:
                        ref_hap_counts.append(variant.ref_count)
                        alt_hap_counts.append(variant.alt_count)
                        other_hap_counts.append(variant.other_count)
                    # Get counts for discordant heterozygous haplotypes
                    else:
                        ref_hap_counts.append(variant.alt_count)
                        alt_hap_counts.append(variant.ref_count)
                        other_hap_counts.append(variant.other_count)
                # or set counts to zero
                else:
                    ref_hap_counts.append(0)
                    alt_hap_counts.append(0)
                    other_hap_counts.append(0)
        # or set counts to zero for homozygotes
        else:
            ref_hap_counts = [0] * len(region_variants)
            alt_hap_counts = [0] * len(region_variants)
            other_hap_counts = [0] * len(region_variants)
        # Convert counts to string and return
        ref_hap_counts = ';'.join(map(str, ref_hap_counts))
        alt_hap_counts = ';'.join(map(str, alt_hap_counts))
        other_hap_counts = ';'.join(map(str, other_hap_counts))
        return(ref_hap_counts, alt_hap_counts, other_hap_counts)

    def write(
        self, test_variants, region_variants, region_start, region_end,
        region_counts, total_counts
    ):
        # Create lines for each file
        for i, outfile in enumerate(self.files):
            # Get variants and metrics
            sample_test_variant = test_variants[i]
            sample_region_variants = [rv[i] for rv in region_variants]
            region_count = region_counts[i]
            total_count = total_counts[i]
            # Create blank output line
            line_list = [
                sample_test_variant.chrom, str(sample_test_variant.position),
                sample_test_variant.id, sample_test_variant.ref,
                sample_test_variant.alt, 'NA', 'NA', str(region_start + 1),
                str(region_end), 'NA', 'NA', 'NA', 'NA', 'NA', 'NA',
                str(region_count), str(total_count)
            ]
            # Add genotype information if test haplotype is known
            if sample_test_variant.haplotype:
                test_genotype = (
                    sample_test_variant.het_prob +
                    (2 * sample_test_variant.alt_prob)
                )
                line_list[5] = '{:.2f}'.format(test_genotype)
                line_list[6] = sample_test_variant.haplotype
            # Remove region variants with unknown haplotype
            filtered_region_variants = [
                srv for srv in sample_region_variants if srv.haplotype
            ]
            # Add region variant data
            if sample_test_variant.haplotype and filtered_region_variants:
                # Add positions
                positions = [
                    frv.position for frv in filtered_region_variants
                ]
                line_list[9] = ';'.join(map(str, positions))
                # Add het probabilities
                line_list[10] = self.get_hetp(filtered_region_variants)
                # Add linkage
                linkage = ['1.00'] * len(filtered_region_variants)
                line_list[11] = ';'.join(linkage)
                # Calculate hap counts
                line_list[12:15] = self.get_haplotype_counts(
                    test_variant=sample_test_variant,
                    region_variants=filtered_region_variants
                )
            # Join line and write to file
            outline = ' '.join(line_list) + '\n'
            outfile.write(outline)


class FilterVariants(object):

    def __init__(
        self, bam_paths, count_paths, out_paths, adjust_hetp
    ):
        # Check arguments
        assert(len(bam_paths) == len(count_paths))
        assert(len(bam_paths) == len(out_paths))
        # Store arguments
        self.bam_paths = bam_paths
        self.count_paths = count_paths
        # Create classes
        self.bam_counts = BamCounts(bam_paths)
        self.variant_counts = VariantCounts(count_paths)
        self.out_files = WriteOutput(
            out_paths, adjust_hetp
        )
        # Create counter
        self.counter = {
            'low_het': 0, 'low_minor': 0, 'low_as_reads': 0,
            'low_total_reads': 0, 'passed': 0
        }

    def filter_all(
        self, min_het, min_minor, min_as_reads, min_total_reads, window
    ):
        # Loop through variants
        for test_variants in self.variant_counts.generator():
            # Skip variants with limited heterozygous genotypes
            test_haplotypes = [tv.haplotype for tv in test_variants]
            het_count = (
                test_haplotypes.count('0|1') +
                test_haplotypes.count('1|0')
            )
            if het_count < min_het:
                self.counter['low_het'] += 1
                continue
            # Skip variants with low minor allele count
            ref_count = sum([tv.ref_count for tv in test_variants])
            alt_count = sum([tv.alt_count for tv in test_variants])
            minor_count = min(ref_count, alt_count)
            if minor_count < min_minor:
                self.counter['low_minor'] += 1
                continue
            # Define region
            chrom = test_variants[0].chrom
            chrom_length = self.bam_counts.chrom_lengths[chrom]
            test_variant_start = test_variants[0].position - 1
            test_variant_end = test_variant_start + len(test_variants[0].ref)
            region_start = max(test_variant_start - window, 0)
            region_end = min(test_variant_end + window, chrom_length)
            # Get counts for window
            region_counts = self.bam_counts.get_counts(
                chrom, region_start, region_end
            )
            if sum(region_counts) < min_total_reads:
                self.counter['low_total_reads'] += 1
                continue
            # Get variants within region
            region_variants = [
                rv for rv in self.variant_counts.generator(
                    chrom, region_start, region_end
                )
            ]
            # Check allele specific counts
            as_counts = 0
            for rv in itertools.chain.from_iterable(region_variants):
                if rv.haplotype in ('0|1', '1|0'):
                    as_counts += (rv.ref_count + rv.alt_count)
            if as_counts < min_as_reads:
                self.counter['low_as_reads'] += 1
                continue
            # Write data to file
            self.counter['passed'] += 1
            self.out_files.write(
                test_variants=test_variants, region_variants=region_variants,
                region_start=region_start, region_end=region_end,
                region_counts=region_counts,
                total_counts=self.bam_counts.total_counts
            )
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
        "--counts", required=True, nargs='+', help=(
            "Input variant count files."
        )
    )
    parser.add_argument(
        "--bams", required=True, nargs='+', help=(
            "Input variant BAM files."
        )
    )
    parser.add_argument(
        "--outfiles", required=True, nargs='+', help=(
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
    parser.add_argument(
        "--adjust_hetp", action='store_true', help=(
            "Distance upstream and downstream of the test variant in which "
            "to exract allele specific and total read counts (default=False)."
        )
    )
    args = parser.parse_args()
    # Filter variants
    variant_filter = FilterVariants(
        bam_paths=args.bams, count_paths=args.counts, out_paths=args.outfiles,
        adjust_hetp=args.adjust_hetp
    )
    variant_filter.filter_all(
        min_het=args.min_het, min_minor=args.min_minor,
        min_as_reads=args.min_as_reads, min_total_reads=args.min_total_reads,
        window=args.window
    )
    variant_filter.close()
