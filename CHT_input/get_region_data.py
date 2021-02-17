import argparse
import collections
import counts
import gzip
import math


class InputReader(object):

    def __init__(
        self, path
    ):
        # Store path and open file
        self.path = path
        if self.path.endswith('.gz'):
            self.infile = gzip.open(self.path, 'rt')
        else:
            self.infile = open(self.path, 'rt')
        # Create named tuple to store targets
        self.target = collections.namedtuple(
            'target', [
                'chrom', 'position', 'ref', 'alt', 'region_start', 'region_end'
            ]
        )

    def get_targets(
        self
    ):
        # Loop through lines in infile and get data
        for line in self.infile:
            # Extract line data
            line_list = line.strip().split('\t')
            chrom = line_list[0]
            position = int(line_list[1])
            ref, alt = line_list[2:4]
            region_start = int(line_list[4]) - 1
            region_end = int(line_list[5])
            # Create tuple and return
            target = self.target(
                chrom=chrom, position=position, ref=ref, alt=alt,
                region_start=region_start, region_end=region_end
            )
            yield(target)

    def close(
        self
    ):
        self.infile.close()


class OutputWriter(object):

    def __init__(
        self, path, adjust_hetp
    ):
        # Store input arguments
        self.path = path
        self.adjust_hetp = adjust_hetp
        # Open outfile
        if self.path.endswith('.gz'):
            self.outfile = gzip.open(self.path, 'wt')
        else:
            self.outfile = open(self.path, 'wt')
        # Open files and write header to each
        self.header = (
            "CHROM TEST.SNP.POS TEST.SNP.ID TEST.SNP.REF.ALLELE "
            "TEST.SNP.ALT.ALLELE TEST.SNP.GENOTYPE TEST.SNP.HAPLOTYPE "
            "REGION.START REGION.END REGION.SNP.POS REGION.SNP.HET.PROB "
            "REGION.SNP.LINKAGE.PROB REGION.SNP.REF.HAP.COUNT "
            "REGION.SNP.ALT.HAP.COUNT REGION.SNP.OTHER.HAP.COUNT "
            "REGION.READ.COUNT GENOMEWIDE.READ.COUNT\n"
        )
        self.outfile.write(self.header)
        # Tuple containing expected haplotypes
        self.haplotypes = set(['0|0', '0|1', '1|0', '1|1'])
        self.heterozygotes = set(['0|1', '1|0'])

    def close(
        self
    ):
        self.outfile.close()

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
        # Adjust heterozygous probabilities or...
        if self.adjust_hetp:
            # Adjust probabilities using total allele counts
            if self.adjust_hetp == 'total_counts':
                hetp = [
                    self.calculate_posterior_hetp(
                        prior=v.het_prob, ref=v.ref_total_count,
                        alt=v.alt_total_count
                    ) for v in region_variants
                ]
            # Adjust probabilities using allele specific allele counts
            elif self.adjust_hetp == 'as_counts':
                hetp = [
                    self.calculate_posterior_hetp(
                        prior=v.het_prob, ref=v.ref_as_count,
                        alt=v.alt_as_count
                    ) for v in region_variants
                ]
        # Or get raw probabailities
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
                        ref_hap_counts.append(variant.ref_as_count)
                        alt_hap_counts.append(variant.alt_as_count)
                        other_hap_counts.append(variant.other_as_count)
                    # Get counts for discordant heterozygous haplotypes
                    else:
                        ref_hap_counts.append(variant.alt_as_count)
                        alt_hap_counts.append(variant.ref_as_count)
                        other_hap_counts.append(variant.other_as_count)
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
        self, test_variant, region_variants, region_start, region_end,
        region_count, total_count
    ):
        # Create blank output line
        line_list = [
            test_variant.chrom, str(test_variant.start + 1),
            test_variant.id, test_variant.ref, test_variant.alt,
            'NA', 'NA', str(region_start + 1), str(region_end),
            'NA', 'NA', 'NA', 'NA', 'NA', 'NA', str(region_count),
            str(total_count)
        ]
        # Add genotype information if test haplotype is known
        if test_variant.haplotype:
            test_genotype = (
                test_variant.het_prob +
                (2 * test_variant.alt_prob)
            )
            line_list[5] = '{:.2f}'.format(test_genotype)
            line_list[6] = test_variant.haplotype
        # Remove region variants with unknown haplotype
        region_variants = [
            rv for rv in region_variants if rv.haplotype
        ]
        # Add region variant data
        if test_variant.haplotype and region_variants:
            # Add positions
            positions = [rv.start + 1 for rv in region_variants]
            line_list[9] = ';'.join(map(str, positions))
            # Add het probabilities
            line_list[10] = self.get_hetp(region_variants)
            # Add linkage
            linkage = ['1.00'] * len(region_variants)
            line_list[11] = ';'.join(linkage)
            # Calculate hap counts
            line_list[12:15] = self.get_haplotype_counts(
                test_variant=test_variant,
                region_variants=region_variants
            )
        # Join line and write to file
        outline = ' '.join(line_list) + '\n'
        self.outfile.write(outline)


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="This program generates target regions for subsequent "
        "analysis using the CHT test. The input count files should be "
        "generated using the 'get_counts.py' script and each count file "
        "should contain data for an identical set of variants."
    )
    parser.add_argument(
        "--variants", required=True, help="Input variant count file."
    )
    parser.add_argument(
        "--bam", required=True, help="Input variant BAM file."
    )
    parser.add_argument(
        "--regions", required=True, help="Input region file."
    )
    parser.add_argument(
        "--outfile", required=True, help="Output file."
    )
    parser.add_argument(
        "--adjust_hetp", default=None, choices=['as_counts', 'total_counts'],
        help="Adjust heterozygous probabilities."
    )
    args = parser.parse_args()
    # Open input file and output files
    variant_file = counts.VariantCounts(args.variants)
    bam_file = counts.BamCounts(args.bam)
    region_file = InputReader(args.regions)
    out_file = OutputWriter(args.outfile, adjust_hetp=args.adjust_hetp)
    # Loop through targets in region file
    for target in region_file.get_targets():
        # Get target variant
        test_variant = variant_file.get_target_variant(
            chrom=target.chrom, position=target.position,
            ref=target.ref, alt=target.alt
        )
        # Get region variants
        region_variants = variant_file.get_region_variants(
            chrom=target.chrom, start=target.region_start,
            end=target.region_end
        )
        # Get counts
        region_count = bam_file.get_counts(
            chrom=target.chrom, start=target.region_start,
            end=target.region_end
        )
        total_count = bam_file.total
        # Write output file
        out_file.write(
            test_variant=test_variant, region_variants=region_variants,
            region_start=target.region_start, region_end=target.region_end,
            region_count=region_count, total_count=total_count
        )
    # Close files
    variant_file.close()
    bam_file.close()
    region_file.close()
    out_file.close()
