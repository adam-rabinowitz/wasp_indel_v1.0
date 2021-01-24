import argparse
import collections
import pysam


class vcf_genotype(object):

    def __init__(self, path, sample):
        # Save input parameters
        self.path = path
        self.sample = sample
        self.vcf = pysam.VariantFile(self.path)
        self.vcf.subset_samples([self.sample])

    def genotypes(self, chrom, starts, ends):
        # Create output list
        genotypes = collections.OrderedDict()
        # Loop through vcf entries in the specified regions
        for start, end in zip(starts, ends):
            for entry in self.vcf.fetch(chrom, start, end):
                # Extract variant data
                variant_id = (
                    entry.chrom, entry.pos, entry.ref, ','.join(entry.alts)
                )
                genotypes[variant_id] = {
                    k: v for k, v in entry.samples[self.sample].items()
                }
                # Convert phred scores to probability
                genotypes[variant_id]['PL'] = [
                    10 ** ((-p) / 10) for p in genotypes[variant_id]['PL']
                ]
                # Add id
                if entry.id is None:
                    genotypes[variant_id]['id'] = '_'.join(
                        map(str, variant_id)
                    )
                else:
                    genotypes[variant_id]['id'] = entry.id
        return(genotypes)

    def close(self):
        self.vcf.close()


class bam_coverage(object):

    def __init__(
        self, path, paired_end
    ):
        # Store parameters and open BAM file
        self.path = path
        self.paired_end = paired_end
        self.bam = pysam.AlignmentFile(self.path)
        # Check pairing
        for read in self.bam.head(100):
            if self.paired_end:
                assert(read.is_paired)
            else:
                assert(not read.is_paired)
        # Get total counts
        if self.paired_end:
            self.total = self.bam.mapped // 2
        else:
            self.total = self.bam.mapped

    def coverage(
        self, chrom, starts, ends
    ):
        # Set processing variables
        reads = set()
        # Get unique read names in region
        for start, end in zip(starts, ends):
            for read in self.bam.fetch(chrom, start, end):
                reads.add(read.query_name)
        # Count reads and return
        count = len(reads)
        return(count)

    def close(
        self
    ):
        self.bam.close()


class allele_counts(object):

    def __init__(
        self, path
    ):
        self.path = path
        self.tabix = pysam.TabixFile(path, parser=pysam.asTuple())

    def counts(
        self, chrom, starts, ends
    ):
        # Create output list
        allele_counts = collections.OrderedDict()
        # Loop through variants in regions
        for start, end in zip(starts, ends):
            for entry in self.tabix.fetch(chrom, start, end):
                # Extract variant data and store
                variant_id = tuple([
                    entry[0], int(entry[1]), entry[2], entry[3]
                ])
                variant_counts = {
                    'ref': int(entry[4]),
                    'alt': int(entry[5]),
                    'other': int(entry[6])
                }
                allele_counts[variant_id] = variant_counts
        return(allele_counts)

    def close(
        self
    ):
        self.tabix.close()


class cht_input(object):

    def __init__(
        self, bam, counts, vcf, sample, paired_end
    ):
        # Create file classes
        self.sample = sample
        self.coverage = bam_coverage(bam, paired_end)
        self.counts = allele_counts(counts)
        self.genotype = vcf_genotype(vcf, self.sample)
        self.regions = {}

    def calc_geno_prob(
        self, probs, n_alt
    ):
        ''' Function designed for multiple alternative alleles'''
        # Check number or provided probalites
        n_prob = sum(range(1, n_alt + 2))
        assert(n_prob == len(probs))
        # Get indices of heterozygous and homozygous alt alleles
        het_indices = [int(i * ((i + 1) / 2)) for i in range(1, 1 + n_alt)]
        alt_indices = [i for i in range(1, n_prob) if i not in het_indices]
        # Get probabilities for the three genotypes and return
        ref_prob = probs[0]
        het_prob = sum([probs[i] for i in het_indices])
        alt_prob = sum([probs[i] for i in alt_indices])
        # Generate probability tuple, check and return
        prob_tup = (ref_prob, het_prob, alt_prob)
        assert(abs(sum(prob_tup) - 1) < 0.01)
        return(prob_tup)

    def region_data(
        self, chrom, starts, ends
    ):
        # Adjust start to account for 0 based index
        starts = [s - 1 for s in starts]
        # Get genotypes of variants
        genotypes = self.genotype.genotypes(chrom, starts, ends)
        counts = self.counts.counts(chrom, starts, ends)
        coverage = self.coverage.coverage(chrom, starts, ends)
        # Check that there are counts for all variants
        assert(set(genotypes.keys()).issubset(set(counts.keys())))
        variants = list(genotypes.keys())
        # Get variant positions
        positions = ';'.join(
            [str(v[1]) for v in variants]
        )
        # Get heterozygous possiblities
        het_prob = []
        for variant in variants:
            n_alt = len(variant[3].split(','))
            het_prob.append(
                self.calc_geno_prob(genotypes[variant]['PL'], n_alt)[1]
            )
        het_prob = ';'.join(
            ['{:.2f}'.format(p) for p in het_prob]
        )
        # Set linkage probability
        linkage = ';'.join(
            ['1.00'] * len(genotypes)
        )
        # Get ref, alt and other counts
        ref_counts = ';'.join(
            [str(counts[v]['ref']) for v in variants]
        )
        alt_counts = ';'.join(
            [str(counts[v]['alt']) for v in variants]
        )
        other_counts = ';'.join(
            [str(counts[v]['other']) for v in variants]
        )
        # Create output and return
        output = [
            ';'.join(map(str, starts)), ';'.join(map(str, starts)),
            positions, het_prob, linkage, ref_counts, alt_counts,
            other_counts, str(coverage), str(self.coverage.total)
        ]
        return(output)

    def var_data(
        self, chrom, position
    ):
        # Get variant data
        genotypes = self.genotype.genotypes(
            chrom, (position - 1,), (position,)
        )
        assert(len(genotypes)) == 1
        variant = list(genotypes.keys())[0]
        variant_data = genotypes[variant]
        # Get alleles
        ref, alt = variant[2:4]
        n_alt = len(alt.split(','))
        # Calculate genotype
        genotype_prob = self.calc_geno_prob(variant_data['PL'], n_alt)
        genotype = '{:.2f}'.format(
            genotype_prob[1] + (genotype_prob[2] * 2)
        )
        # Get haplotype
        haplotype = '/'.join(
            [str(min(g, 1)) for g in variant_data['GT']]
        )
        # Create output and return
        output = [
           chrom, str(position), variant_data['id'], ref, alt,
           genotype, haplotype
        ]
        return(output)

    def generate_line(
        self, chrom, position, region_starts, region_ends
    ):
        # Generate variant data
        var_data = self.var_data(chrom, position)
        # Geberate region data
        region_id = (chrom, region_starts, region_ends)
        if region_id in self.regions:
            region_data = self.regions[region_id]
        else:
            region_data = self.region_data(
                chrom, region_starts, region_ends
            )
            self.regions[region_id] = region_data
        # Format data and return
        line = '\t'.join(var_data + region_data) + '\n'
        return(line)

    def close(
        self
    ):
        # Cloe all input file
        self.coverage.close()
        self.counts.close()
        self.genotype.close()


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="This program makes the input file for the combined "
        "haplotype test. The input BAM and counts file should be generated "
        "using the mapping scripts. The VCF file should be identical to, or "
        "a subset of, the VCF used to generate the count file. The region "
        "file sould contain the following 4 tabe delimite columns: "
        "chromosome, variant start, region start and region end. If the "
        "region is non-contiguous then the multiple starts and ends should "
        "be seperated by a ';' symbol. The genomic coordinates should be "
        "1-based so that region 1-3 is the first three bases of a chromsome."
    )
    parser.add_argument(
        "--bam", required=True, help=(
            "Coordinate sorted and indexed BAM file"
        )
    )
    parser.add_argument(
        "--counts", required=True, help=(
            "Coordinate sorted and tabix indexed count file"
        )
    )
    parser.add_argument(
        "--vcf", required=True, help=(
            "Coordinate sorted and tabix indexed VCF file"
        )
    )
    parser.add_argument(
        "--sample", required=True, help=(
            "Name of sample to be parsed from within VCF file"
        )
    )
    parser.add_argument(
        "--regions", required=True, help=(
            "Text file listing variants and their associated regions"
        )
    )
    parser.add_argument(
        "--cht", required=True, help=(
            "Output combined haplotype test file"
        )
    )
    parser.add_argument(
        "--paired_end", action='store_true', default=False, help=(
            "Indicates that reads are paired-end (default is single)."
        )
    )
    args = parser.parse_args()
    # Create cht_input object to generate cht input
    cht_generator = cht_input(
        bam=args.bam, counts=args.counts, vcf=args.vcf, sample=args.sample,
        paired_end=args.paired_end
    )
    # Open input and output files and loop through input file
    with open(args.regions) as infile, open(args.cht, 'wt') as outfile:
        for inline in infile:
            # Process line
            chrom, position, region_starts, region_ends = (
                inline.strip().split('\t')
            )
            position = int(position)
            region_starts = tuple([int(x) for x in region_starts.split(';')])
            region_ends = tuple([int(x) for x in region_ends.split(';')])
            # Create output line and write to file
            outline = cht_generator.generate_line(
                chrom=chrom, position=position, region_starts=region_starts,
                region_ends=region_ends
            )
            outfile.write(outline)
    # Tidy up
    cht_generator.close()
