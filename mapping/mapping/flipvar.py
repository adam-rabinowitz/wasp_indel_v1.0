import collections
import itertools


class HaplotypeError(Exception):
    pass


class FlipVar(object):

    def __init__(self, reads, variants, samples=None, offset=33):
        # Check arguments
        assert(len(reads) == len(variants))
        # Check read pairing
        if len(reads) == 1:
            assert(not reads[0].is_paired)
            self.paired = False
        else:
            assert(len(reads) == 2)
            assert(reads[0].is_read1 and reads[1].is_read2)
            assert(reads[0].query_name == reads[1].query_name)
            assert(reads[0].reference_id == reads[1].reference_id)
            self.paired = True
        # Extract read data
        self.name = reads[0].query_name
        self.read_data = [
            (read.query_sequence, list(read.query_qualities), read.is_reverse)
            for read in reads
        ]
        self.position = [reads[0].reference_id]
        for read in reads:
            self.position.extend([read.reference_start, read.reference_end])
        # Extract variant data
        self.variants = variants
        self.common_variants = []
        if self.paired:
            for i1, variant in enumerate(variants[0]):
                try:
                    i2 = variants[1].index(variant)
                except ValueError:
                    continue
                self.common_variants.append((i1, i2))
        # Set variables for allele flipping
        self.haplotypes = None
        self.flipped_reads = None
        # Set variable for flipping alleles
        self.complement = {
            'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'
        }
        self.offset = offset

    def get_unique_variants(self):
        unique_variants = set(itertools.chain.from_iterable(self.variants))
        return(unique_variants)

    def max_variant_count(self):
        '''Get maximum number of variants overlapping a single read'''
        max_count = max([len(v) for v in self.variants])
        return(max_count)

    def variants_overlap(self):
        '''Function checks if any variants overlap'''
        for read_variants in self.variants:
            current_end = 0
            for variant in read_variants:
                if variant.start < current_end or '*' in variant.alleles:
                    return(True)
                current_end = max(current_end, variant.end)
        return(False)

    def conflicting_alleles(self):
        '''Check for conflicting alleles between read pairs'''
        for i1, i2 in self.common_variants:
            read1_variant = self.variants[0][i1]
            read2_variant = self.variants[1][i2]
            if read1_variant.read_allele != read2_variant.read_allele:
                return(True)
        return(False)

    def count_variant_alleles(self):
        ''' Function counts matches between read alleles and vcf alleles'''
        # Loop through varaints and extract read allele
        counts = collections.Counter()
        for variant in itertools.chain.from_iterable(self.variants):
            try:
                allele_index = variant.alleles.index(variant.read_allele)
            except ValueError:
                allele_index = None
            # Count match
            if allele_index is None:
                counts['other'] += 1
            elif allele_index == 0:
                counts['ref'] += 1
            else:
                counts['alt'] += 1
        return(counts)

    def __get_sample_haplotypes(self, samples):
        # Loop through read variants and get read haplotypes
        haplotypes = []
        for read_variants in self.variants:
            read_haplotypes = set()
            # Add initial haplotype as none tuple
            initial_read_haplotype = tuple(
                [None] * len(read_variants)
            )
            read_haplotypes.add(initial_read_haplotype)
            # Loop through samples and get phased haplotypes
            for sample in samples:
                # Extract alleles for each sample
                sample_alleles = [
                    range(len(v.genotypes[sample]))
                    for v in read_variants
                ]
                # Convert alleles to haplotypes
                for sample_haplotype in itertools.zip_longest(*sample_alleles):
                    if None in sample_haplotype:
                        raise HaplotypeError('missing alleles')
                    read_haplotypes.add(sample_haplotype)
            # Store haplotypes
            read_haplotypes = list(read_haplotypes)
            haplotypes.append(read_haplotypes)
        return(haplotypes)

    def __get_possible_haplotypes(self):
        # Loop through variants for each read
        haplotypes = []
        for read_variants in self.variants:
            read_alleles = []
            for variant in read_variants:
                # Find and store all possible variant alleles
                variant_alleles = [None]
                for index, allele in enumerate(variant.alleles):
                    if allele != variant.read_allele:
                        variant_alleles.append(index)
                read_alleles.append(variant_alleles)
            # Get all possible haplotypes and store
            read_haplotypes = list(itertools.product(*read_alleles))
            haplotypes.append(read_haplotypes)
        return(haplotypes)

    def __compatible_haplotypes(self, hap1, hap2):
        for i1, i2 in self.common_variants:
            if hap1[i1] != hap2[i2]:
                return(False)
        return(True)

    def get_haplotypes(self, samples):
        # Get haplotypes for each read
        if samples:
            haplotypes = self.__get_sample_haplotypes(samples)
        else:
            haplotypes = self.__get_possible_haplotypes()
        # Further process paired read haplotypes or...
        if self.paired and haplotypes[0] and haplotypes[1]:
            self.haplotypes = [[], []]
            for hap1, hap2 in itertools.product(*haplotypes):
                if self.__compatible_haplotypes(hap1, hap2):
                    self.haplotypes[0].append(hap1)
                    self.haplotypes[1].append(hap2)
        # ...store haplotypes for single reads
        else:
            self.haplotypes = haplotypes

    def count_haplotypes(self):
        n_hap = len(self.haplotypes[0])
        return(n_hap)

    def __generate_flipped_read(
        self, read_data, variants, haplotype
    ):
        """Generate set of reads with all possible combinations"""
        # Extract sequence data
        sequence, quality, reverse = read_data
        # Loop through variants and haplotype in reverse
        for variant, allele_index in zip(variants[::-1], haplotype[::-1]):
            # Keep data unchanged if allele index is None
            if allele_index is None:
                continue
            # Get new and old alleles and check
            new_allele = variant.alleles[allele_index]
            old_allele = sequence[variant.read_start:variant.read_end]
            assert(old_allele == variant.read_allele)
            # Keep data unchanged if old allele equals new allele
            if new_allele == old_allele or new_allele == '*':
                continue
            # Update sequence
            sequence = (
                sequence[:variant.read_start] +
                new_allele +
                sequence[variant.read_end:]
            )
            # Keep same quality if length is unchanged or...
            if len(old_allele) == len(new_allele):
                continue
            # ...modify quality
            else:
                # Calculate mean quality across variant
                qualities = quality[
                    variant.read_start:variant.read_end
                ]
                mean_quality = sum(qualities) // len(qualities)
                # Update quality
                quality = (
                    quality[:variant.read_start] +
                    [mean_quality] * len(new_allele) +
                    quality[variant.read_end:]
                )
        # Flip reads to account for strand
        if reverse:
            quality = ''.join([chr(q + self.offset) for q in quality[::-1]])
            sequence = ''.join([self.complement[s] for s in sequence[::-1]])
        else:
            quality = ''.join([chr(q + self.offset) for q in quality])
        # Return modified read as a flipped read
        new_read = (sequence, quality)
        return(new_read)

    def generate_flipped_reads(self):
        self.flipped_reads = []
        # Generate interleaved flipped sequences for paired end reads
        if self.paired:
            read1_data, read2_data = self.read_data
            read1_variants, read2_variants = self.variants
            for read1_haplotype, read2_haplotype in zip(*self.haplotypes):
                flipped = (
                    self.__generate_flipped_read(
                        read1_data, read1_variants, read1_haplotype
                    ),
                    self.__generate_flipped_read(
                        read2_data, read2_variants, read2_haplotype
                    )
                )
                self.flipped_reads.append(flipped)
        # Generate flipped sequence for signle end reads
        else:
            read1_data = self.read_data[0]
            read1_variants = self.variants[0]
            for read1_haplotype in self.haplotypes[0]:
                flipped = (
                    self.__generate_flipped_read(
                        read1_data, read1_variants, read1_haplotype
                    ),
                )
                self.flipped_reads.append(flipped)

    def min_length(self):
        '''Calculate min length of reads after flipping'''
        flipped_iter = itertools.chain.from_iterable(self.flipped_reads)
        min_length = min([len(f[0]) for f in flipped_iter])
        return(min_length)

    def generate_fastq(self):
        # set fastq string and count read number
        fastq = ''
        n_seq = len(self.flipped_reads)
        position = '-'.join(map(str, self.position))
        # Process new reads sequentially
        for i, reads in enumerate(self.flipped_reads):
            # Generate identifier
            identifier = '@{}.{}.{}.{:06d}'.format(
                self.name, position, n_seq, i
            )
            # Loop through paired reads
            for sequence, quality in reads:
                fastq_entry = '{}\n{}\n+\n{}\n'.format(
                    identifier, sequence, quality
                )
                fastq += fastq_entry
        # Return fastq
        return(fastq)
