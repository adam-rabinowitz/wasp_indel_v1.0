import argparse
import collections
import gzip
import itertools
import pysam
import vartree


class BamGenerator():

    def __init__(
        self, bam, min_mapq
    ):
        # Store parameters and open BAM file
        self.bam_path = bam
        self.min_mapq = min_mapq
        self.bam = pysam.AlignmentFile(self.bam_path)
        # Check pairing
        paired = [read.is_paired for read in self.bam.head(100)]
        if any(paired):
            assert(all(paired))
            self.paired = True
        else:
            self.paired = False
        # Get bam metrics
        self.chromosomes = self.bam.references
        self.total = self.bam.mapped + self.bam.unmapped
        self.nocoordinate = self.bam.nocoordinate
        # Create counter
        self.counter = {
            'total': 0,
            'secondary': 0,
            'supplementary': 0,
            'unmapped': 0,
            'mate_unmapped': 0,
            'different_chromosomes': 0,
            'improper_pair': 0,
            'low_mapq': 0,
            'passed': 0
        }

    def close(
        self
    ):
        self.bam.close()

    def get_reads(
        self, chromosome
    ):
        # Create cache to store unfound reads
        read_pair_cache = {}
        # Find reads properly mapped to the chromosome
        for read in self.bam.fetch(contig=chromosome):
            self.counter['total'] += 1
            mapped_reads = None
            # Count and skip secondary reads
            if read.is_secondary:
                self.counter['secondary'] += 1
            # Count and skip supplementary reads
            elif read.is_supplementary:
                self.counter['supplementary'] += 1
            # Count and skip unmapped reads
            elif read.is_unmapped:
                self.counter['unmapped'] += 1
            # Store single reads for further filtering...
            elif not self.paired:
                assert(not read.is_paired)
                mapped_reads = [read]
            # or process paired to get read pairs
            else:
                # Count mate unmapped reads
                if read.mate_is_unmapped:
                    self.counter['mate_unmapped'] += 1
                # Count reads whose mate is mapped to different chromosome
                elif read.reference_id != read.next_reference_id:
                    self.counter['different_chromosomes'] += 1
                # Count reads in improper pairs
                elif not read.is_proper_pair:
                    self.counter['improper_pair'] += 1
                # Get paired reads if mate has been cached...
                elif read.query_name in read_pair_cache:
                    cached_read = read_pair_cache.pop(read.query_name)
                    if read.is_read1:
                        assert(cached_read.is_read2)
                        mapped_reads = [read, cached_read]
                    else:
                        assert(cached_read.is_read1)
                        assert(read.is_read2)
                        mapped_reads = [cached_read, read]
                # or store first in pair
                else:
                    read_pair_cache[read.query_name] = read
            #  Further filter properly mapped reads
            if mapped_reads is not None:
                # Count and skip reads with low mapping quality
                if any(
                    [r.mapping_quality < self.min_mapq for r in mapped_reads]
                ):
                    self.counter['low_mapq'] += len(mapped_reads)
                else:
                    self.counter['passed'] += len(mapped_reads)
                    yield(mapped_reads)
        # Check read pair cache is empty
        assert(len(read_pair_cache) == 0)


class FastqWriter(object):

    def __init__(
        self, path, offset=33
    ):
        # Store variables
        self.path = path
        self.offset = offset
        # Open output file
        self.fastq = gzip.open(self.path, 'wt')

    def close(
        self
    ):
        self.fastq.close()

    def rev_comp(
        self, sequence
    ):
        # Reverse sequence and make complement
        rev_sequence = sequence[::-1]
        translation = rev_sequence.maketrans(
            'ATCGMRWSYKNatcgmrwsykn',
            'TAGCKYWSRMNtagckywsrmn'
        )
        revcomp_sequence = rev_sequence.translate(translation)
        return(revcomp_sequence)

    def quality_string(
        self, quality, reverse
    ):
        # Convert to characters
        if reverse:
            characters = [chr(q + self.offset) for q in quality[::-1]]
        else:
            characters = [chr(q + self.offset) for q in quality]
        # Join and return
        quality_string = ''.join(characters)
        return(quality_string)

    def write_fastq(
        self, read, new_reads
    ):
        # Check reads
        assert(read.query_sequence == new_reads[0].sequence)
        # Get read name
        read_name = read.query_name
        reference_id = read.reference_id
        # Get position of original mapped reads and number of new reads
        position = "{id}-{start}-{end}".format(
            id=reference_id, start=read.reference_start,
            end=read.reference_end
        )
        n_seq = len(new_reads)
        # Process new reads sequentially
        for i, new_read in enumerate(new_reads):
            # Get new name
            fastq_name = "{}.{}.{}.{:06d}".format(
                read_name, position, n_seq, i
            )
            # Get new read sequence
            if read.is_reverse:
                fastq_seq = self.rev_comp(new_read.sequence)
                fastq_qual = self.quality_string(
                    new_read.quality, reverse=True
                )
            else:
                fastq_seq = new_read.sequence
                fastq_qual = self.quality_string(
                    new_read.quality, reverse=False
                )
            assert(len(fastq_seq) == len(fastq_qual))
            # Write read to file
            fastq_entry = '@{name}\n{seq}\n+\n{quality}\n'.format(
                name=fastq_name, seq=fastq_seq, quality=fastq_qual
            )
            self.fastq.write(fastq_entry)

    def write_pair_fastq(
        self, read1, read2, new_read_pairs
    ):
        # Check reads
        assert(read1.is_read1 & read2.is_read2)
        assert(read1.query_name == read2.query_name)
        assert(read1.reference_id == read2.reference_id)
        assert(read1.query_sequence == new_read_pairs[0][0].sequence)
        assert(read2.query_sequence == new_read_pairs[0][1].sequence)
        # Check name and reference id of original reads
        read_name = read1.query_name
        reference_id = read1.reference_id
        # Get position of original mapped reads and number of new reads
        position = "{id}-{start1}-{end1}-{start2}-{end2}".format(
            id=reference_id, start1=read1.reference_start,
            end1=read1.reference_end, start2=read2.reference_start,
            end2=read2.reference_end
        )
        n_pair = len(new_read_pairs)
        # Process new pairs sequentially
        for i, (new_read1, new_read2) in enumerate(new_read_pairs):
            # Generate read pair name
            fastq_name = "{}.{}.{}.{:06d}".format(
                read_name, position, n_pair, i
            )
            # Process each read in pair sequentially
            for old_read, new_read in [(read1, new_read1), (read2, new_read2)]:
                # Extract sequence and quality for read 1
                if old_read.is_reverse:
                    fastq_seq = self.rev_comp(new_read.sequence)
                    fastq_qual = self.quality_string(
                        new_read.quality, reverse=True
                    )
                else:
                    fastq_seq = new_read.sequence
                    fastq_qual = self.quality_string(
                        new_read.quality, reverse=False
                    )
                assert(len(fastq_seq) == len(fastq_qual))
                # Write fastq
                fastq_entry = '@{name}\n{seq}\n+\n{quality}\n'.format(
                    name=fastq_name, seq=fastq_seq, quality=fastq_qual
                )
                self.fastq.write(fastq_entry)


class CreateLog(object):

    def __init__(self, path):
        self.path = path

    def write_counts(
        self, alignment_counts, variant_counts
    ):
        # Get alignment counts log
        alignment_log = (
            'total reads: {total}\n'
            'alignment filter:\n'
            '  secondary: {secondary}\n'
            '  supplementary: {supplementary} \n'
            '  unmapped: {unmapped} \n'
            '  mate umapped: {mate_unmapped}\n'
            '  different chromosomes: {different_chromosomes}\n'
            '  improper pair: {improper_pair}\n'
            '  low mapping quality: {low_mapq}\n'
            '  passed: {passed}\n'
        ).format(**alignment_counts)
        # Get variant counts log
        variant_count_log = (
            'read variant types:\n'
            '  none: {variants_absent}\n'
            '  unknown: {abnormal_alignment}\n'
            '  biallelic snv: {biallelic_snv}\n'
            '  biallelic: {biallelic}\n'
            '  snv: {snv}\n'
            '  mixed: {mixed}\n'
        ).format(**variant_counts)
        # Get allele count log
        allele_count_log = (
            'read allele counts:\n'
            '  reference: {ref_count}\n'
            '  alternative: {alt_count}\n'
            '  other: {other_count}\n'
        ).format(**variant_counts)
        # Get variant filter log
        variant_filter_log = (
            'variant filter:\n'
            '  unwanted variants: {unwanted_variants}\n'
            '  excess variants: {excess_variants}\n'
            '  overlapping variants: {overlapping_variants}\n'
            '  conflicting alleles: {conflicting_alleles}\n'
            '  truncated reads: {too_short}\n'
            '  excess reads: {excess_reads}\n'
            '  to remap: {to_remap}\n'
        ).format(**variant_counts)
        # Write logs to file
        with open(self.path, 'wt') as outfile:
            outfile.write(alignment_log)
            outfile.write(variant_count_log)
            outfile.write(allele_count_log)
            outfile.write(variant_filter_log)


class ProcessAlignments(object):

    def __init__(
        self, bam, vcf, out_prefix, min_mapq, max_vars, max_seqs, min_len,
        only_snv, only_biallelic
    ):
        # Store input parameters
        self.inbam_path = bam
        self.vcf_path = vcf
        self.out_prefix = out_prefix
        self.min_mapq = min_mapq
        self.max_vars = max_vars
        self.max_seqs = max_seqs
        self.min_len = min_len
        self.only_snv = only_snv
        self.only_biallelic = only_biallelic
        # Create output paths
        self.outbam_path = self.out_prefix + '.no_variants.bam'
        self.fastq_path = self.out_prefix + '.allele_flipped.fq.gz'
        self.log_path = self.out_prefix + '.first_alignment_log.txt'
        # Initialise obejcts and open files
        self.bam_generator = BamGenerator(
            self.inbam_path, min_mapq=self.min_mapq
        )
        self.var_tree = vartree.VarTree(self.vcf_path)
        self.outbam = pysam.AlignmentFile(
            self.outbam_path, mode='wb', template=self.bam_generator.bam
        )
        self.fastq = FastqWriter(self.fastq_path)
        self.log = CreateLog(self.log_path)
        # Create counter
        self.counter = {
            'abnormal_alignment': 0,
            'variants_absent': 0,
            'biallelic_snv': 0,
            'biallelic': 0,
            'snv': 0,
            'mixed': 0,
            'ref_count': 0,
            'alt_count': 0,
            'other_count': 0,
            'unwanted_variants': 0,
            'excess_variants': 0,
            'overlapping_variants': 0,
            'conflicting_alleles': 0,
            'too_short': 0,
            'excess_reads': 0,
            'to_remap': 0
        }
        # Create named tuple to contain allele flipped reads
        self.FlippedRead = collections.namedtuple(
            'FlippedRead', ['sequence', 'quality', 'edits']
        )

    def close(
        self
    ):
        # Close all open objects and files
        self.bam_generator.close()
        self.outbam.close()
        self.fastq.close()

    def count_ref_alt_matches(
        self, variants
    ):
        ''' Function counts matches between read alleles and vcf alleles'''
        # Loop through varaints and extract read allele
        for variant in variants:
            try:
                allele_index = variant.alleles.index(variant.read_allele)
            except ValueError:
                allele_index = None
            # Count match
            if allele_index is None:
                self.counter['other_count'] += 1
            elif allele_index == 0:
                self.counter['ref_count'] += 1
            else:
                self.counter['alt_count'] += 1

    def variants_overlap(
        self, variant_list
    ):
        '''Function checks if any variants overlap'''
        for variants in variant_list:
            start_position = 0
            for variant in variants:
                if variant.start < start_position:
                    return(True)
                start_position = variant.end
        return(False)

    def conflicting_alleles(
        self, variant_list
    ):
        """Checks if paired reads have conflicting alleles"""
        # Only process paired reads
        if self.bam_generator.paired:
            common_variants = set(variant_list[0]).intersection(
                set(variant_list[1])
            )
            for variant in common_variants:
                var1 = variant_list[0][variant_list[0].index(variant)]
                var2 = variant_list[1][variant_list[1].index(variant)]
                if var1.read_allele != var2.read_allele:
                    return(True)
        return(False)

    def get_mean_quality(
        self, read
    ):
        # Get qualities
        aligned_base_qualities = read.query_qualities[
            read.query_alignment_start:read.query_alignment_end
        ]
        mean_aligned_base_quality = (
            sum(aligned_base_qualities) // len(aligned_base_qualities)
        )
        return(mean_aligned_base_quality)

    # Function generates allele flipped reads for supplied variants
    def generate_flipped_reads(
        self, read, read_variants
    ):
        """Generate set of reads with all possible combinations"""
        # Create inital read and get mean base quality across aligned segment
        initial_read = self.FlippedRead(
            sequence=read.query_sequence, quality=list(read.query_qualities),
            edits={}
        )
        mean_quality = self.get_mean_quality(read)
        # Create lists to store current and newly generated reads
        current_reads = [initial_read]
        new_reads = []
        # Loop though variants end to start
        read_variants = sorted(
            read_variants, key=lambda x: x.start, reverse=True
        )
        for variant in read_variants:
            # Merge possible reference and alternative alleles
            for new_allele in variant.alleles:
                for old_read in current_reads:
                    # Skip identical alleles
                    old_allele = old_read.sequence[
                        variant.read_start:variant.read_end
                    ]
                    assert(old_allele == variant.read_allele)
                    # Skip '*' marking deletions spanning variant position
                    if new_allele == '*':
                        continue
                    # Create new sequence
                    new_sequence = (
                        old_read.sequence[:variant.read_start] +
                        new_allele +
                        old_read.sequence[variant.read_end:]
                    )
                    # Create new quality
                    if len(old_allele) == len(new_allele):
                        new_quality = old_read.quality
                    else:
                        new_quality = (
                            old_read.quality[:variant.read_start] +
                            [mean_quality] * len(new_allele) +
                            old_read.quality[variant.read_end:]
                        )
                    assert(len(new_sequence) == len(new_quality))
                    # Create new edits
                    new_edits = old_read.edits
                    new_edits[variant.start] = new_allele
                    # Store modified read
                    new_read = self.FlippedRead(
                        sequence=new_sequence, quality=new_quality,
                        edits=new_edits
                    )
                    new_reads.append(new_read)
            # update current reads with new read versions
            current_reads.extend(new_reads)
            new_reads = []
        # Return all reads
        return(current_reads)

    def too_short(
        self, flipped_list
    ):
        '''Determines if flipped reads are of sufficient length'''
        for flipped_reads in flipped_list:
            for read in flipped_reads:
                if len(read.sequence) < self.min_len:
                    return(True)
        return(False)

    def read_pair_combos(
        self, read1_list, read2_list
    ):
        """Generate read pairs with matched edits"""
        # Set output variables
        pair_list = []
        # Loop through all combinations of reads
        for read1, read2 in itertools.product(read1_list, read2_list):
            # Check alleles are identical across common positions
            common = list(read1.edits.keys() & read2.edits.keys())
            if len(common) > 0:
                # Skip read pairs where common positions differ
                read1_alleles = [read1.edits[c] for c in common]
                read2_alleles = [read2.edits[c] for c in common]
                if read1_alleles != read2_alleles:
                    continue
            # Add new alleles and check pait list length
            pair_list.append((read1, read2))
            if len(pair_list) > self.max_seqs:
                break
        # Return pair list and errors
        return(pair_list)

    def filter_reads(
        self
    ):
        # Get variants for each chromosome
        for chromosome in self.bam_generator.chromosomes:
            self.var_tree.read_vcf(chromosome)
            # Loop through reads on chromosome
            for read_list in self.bam_generator.get_reads(chromosome):
                read_no = len(read_list)
                # Get variants for each read...
                try:
                    variant_list = [
                        self.var_tree.get_read_variants(read, partial=True) for
                        read in read_list
                    ]
                # or count abnormal alignments and continue
                except AssertionError:
                    self.counter['abnormal_alignment'] += read_no
                    continue
                # Get variant types
                variant_counts = [len(v) for v in variant_list]
                total_variants = sum(variant_counts)
                all_snv = all(
                    [v.is_snv() for v in itertools.chain(*variant_list)]
                )
                all_biallelic = all(
                    [v.is_biallelic() for v in itertools.chain(*variant_list)]
                )
                # Count variant types...
                if total_variants > 0:
                    if all_snv:
                        if all_biallelic:
                            self.counter['biallelic_snv'] += read_no
                        else:
                            self.counter['snv'] += read_no
                    else:
                        if all_biallelic:
                            self.counter['biallelic'] += read_no
                        else:
                            self.counter['mixed'] += read_no
                # Or save reads with no variants
                else:
                    self.counter['variants_absent'] += read_no
                    for read in read_list:
                        self.outbam.write(read)
                    continue
                # Count reference, alternative and other matches
                for variants in variant_list:
                    self.count_ref_alt_matches(variants)
                # Count reads with unwanted variants
                if not all_snv and self.only_snv:
                    self.counter['unwanted_variants'] += read_no
                    continue
                if not all_biallelic and self.only_biallelic:
                    self.counter['unwanted_variants'] += read_no
                    continue
                # Count reads with excess variants and continue
                if any([count > self.max_vars for count in variant_counts]):
                    self.counter['excess_variants'] += read_no
                    continue
                # Count reads containing overlapping variants and continue
                if self.variants_overlap(variant_list):
                    self.counter['overlapping_variants'] += read_no
                    continue
                # Count paired reads with conflicting alleles and continue
                if read_no == 2:
                    if self.conflicting_alleles(variant_list):
                        self.counter['conflicting_alleles'] += read_no
                        continue
                # Generate flipped reads
                flipped_list = [
                    self.generate_flipped_reads(read, variants) for
                    read, variants in zip(read_list, variant_list)
                ]
                # Count reads with truncated flipped alleles and continue
                if self.too_short(flipped_list):
                    self.counter['too_short'] += len(read_list)
                    continue
                # Process paired reads
                if read_no == 1:
                    # Check read number if not excessive
                    new_reads = flipped_list[0]
                    if len(new_reads) > self.max_seqs:
                        self.counter['excess_reads'] += read_no
                        continue
                    # Write paired reads to file
                    self.counter['to_remap'] += read_no
                    self.fastq.write_fastq(
                        read=read_list[0], new_reads=new_reads
                    )
                else:
                    # Get read pairs and check nor wxcessive
                    new_read_pairs = self.read_pair_combos(*flipped_list)
                    if len(new_read_pairs) > self.max_seqs:
                        self.counter['excess_reads'] += read_no
                        continue
                    # Write paired reads to file
                    self.counter['to_remap'] += read_no
                    self.fastq.write_pair_fastq(
                        read1=read_list[0], read2=read_list[1],
                        new_read_pairs=new_read_pairs
                    )
        # Get counts, adjust and write to file
        alignment_counts = self.bam_generator.counter
        alignment_counts['total'] += self.bam_generator.nocoordinate
        alignment_counts['unmapped'] += self.bam_generator.nocoordinate
        variant_counts = self.counter
        self.log.write_counts(
            alignment_counts=alignment_counts, variant_counts=variant_counts
        )


# Run script
if __name__ == '__main__':
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="Looks for variants (SNPs & indels) overlapping reads. "
        "Three output files are generated with the provided prefix and the "
        "following suffixes. 1) '.invariant.bam' - BAM file containing reads "
        "not overlapping variants. 2) '.remap.fq.gz' - a gzipped FASTQ file "
        "containing the original and variant flipped reads for realignment. "
        "If the input BAM file contains paired end reads then the pairs are "
        "interlevead in the FASTQ file. 3) '.variant_log.txt' text file "
        "containing read processing metrics."
    )
    parser.add_argument(
        "--bam", required=True, help=(
            "Coordinate-sorted and indexed input BAM file."
        )
    )
    parser.add_argument(
        "--vcf", required=True, help=(
            "Coordinate sorted and tabix indexed VCF file."
        )
    )
    parser.add_argument(
        "--out_prefix", required=True, help=(
            "Prefix of output files."
        )
    )
    parser.add_argument(
        "--min_mapq", type=int, default=0, help=(
            "Minimum read mapping quality (default=0)."
        )
    )
    parser.add_argument(
        "--max_vars", type=int, default=6, help=(
            "The maximum number of variants allowed to overlap a read before "
            "discarding the read (default=6)."
        )
    )
    parser.add_argument(
        "--max_seqs", type=int, default=64, help=(
            "The maximum number of allele flipped reads (or read pairs) "
            " to consider remapping. Reads with more allelic combinations "
            "than the specified value are discarded (defult=64)."
        )
    )
    parser.add_argument(
        "--min_len", type=int, default=20, help=(
            "The minimum length of allele flipped reads for remapping "
            "(default=20)"
        )
    )
    parser.add_argument(
        "--snv_only", action='store_true', help=(
            "Discard reads overlapping non snv variants (deafult=False)."
        )
    )
    parser.add_argument(
        "--biallelic_only", action='store_true', help=(
            "Discard reads overlapping non biallelic variants (deafult=False)."
        )
    )
    args = parser.parse_args()
    # Run program
    process_alignments = ProcessAlignments(
        bam=args.bam, vcf=args.vcf, out_prefix=args.out_prefix,
        min_mapq=args.min_mapq, max_vars=args.max_vars, max_seqs=args.max_seqs,
        min_len=args.min_len, only_snv=args.snv_only,
        only_biallelic=args.biallelic_only
    )
    process_alignments.filter_reads()
    process_alignments.close()
