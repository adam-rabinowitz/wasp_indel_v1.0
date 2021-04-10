import collections
import pysam


class FirstBamGenerator():

    def __init__(self, bam, min_mapq):
        # Store parameters and open BAM file
        self.bam_path = bam
        self.bam = pysam.AlignmentFile(self.bam_path)
        self.min_mapq = min_mapq
        # Check pairing
        paired = [read.is_paired for read in self.bam.head(100)]
        if any(paired):
            assert(all(paired))
            self.paired = True
        else:
            self.paired = False
        # Get genome and alignment metrics
        self.index_statistics = self.bam.get_index_statistics()
        self.nocoordinate = self.bam.nocoordinate
        # Create counter
        self.counts = collections.Counter()

    def abnormal_cigar(self, read):
        if 'I' in read.cigarstring:
            if 'M' not in read.cigarstring:
                return(True)
            if read.cigarstring.index('I') < read.cigarstring.index('M'):
                return(True)
            if read.cigarstring.rindex('I') > read.cigarstring.rindex('M'):
                return(True)
        return(False)

    def get_reads(self, chromosome):
        # Create cache to store unfound reads
        read_pair_cache = {}
        # Loop though reads on chromosome
        for read in self.bam.fetch(contig=chromosome):
            # Count total reads
            self.counts['total_align'] += 1
            # Count and skip secondary reads
            if read.is_secondary:
                self.counts['secondary_align'] += 1
                continue
            # Count and skip supplementary reads
            if read.is_supplementary:
                self.counts['supplementary_align'] += 1
                continue
            # Count and skip unmapped reads
            if read.is_unmapped:
                self.counts['unmapped_align'] += 1
                continue
            # Further process paired end reads
            if self.paired:
                assert(read.is_paired)
                # Count and skip unmapped mates
                if read.mate_is_unmapped:
                    self.counts['mate_unmapped_align'] += 1
                    continue
                # Count and skip mates mapped to different chromosome
                if read.reference_id != read.next_reference_id:
                    self.counts['different_chrom_align'] += 1
                    continue
                # Count and skip imporperly paired mates
                if not read.is_proper_pair:
                    self.counts['improper_pair_align'] += 1
                    continue
                # Get paired reads if mate has been cached...
                if read.query_name in read_pair_cache:
                    cached_read = read_pair_cache.pop(read.query_name)
                    if read.is_read1:
                        assert(cached_read.is_read2)
                        mapped_reads = [read, cached_read]
                    else:
                        assert(cached_read.is_read1)
                        assert(read.is_read2)
                        mapped_reads = [cached_read, read]
                # ...or store first in pair
                else:
                    read_pair_cache[read.query_name] = read
                    mapped_reads = None
            # Process single end reads
            else:
                assert(not read.is_paired)
                mapped_reads = [read]
            # Further filter properly mapped reads
            if mapped_reads:
                # Get mapping quality of reads
                mapq = [read.mapping_quality for read in mapped_reads]
                # Count and skip reads with low mapping quality
                if min(mapq) < self.min_mapq:
                    self.counts['low_mapq_align'] += len(mapped_reads)
                    continue
                # Check cigarstring of reads
                abnormal_cigar = [
                    self.abnormal_cigar(read) for read in mapped_reads
                ]
                # Count and skip reads with abnormal cigars
                if any(abnormal_cigar):
                    self.counts['abnormal_cigar_align'] += len(mapped_reads)
                    continue
                # ...or count and yield reads with acceptable quality
                else:
                    self.counts['passed_align'] += len(mapped_reads)
                    yield(mapped_reads)
        # Check read pair cache is empty
        assert(len(read_pair_cache) == 0)

    def close(self):
        self.bam.close()
