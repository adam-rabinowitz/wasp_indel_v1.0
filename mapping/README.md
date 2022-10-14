Pipeline for mappability filtering
==================================

This directory contains scripts that can be used to eliminate mapping bias from
mapped allele-specific reads.  First, reads are mapped normally using a mapper
chosen by the user (must output BAM or SAM format). Then mapped reads that
overlap variants are identified. Allele flipped versions of the overlapping
reads are generated and re-mapped to the genone. Re-mapped reads that fail to
map to exactly the same location in the genome are discarded.

SnakeMake
---------

We now provide a Snakemake workflow that can be used to run the entire
mappability filtering pipeline along with a config file containing the parameters that should be supplied for the Snakemake workflow to function. Additionally, there is a environment yml file for instaling the required packages using conda.

Step 1 & 2:
-----------

Initial alignment of FASTQ file(s) to refrence and generation of a coordinate
sorted BAM file. Can use single-end or paired-end FASTQ files but not a mix of
the two. The BAM file should not be filtered and is expected to contain all the
reads from the orginal FASTQ files.

If using BWA for alignment use the -M flag to flag shorter split hits as
secondary.

Step 3:
-------

Parse reads in the initial BAM file using the script
'generate_variant_reads.py'. The following alignments are discarded: unmapped
reads, secondary alignments, supplementary alignments and alignments with a
mapping quality below the specified value. For paired-end reads if one alignment
in the pair is discarded then the other alignment is discarded too.
Additionally, paired-alignments aligned to different chromosomes or forming an
improper pair are also discarded.

The remaining alignments are then split by whether they overlap a variant in the
supplied VCF file. Alignments not overlapping variants are written to a BAM file
with the suffix '.invariant.bam'. For alignments overlapping variants then the
original read sequence is written to a FASTQ file, with the suffix
'.remap.fq.gz', along with all possible allele-flipped versions of the read(s).
For paired-end alignments the read-pairs are interleaved in this FASTQ file.

Please note the allele flipping can be done in a phase specific manner if samples within the VCF are specified using the "--samples" argument. If samples are specified the script checks to determine if the VCF contains phased data for relevant samples. This check for phased data can be ignored by adding the --assume phased flag.

Step 4 & 5:
-----------

Second alignment of allele-flipped FASTQ file to reference and generation of a
name sorted BAM file. The BAM file should not be filtered and is expected to
contain all the reads from the allele-flipped FASTQ file.

Step 6:
-------

Parse reads in the allele-flipped BAM file using the script
'filter_remapped_reads.py'. Script checks that all allele-flipped versions of
the original sequence map to the same location as the orginal sequence. Reads
with inconsistent mapped locations are discarded while the original read of
consistently mapped alignments are written to a BAM file with the suffix
'.consistent.bam'.

The argument "--wobble" allows a spcified by distance between the originally mapped reads and the remapped flipped reads. The default is to allow zero wobble. 

Step 7 & 8:
-----------

Merge the '.invariant.bam' BAM file from step 3 and the '.consistent.bam' BAM
file from step 6. Coordinate sort the merged BAM file.

Step 9:
-------

Remove duplicates from the merged BAM file. To ensure the reference allele is
not favoured over alternative alleles the reads to be kept should be randomly
selected from the duplicates. This random removal of duplicates can be specified
in the MarkDuplicates tool from Picard.

Step 10:
--------

Generate variant metrics from the filtered and deduplicated BAM file using the
script 'get_counts.py'. The count metrics file is bgzip compressed and tabix
indexed. This tabix file is a key input for generating data for the CHT test.
