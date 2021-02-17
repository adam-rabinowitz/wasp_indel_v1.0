Pipeline for generating CHT input
=================================

This directory contains scripts that can be used to generate input files for the combined haplotype test (CHT).

Step 1:
-------

Identify test variants and target regions using the script
'get_target_regions.py'. The script requires two classes of input files: 1)
--variants - variant count files generated using the 'get_counts.py' script;
2) --bams - alignment files from which the variant count files were generated.
There should be a single variant count file and BAM file for each
sample/individual in the study.

There are four parameters for filtering the test variants and associated target
regions: 1) --min_het - mimimum number of heterozygous test variantss across all
samples; 2) --min_minor - minimum number of minor allele counts for test variant
across all samples; 3) --min_as_reads - minimum number of allele specific reads
for test region across all samples; 4) --min_total_reads - minimum number of
total reads for test region across all samples.

The test variants and their associated target regions can be defined using two
seperate approaches. If a window is specified via the --window argument then the
test regions are defined as the interval +/- the specified window around the
test variant. Alternatively the locations in which to identify accetable test
variants and target regions can be supplied directly via the --regions argument.
The text file should contain 5 tab seperated columns: 1) chromosome; 2) start
and end of the region in which to identify test variants. 4 & 5) start and end
of test region. The coordinates are 1-based so that the first 100 bases of a
chromosome would be specified with start and end of 1 and 100.

Step 2:
-------

Generate the CHT input file for a single sample/individual using the
'get_region_data.py' script. The script requires four arguments: 1) --variants -
variant count file fot the individual; 2) --bam alignment file from which the
variant count file was generated; 3) --regions - the file containining test
variants and test regions generated in Step 1; 4) --outfile - The output file
that can be used as CHT input.

The heterozygous probabilities of varaints within the target regions can be
corrected using one of two sets of counts. If --adjust_hetp is specified as
'as_counts' then the probabilities are corrected using only allele specific
counts. If --adjust_hetp is specified as 'total_counts' then the probabilities
are corrected using total allele counts. 

The output CHT file is specified using the --outfile argument.
