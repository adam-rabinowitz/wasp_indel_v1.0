# Pipeline for mappability filtering

## Contents

This directory contains scripts that can be used to eliminate mapping bias from mapped allele-specific reads. First, reads are mapped normally using a mapper chosen by the user. Then mapped reads that overlap variants are identified. Allele-flipped versions of the overlapping reads are generated and re-mapped to the genome. Re-mapped reads that fail to map to exactly the same location in the genome are discarded. The required input for the mappability filtering is FASTQ sequence data and a sorted and TABIX indexed VCF file containing the variant information.</br></br> 

## Steps to generate unbiased allele specific alignments

The following steps should be taken to generate an unbiased allele-specific BAM file:</br></br>

1. Align FASTQ file(s) to reference genome and generate a coordinate sorted BAM file from the initial alignment. Can use single-end or paired-end FASTQ files but not a mix of the two. The BAM file should not be filtered and is expected to contain all the reads from the orginal FASTQ files. If using BWA for alignment use the -M flag to flag shorter split hits as secondary.</br></br>
   
2. Parse reads in the initial BAM file using the script 'generate_variant_reads.py'. The script requires three parameters which are a BAM file (from step 1), the variants (in VCF format) and the prefix of the output files which are specified with the --bam, --vcf and --out-prefix arguments. The following alignments are discarded in the inital processing of the BAM file: unmapped reads, secondary alignments, supplementary alignments and alignments with a mapping quality below the value specified with the argument --min-mapq. Additionally, paired-alignments aligned to different chromosomes or forming an improper pair are also discarded. For paired-end reads if one alignment in the pair is discarded then the other alignment is discarded too.</br></br>The remaining alignments are then split by whether they overlap a variant in the supplied VCF file. Alignments not overlapping variants are written to a BAM file with the suffix '.invariant.bam'. For alignments overlapping variants then the original read sequence is written to a FASTQ file, with the suffix '.remap.fq.gz', along with all possible allele-flipped versions of the read(s). For paired-end alignments the read-pairs are interleaved in this FASTQ file.</br></br>The default activity of the script is to generate allele-flipped version of the read(s) containing all possible combinations of variants within the VCF. To limit processing time the maximum number of variants overlapping each read and the maximum number of allele-flipped reads (or read-pairs) can be speciffied using the arguments --max-vars and --max-seqs, respectively. If theses values are exceeded the read (or read-pair) is discarded. The allele flipping can be done in a phase specific manner if the samples within the VCF are specified using the "--samples" argument. If samples are specified the script checks to determine if the VCF contains phased data for relevant samples. This check for phased data can be ignored by adding the --assume-phased flag.</br></br>

3. Realign the reads in the '.remap.fq.gz' to the reference genome and generate a name sorted BAM from the initial alignment. The BAM file should not be filtered and is expected to contain all the reads from the allele-flipped FASTQ file. If using BWA for alignment use the -M flag to flag shorter split hits as secondary.</br></br>

4. Parse reads in the allele-flipped BAM file using the script 'filter_remapped_reads.py'. Script checks that all allele-flipped versions of the original sequence map to the same location as the orginal sequence. Reads with inconsistent mapped locations are discarded while the original read of consistently mapped alignments are written to a BAM file with the suffix '.consistent.bam'. The argument "--wobble" allows a specified distance between the originally mapped reads and the remapped flipped reads. The default is to allow zero wobble.</br></br>

5. Merge the '.invariant.bam' BAM file from step 2 and the '.consistent.bam' BAM file from step 4 and coordinate sort the merged BAM file.</br></br>

6. Remove duplicates from the merged BAM file. **Importantly**, to ensure the reference allele is not favoured over alternative alleles the reads to be kept should be randomly selected from the duplicate set. This random removal of duplicates can be specified in the MarkDuplicates tool from Picard by supplying the value "RANDOM" to the --DUPLICATE_SCORING_STRATEGY argument.</br></br>

7. This step collates metrics from the unbiased and deduplicated BAM file (from step 6) for the variants within the VCF file. This step is only required if you will go on to perform the combined haplotype test. The metrics are generated using the script 'get_counts.py' and the output metrics file is bgzip compressed and tabix indexed. The script generates metrics for an individual sample within the VCF file specified with the argument --sample flag.</br></br>The output text file contains 15 columns. The first 6 columns are
chromosome, position, variant ID, reference allele, alternative allele and haplotype. Columns 7-9 are the probabilities that the genotype is homozygous reference allele, heterzoygous, and homozygous alternative allele, respectively. The data for columns 1-9 are extracted from the VCF file. Columns 10-12 contain the allele specific counts for the reference, alternative and other alleles. The allele specific counts are generated by randomly selecting a biallelic heterozygous variant contained within each alignment and determing the allele identity. Columns 13-15 contain the total reference, alternative and other alleles counts. These counts are generated by determining the variant allele for all variants contained in each alignment. Columns 10-12 are used directly in the combined haplotype test while columns 13-15 can be used to correct for genotyping errors.</br></br>

## Dependencies

The scripts in the mapping portion of the WASP-indel package have been developed with following python versions and packages:

* python version 3.6

* [intervaltree](https://github.com/chaimleib/intervaltree) version 3.1.0

* [pysam](https://github.com/pysam-developers/pysam) version 0.15.3

It is likely that alternative versions of these tools will work but they have not been tested.