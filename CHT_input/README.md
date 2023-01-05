# Pipeline for generating input for CHT

## Contents

This directory contains scripts that can be used to generate the input files for the combined haplotype test. The required input for these steps are the unbiased BAM files (generated in step 6 of the mapping pipeline), variant count files (generated in step 7 of the mapping pipeline) and an indexed genome FASTA file. There should be a single BAM file and count file for each individual/genotype in the study.</br></br>

The combined haplotype type requires a set of query variants and, for each variant, a test region in which to perform the combined haplotype test. The read depth and allelic imbalance in the test region is then associated with the genotype of the query variant. The query variants will be a subset of those contained in the VCF file used to generate the variant count files in the mapping pipeline. As explained below, the associated test regions can be set as the region surrounding the query variants or defined manually.</br></br> 

## Steps to generate the combined haplotype test input

The following steps should be taken to generate the input data for the combined haplotype test. The scripts in steps 1 and 3 are run using the data of all individuals/genotypes in combination. The script in step 2 is run on each individual/genotype in isolation.</br></br>

1. Identify potentially informative variants and associated test regions using the script 'get_target_regions.py'. The script requires variant count files generated in the mapping pipeline along with the BAM files from which the variant count files were generated. These files are supplied with the --variants and --bams arguments, respectively. The user must also define the variants and their associated test regions. The variants and regions can be defined in one of two ways. If a window is specified via the --window argument, then the putative query variants are all the variants within the variant count files and the associated target regions is the interval +/- the specified window around the query variants. Alternatively, regions in which to extract query variants and their target regions can be specified by providing a text file to the --regions argument. The text file should contain 5 tab seperated columns: chromosome, start of region in which to identify query variants, end of region in which to identify query variants, start of test region, and end of test region. The script will then associate each query variant found in the specifoed query variant region with the same test region. The coordinates in this file are 1-based so that the first 100 bases of a chromosome would be specified with start and end of 1 and 100. This region file should not contain a header. </br></br>There are four parameters for filtering the query variants and associated test regions. First, is the mimimum number of heterozygous genotypes of the variant across all genotypes in the study. Second, is the minimum number of minor allele counts for the variant across all individuals. Third, is the minimum number of allele specific reads for the test region across all individuals. Fourth, is the minimum number of total reads for test region across all individuals. These values can be supplied by the the --min-het, min-minor, --min-as-reads and --min-total-reads arguments, respectively.</br></br>The output of this script, is a text file containing the following six columns: chromsome, position of variant, reference allele of variant, alternative allele of variant, start of test region, and end of test region. Again, the coordinates are 1 based.</br></br>

2. Generate the CHT input file for a single individual/genotype using the 'get_region_data.py' script. The script requires the variant count file and BAM alignment file generated in the mapping workflow and used in step 1 above. Additionaly, the query variants and associated test regions, generated in the first step, should be supplied with the --regions argument.</br></br>The heterozygous probabilities of varaints within the target regions are corrected using one of two methodolgies. If the --adjust-hetp argument is specified as 'as_counts' then the probabilities are corrected using only allele specific counts. If --adjust-hetp is specified as 'total_counts' then the probabilities are corrected using total allele counts. The original version of WASP only used the allele specific counts to correct the heterozygous probabilities but this disregards additional information present in the BAM files. We therefore recommend adjusting heterozygous probabilities using total counts.</br></br>The output text file contains 17 columns and these are all explained below. **Notes** The text "SNP" in the column headers is a carry over over from the original version of WASP. The query variants and region variants can be indels and multinucleotide polymorhisms as well as SNPs. If the query variant is homozygous for the particular genotype/individual then all the region haplotype counts will be zero. If the variant in the test region is homozygous for the particular genotype/individual then the region haplotype counts will be zero for that region variant. Linkage probability is always set to 1.</br>
    - CHROM - Chromosome of the query variant and test region
    - TEST.SNP.POS - Position of the query variant
    - TEST.SNP.ID - Identifier of the query variant
    - TEST.SNP.REF.ALLELE - Reference allele of the query variant
    - TEST.SNP.ALT.ALLELE - Alternative allele of the query variant
    - TEST.SNP.GENOTYPE - Genotype of the query variant
    - TEST.SNP.HAPLOTYPE - Haplotype of the query variant
    - REGION.START - Start of the test region
    - REGION.END - End of the test region
    - REGION.SNP.POS - Position of the variants in the test region
    - REGION.SNP.HET.PROB - Probability of the test region variants being heterozygous
    - REGION.SNP.LINKAGE.PROB - Probability of the test region variants being linked to query region.
    - REGION.SNP.REF.HAP.COUNT - Count of allele specific reads of the test region variants that are the same haplotype as the reference allele of the query variant 
    - REGION.SNP.ALT.HAP.COUNT - Count of allele specific reads of the test region variants that are the same haplotype as the alternative allele of the query variant 
    - REGION.SNP.OTHER.HAP.COUNT - Count of allele specific reads of the test region variants that are of an unexpected genotype.
    - REGION.READ.COUNT - Total number of reads spanning the test region.
    - GENOMEWIDE.READ.COUNT - Total number of reads in the genome.</br></br>

3. The final step in generating the CHT input files is updating the read counts of the test regions (column 16 of the file outlined above) taking into account the regions GC content and the total counts for each individual/genotype. The GC content is extracted from the indexed FASTA file supplied using the --fasta argument. The script takes in all the CHT input files generated in step 2 and supplied using the --infiles argument. A single output file is generated for each of the input files and these output files are specified using the --outfiles argument. The order of the indivduals/genotypes in the input files and output files should be identical.</br></br>The script generates a fit between GC content and read count and then adjusts the input counts according to this fit. To save the fit to file one can supply a file path to the --fit-outfile argument. A previously generated fit can be supplied to the script using the --fit-infile argument. Two further arguments can help to speed up the fitting. The argument --sample allows the specification of a finite number of regions to be sampled for calculation of the fit. Supplying a value to the argument --min-counts will cause regions whose total counts, across all samples, is below the specified value to be skipped.</br></br>The output file for this step is of the same format as the previous step but the values in column 16 are updated and column 17 is removed. This file can be used as the input for the combined haplotype test.   

## Dependencies

The scripts in the mapping portion of the WASP-indel package have been developed with following python versions and packages:

* python version 3.6

* [intervaltree](https://github.com/chaimleib/intervaltree) version 3.1.0

* [pysam](https://github.com/pysam-developers/pysam) version 0.15.3

* [numpy](https://numpy.org/) version ?

* [scipy](https://scipy.org/) version ?

It is likely that alternative versions of these tools will work but they have not been tested.
