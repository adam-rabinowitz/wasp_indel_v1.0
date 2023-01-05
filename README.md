# WASP-Indel: Indel sensitive pipeline for unbiased read mapping and molecular QTL discovery

## Introduction

WASP-Indel is a suite of tools for unbiased allele-specific read mapping and
discovery of molecular QTLs. WASP-Indel is based on the software tool WASP described in
[the paper](http://biorxiv.org/content/early/2014/11/07/011221): van
de Geijn B\*, McVicker G\*, Gilad Y, Pritchard JK. "WASP:
allele-specific software for robust discovery of molecular
quantitative trait loci".

The main practical advance of WASP-Indel, over the original WASP, is the ability to determine the impact of indels, as well as SNPs, on molecular traits. WASP-Indel also benefits from a simplification of the required input files by eliminating the requirement of converting the input genome FASTA files and VCF to the HDF5 format. There are additional advances in identifying potentially mis-genotyped variants and in the processing of large number of phenotype/genotype combinations.

WASP-Indel has two parts, which can be used independently of each
other:

1. Read filtering tools that correct for biases in allele-specific
   mapping. 

2. A Combined Haplotype Test (CHT) that tests for genetic association
   with a molecular trait using counts of mapped and allele-specific
   reads.

The following directories and files are included with WASP.
Each directory contains its own README file:

* [mapping](./mapping) - Mappability filtering pipeline for correcting allelic mapping biases

* [CHT_input](./CHT_input) - Generate the input for the Combined Haplotype Test from the unbiased allele specific mapping.

* [CHT](./CHT) - Code for running the Combined Haplotype Test



