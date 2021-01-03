Pipeline for mappability filtering
==================================

This directory contains scripts that can be used to eliminate mapping
bias from mapped allele-specific reads.  First, reads are mapped
normally using a mapper chosen by the user (must output BAM or SAM
format). Then mapped reads that overlap variants are identified. 
Allele flipped versions of the overlapping reads are generated and
re-mapped to the genone. Re-mapped reads that fail to map to exactly the
same location in the genome are discarded.

SnakeMake
---------

We now provide a Snakemake workflow that can be used to run the entire
mappability filtering pipeline.

Step 1:
-------

Step 2:
-------

Things to do:

complete documentation
create test suite
think about adding indel realignment
