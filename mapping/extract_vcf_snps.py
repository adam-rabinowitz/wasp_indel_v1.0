import argparse
import gzip
import os


def extract_variant(
    vcf, outdir
):
    # Set variable for processing vcf
    outfile = None
    current_chrom = None
    previous_chroms = set()
    # Open input vcf and loop through file
    infile = gzip.open(vcf, 'rt')
    for line in infile:
        # Skip comments
        if line.startswith('#'):
            continue
        # Extract variant data
        line_data = line.strip().split('\t')
        chrom, position = line_data[0:2]
        ref, alt = line_data[3:5]
        # Process novel chromsomes
        if chrom != current_chrom:
            # Check new chromosome has not been previously observed
            current_position = 0
            previous_chroms.add(current_chrom)
            current_chrom = chrom
            if current_chrom in previous_chroms:
                raise ValueError('vcf is not sorted by chromosome')
            # Close previous file and create new one
            if outfile is not None:
                outfile.close()
            outpath = os.path.join(outdir, chrom + '.txt.gz')
            outfile = gzip.open(outpath, 'wt')
        # Check position
        if int(position) < current_position:
            raise ValueError('vcf is not sorted by position')
        current_position = int(position)
        # Write variant to output file
        output = '\t'.join([position, ref, alt]) + '\n'
        outfile.write(output)
    # Close files
    infile.close()
    if outfile is not None:
        outfile.close()


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(
        description='creates chromosome specific variant files from vcf'
    )
    parser.add_argument('vcf', help='path to vcf file')
    parser.add_argument('outdir', help='path to output directory')
    args = parser.parse_args()
    # Extract chromosome specific variants
    extract_variant(vcf=args.vcf, outdir=args.outdir)
