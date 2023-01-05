import argparse
import gzip
import math
import numpy as np
import pysam
from scipy.optimize import fmin


def open_files(file_list, mode):
    files = []
    for filename in file_list:
        if filename.endswith('.gz'):
            files.append(gzip.open(filename, mode))
        else:
            files.append(open(filename, mode))
    return(files)


def get_gc_ratio(fasta, chrom, starts, ends):
    # Set intial counts
    at_count, gc_count = 0, 0
    # Loop trhough starts and ends
    for start, end in zip(starts, ends):
        # Get sequence
        sequence = fasta.fetch(
            reference=chrom, start=start, end=end
        )
        sequence = sequence.upper()
        # Add counts
        at_count += sequence.count('A') + sequence.count('T')
        gc_count = sequence.count('C') + sequence.count('G')
    # Return base counts
    gc_ratio = gc_count / (at_count + gc_count)
    return(gc_ratio)


def get_region_tuple(lines):
    # Get regions and check identity
    region_tuples = [(line[0], line[7], line[8]) for line in lines]
    for rt in region_tuples[1:]:
        if rt != region_tuples[0]:
            raise ValueError('input files have differing regions')
    # Return region tuple
    region_tuple = region_tuples[0]
    return(region_tuple)


def load_data(inlist, fasta_path, min_counts):
    # Open files and skip header
    infiles = open_files(inlist, "rt")
    _ = [file.readline() for file in infiles]
    # Open fasta file
    fasta = pysam.FastaFile(fasta_path)
    # Process files
    region_counts = {}
    fit_table = []
    # Loop through files
    while True:
        # Read lines from files
        lines = [file.readline().strip().split() for file in infiles]
        if not lines[0]:
            assert not any(lines)
            break
        # Get region data
        region_tuple = get_region_tuple(lines)
        # Skip region if counts have been previously generated
        if region_tuple in region_counts:
            continue
        # Get gc ratio
        chrom = region_tuple[0]
        starts = [int(s) - 1 for s in region_tuple[1].split(';')]
        ends = [int(e) for e in region_tuple[2].split(';')]
        try:
            gc_ratio = get_gc_ratio(
                fasta=fasta, chrom=chrom, starts=starts, ends=ends
            )
        except ZeroDivisionError:
            region_counts[region_tuple] = None
            continue
        # Get counts for each file and the total
        try:
            line_counts = [int(c[15]) for c in lines]
        except ValueError:
            region_counts[region_tuple] = None
            continue
        count_total = sum(line_counts)
        # Store counts above minimum
        if count_total >= min_counts:
            count_list = [gc_ratio, count_total] + line_counts
            fit_table.append(count_list)
            region_counts[region_tuple] = count_list
        # Or dicard counts
        else:
            region_counts[region_tuple] = None
    # Close input file
    for infile in infiles:
        infile.close()
    # Convert fit table to array and return
    fit_table = np.array(fit_table, dtype=np.float64)
    return(region_counts, fit_table)


def read_splines(fit_in_file):
    spline_file = open(fit_in_file, "rt")
    coefs = []
    for line in spline_file:
        coefs.append([float(x) for x in line.strip().split()])
    spline_file.close()
    return coefs


def write_splines(coefs, fit_out_file):
    spline_file = open(fit_out_file, "wt")
    for ind in coefs:
        spline_file.write("\t".join([str(x) for x in ind]) + "\n")
    spline_file.close()
    return


def fit_splines(count_table):
    coefs = []
    # Get column sums
    col_sums = np.sum(count_table, 0)
    # Check column shape
    if len(count_table.shape) == 1:
        raise ValueError('count table contains no data')
    for col in range(2, count_table.shape[1]):
        slope_start = col_sums[col] / col_sums[1]
        print(str(slope_start))

        arg_list = (count_table[:, 0], count_table[:, 1],
                    count_table[:, col])
        slope = fmin(splineone, [slope_start], args=arg_list,
                     maxiter=500000, maxfun=500000)

        # Fit a slope only model (total read depth is different
        # between individuls within an individual reads per region
        # is linear
        coefs_slope = [0, 0, 0, 0, 0, 0, slope[0], 0, 0, 0]
        print(str(coefs_slope))
        # Read shift only (allowing more or fewer reads to fall
        # in peaks)
        arg_list = ([0, 0, 0, 0, 0], count_table[:, 0],
                    count_table[:, 1], count_table[:, col])
        coefs_depth = fmin(splineread, [0, slope[0], 0, 0, 0],
                           args=arg_list, maxiter=500000,
                           maxfun=500000)
        print(str(coefs_depth))

        # GC correct and total count only
        arg_list = ([0, slope[0], 0, 0, 0], count_table[:, 0],
                    count_table[:, 1], count_table[:, col])

        coefs_gc = fmin(splinegc, [0, 0, 0, 0, 0],
                        args=arg_list, maxiter=500000, maxfun=500000)
        print(str(coefs_gc))

        # GC correct plus count shift
        arg_list = (count_table[:, 0], count_table[:, 1],
                    count_table[:, col])
        coefs_full = fmin(splinefit,
                          list(coefs_gc) + list(coefs_depth),
                          args=arg_list,
                          maxiter=500000, maxfun=500000)
        print(str(coefs_full))
        coefs.append(coefs_full)

    return(coefs)


def splineone(slope, x1, x2, y):
    return splinefit([0, 0, 0, 0, 0, 0, slope[0], 0, 0, 0], x1, x2, y)


def splinegc(gcargs, otherargs, x1, x2, y):
    a = np.hstack((gcargs, otherargs))
    return splinefit(a, x1, x2, y)


def splineread(otherargs, gcargs, x1, x2, y):
    a = np.hstack((gcargs, otherargs))
    return splinefit(a, x1, x2, y)


def splinefit(arg, x1, x2,  y):
    expecteds = [
        max(0.000001, x) for x in calc_adjusted_totals(x1, x2, arg)
    ]
    # resids = y - expecteds
    # sse = sum(resids**2)
    loglike = -sum(y * [math.log(x) for x in expecteds] - expecteds)
    return loglike


def array_exp(x):
    if hasattr(x, "__len__"):
        return [math.exp(i) for i in x]
    else:
        return math.exp(x)


def calc_adjusted_totals(gc, tot, coefs):
    return array_exp(coefs[0] + gc*coefs[1] + gc**2*coefs[2]
                     + gc**3*coefs[3] + gc**4*coefs[4]) * \
                     (0 + tot*coefs[6] + tot**2*coefs[7] +
                      tot**3*coefs[8] + tot**4*coefs[9])


# Update read depths
def update_totals(
    inlist, outlist, counts, coefs_table
):
    # Create dictionary to store counts
    adj_totals = {}
    # Open files
    infiles = open_files(inlist, "rt")
    outfiles = open_files(outlist, "wt")
    # Write the header line
    for ind in range(len(infiles)):
        line = infiles[ind].readline()
        outfiles[ind].write(line)
    # Loop through files
    while True:
        # Read lines from files
        lines = [file.readline().strip().split() for file in infiles]
        if not lines[0]:
            assert not any(lines)
            break
        # Get region tuple and counts
        region_tuple = get_region_tuple(lines)
        region_counts = counts[region_tuple]
        # Skip lines for which there are no counts
        if region_counts is None:
            continue
        # Get adjusted counts or...
        try:
            adj_tot_list = adj_totals[region_tuple]
        # calculate adjusted counts
        except KeyError:
            adj_tot_list = []
            for ind in range(len(infiles)):
                # Get adjusted total and store
                adj_tot = calc_adjusted_totals(
                    gc=region_counts[0],
                    tot=region_counts[1],
                    coefs=coefs_table[ind]
                )
                adj_tot_list.append(max(adj_tot, -1000000))
            adj_totals[region_tuple] = adj_tot_list
        # Write adjusted totals to file
        for ind in range(len(infiles)):
            line = lines[ind]
            adj_tot_ind = adj_tot_list[ind]
            outfiles[ind].write(
                '\t'.join(line[:16]) + '\t' + str(adj_tot_ind) + '\n'
            )
            outfiles[ind].flush()
    # Close files
    for infile in infiles:
        infile.close()
    for outfile in outfiles:
        outfile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Adjusts the expected read depths for each "
        "region in the CHT input files. Expected read depths take "
        "into account the GC content and total read depth of each "
        "region (accross individuals). The adjusted read depths "
        "are computed by fitting quartic functions for these values "
        "for each individual."
    )
    parser.add_argument(
        "--infiles", nargs="+", required=True, help="Input CHT files"
    )
    parser.add_argument(
        "--outfiles", nargs="+", help="Output CHT files"
    )
    parser.add_argument(
        '--fasta', required=True, help='faidx indexed fasta file'
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--fit_infile', default=None, help=(
            'Read coefficients from specified file instead of estimating them.'
        )
    )
    group.add_argument(
        "--fit_outfile", default=None, help=(
            "Estimate coefficients and write them to specified file, "
            "but do not adjust read counts."
        )
    )
    parser.add_argument(
        "--min_counts", type=int, default=0, help=(
            "only use rows with at least min_counts for fitting"
        )
    )
    parser.add_argument(
        "--sample", type=int, default=100000, help=(
            "randomly sample this many rows and use them for fitting "
            "coefficients. Specify 0 if all rows are to be used. "
            "(default=100000)"
        )
    )
    args = parser.parse_args()
    # Check arguments
    if args.outfiles:
        if len(args.infiles) != len(args.outfiles):
            raise ValueError('must be one outfile for each infile')
    # Create count table
    print("Loading count table")
    counts, fit_table = load_data(
        inlist=args.infiles, fasta_path=args.fasta, min_counts=args.min_counts
    )
    # Read fit from file or generate new fit
    print("Generating coefficients")
    if args.fit_infile:
        coefs_list = read_splines(args.fit_infile)
    else:
        # Sample data for the fit
        if args.sample:
            if args.sample < fit_table.shape[0]:
                rows = np.arange(fit_table.shape[0])
                np.random.seed(42)
                sampled_rows = np.random.choice(
                    rows, size=args.sample, replace=False
                )
                fit_table = fit_table[sampled_rows]
        # Perform fit and write to file
        coefs_list = fit_splines(fit_table)
        write_splines(coefs_list, args.fit_outfile)
    # or update total
    if args.outfiles:
        print("Updating depth")
        update_totals(
            args.infiles, args.outfiles, counts, coefs_list
        )
