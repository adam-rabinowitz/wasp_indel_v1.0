import collections
import gzip
import intervaltree
import sys
import util


class SNPTable(object):

    def __init__(self):
        self.clear()

    def clear(self):
        self.alleles = None

    def read_file(self, filename):
        """read in SNPs and indels from text input file"""
        # Open input file for chromosome
        try:
            if util.is_gzipped(filename):
                f = gzip.open(filename, "rt")
            else:
                f = open(filename, "rt")
        except IOError:
            sys.stderr.write(
                "WARNING: unable to read from file '%s', "
                "assuming no SNPs for this chromosome\n" %
                filename
            )
            self.clear()
            return
        # Create variables to store variants
        allele_dict = collections.defaultdict(list)
        # Loop through lines of input file
        for line in f:
            # Extract position and variant from line
            pos, a1, a2 = line.strip().split('\t')
            pos = int(pos) - 1
            a1 = a1.upper().replace("-", "")
            a2 = a2.upper().replace("-", "")
            # Check position and variants
            if pos < 0:
                raise ValueError("SNP position < 0: {}".format(line))
            if len(a1) < 1:
                raise ValueError("absent reference allele: {}".format(line))
            # Add interval and alleles to dictionary
            interval = (pos, pos + len(a1))
            alleles = (a1, a2)
            allele_dict[interval].append(alleles)
        f.close()
        # Create interval tree from all potential alleles
        interval_list = [
            intervaltree.Interval(interval[0], interval[1], tuple(alleles))
            for
            interval, alleles in allele_dict.items()
        ]
        self.alleles = intervaltree.IntervalTree(
            interval_list
        )

    def get_overlapping_alleles(self, read):
        # Get position of alleles
        alleles = self.alleles.overlap(
            read.reference_start, read.reference_end
        )
        return(alleles)
