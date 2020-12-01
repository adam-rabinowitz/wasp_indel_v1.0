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
            pos, ref, alt = line.strip().split()
            pos = int(pos) - 1
            ref = ref.upper().replace("-", "")
            alt = alt.upper().replace("-", "")
            # Check position and variants
            if pos < 0:
                raise ValueError("SNP position < 0: {}".format(line))
            if len(ref) < 1:
                raise ValueError("absent reference allele: {}".format(line))
            # Add interval and alleles to dictionary
            interval = (pos, pos + len(ref))
            if interval in allele_dict:
                raise ValueError('duplicated reference bases')
            else:
                allele_dict[interval] = (ref, alt)
        f.close()
        # Create interval tree from all potential alleles
        interval_list = [
            intervaltree.Interval(start, end, (ref, alt)) for
            (start, end), (ref, alt) in allele_dict.items()
        ]
        self.alleles = intervaltree.IntervalTree(
            interval_list
        )

    def get_overlapping_alleles(self, read):
        # Get intervals contained entirely within mapped positions
        alleles = self.alleles.envelop(
            read.reference_start, read.reference_end
        )
        # Convert to tuples sort and return
        allele_tuple = [
            (x.begin, x.end, x.data[0], x.data[1]) for x in alleles
        ]
        allele_tuple = sorted(allele_tuple)
        return(allele_tuple)
