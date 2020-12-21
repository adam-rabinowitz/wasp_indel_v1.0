import collections
import gzip
import intervaltree
import os
import sys


class SNPTable(object):

    def __init__(self):
        self.clear()

    def clear(self):
        self.alleles = None

    def read_file(self, filename):
        """read in SNPs and indels from text input file"""
        # Open file if it exists...
        if os.path.isfile(filename):
            if filename.endswith('.gz'):
                f = gzip.open(filename, "rt")
            else:
                f = open(filename, "rt")
        # or print warning message if it doesn't
        else:
            warning = (
                "WARNING: cannot find file '{}', assuming no SNPs "
                "for this chromosome\n"
            ).format(filename)
            sys.stderr.write(warning)
            self.clear()
            return
        # Create variant tuple
        variant = collections.namedtuple(
            'variant', ['start', 'end', 'ref', 'alt']
        )
        # Create variables to store variants
        interval_list = []
        # Loop through lines of input file
        for line in f:
            # Extract variant start and end and ref and alt alleles
            start, ref, alt = line.strip().split()
            start = int(pos) - 1
            ref = ref.upper().replace("-", "")
            alt = alt.upper().replace("-", "")
            end = start + len(ref)
            # Check position and variants
            if start < 0:
                raise ValueError("variant start < 0: {}".format(line))
            if len(ref) < 1:
                raise ValueError("absent reference allele: {}".format(line))
            # Create intervaltree interval and add to list
            interval = intervaltree.Interval(
                start, end, variant(start, end, ref, alt)
            )
        f.close()
        # Create intervaltree IntervalTree from all intervals
        self.alleles = intervaltree.IntervalTree(interval_list)

    def get_overlapping_variants(self, start, end):
        # Get read positions and overlapping allele)
        intervals = self.alleles.overlap(start, end)
        variants = [x.data for x in intervals]
        return(variants)


#             # Count the number of ends found
#             found = [x is not None for x in (read_begin, read_end)]
#             found_count = sum(found)
#             # Break and store error where only 1 end is localised
#             if found_count == 1:
#                 error = 'partial'
#                 variants.clear()
#                 break
#             # Store variants within read
#             variants[begin] = (read_begin, read_end, ref, alt)
#         # Return errors and passed variants
#         return((error, variants))

# # Check that either both or neither start and ends are present
#         for begin, end, (ref, alt) in alleles:
#             # Check current allele does not overlap previous
#             if begin < current_position:
#                 error = 'overlapping'
#                 variants.clear()
#                 break
#             else:
#                 current_position = end
#             # Get start position of allele within the read
#             try:
#                 read_begin = positions.index(begin)
#             except ValueError:
#                 try:
#                     read_begin = positions.index(begin - 1) + 1
#                 except ValueError:
#                     read_begin = None
#             # Get end position of allele within the read
#             try:
#                 read_end = positions.index(end)
#             except ValueError:
#                 try:
#                     read_end = positions.index(end - 1) + 1
#                 except ValueError:
#                     read_end = None