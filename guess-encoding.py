"""
guess-encoding.py: (c) Brent Pedersen

usage: awk 'NR % 4 == 0' your.fastq | python %prog [options]

guess the encoding of a stream of qual lines.
"""
import sys
import optparse

RANGES = {
    'Sanger, Phred+33': (33, 93),
    'Solexa, Solexa+64': (59, 104),
    'Illumina-1.3, Phred+64': (64, 104),
    'Illumina-1.5, Phred+64': (67, 104)
}

BOWTIE_PARAMS = {
    "Sanger, Phred+33": '--phred33-quals',
    "Solexa, Solexa+64": '--solexa-quals',
    "Illumina-1.3, Phred+64": '--phred64-quals',
    "Illumina-1.5, Phred+54": '--phred64-quals'
}

def get_qual_range(qual_str):
    """
    >>> get_qual_range("DLXYXXRXWYYTPMLUUQWTXTRSXSWMDMTRNDNSMJFJFFRMV")
    (68, 89)
    """

    vals = [ord(c) for c in qual_str]
    return min(vals), max(vals)

def get_bowtie_option(encoding):
    return BOWTIE_PARAMS[encoding]

def get_encodings_in_range(rmin, rmax, ranges=RANGES):
    valid_encodings = []
    for encoding, (emin, emax) in ranges.items():
        if rmin >= emin and rmax <= emax:
            valid_encodings.append(encoding)
    return valid_encodings

def main():
    p = optparse.OptionParser(__doc__)
    p.add_option("-n", dest="n", help="number of qual lines to test default:-1"
                 " means test until end of file or until it it possible to "
                 " determine a single file-type",
                 type='int', default=-1)
    p.add_option("-b", dest="b", help="make output the option to use in bowtie"
                 " corresponding to the quality type identified",
                 action="store_true", default=False)

    opts, args = p.parse_args()
    ##print >>sys.stderr, "# reading qualities from stdin"
    gmin, gmax  = 99, 0
    valid = []
    for i, line in enumerate(sys.stdin):
        lmin, lmax = get_qual_range(line.rstrip())
        if lmin < gmin or lmax > gmax:
            gmin, gmax = min(lmin, gmin), max(lmax, gmax)
            valid = get_encodings_in_range(gmin, gmax)
            if len(valid) == 0:
                ##print >>sys.stderr, "no encodings for range: %s" % str((gmin, gmax))
                sys.exit()
            if len(valid) == 1 and opts.n == -1:
                if (opts.b):
                    print get_bowtie_option(valid[0])
                else: 
                    print "\t".join(valid) + "\t" + str((gmin, gmax))
                sys.exit()

        if opts.n > 0 and i > opts.n:
            if (opts.b):
                print get_bowtie_option(valid[0])
            else: 
                print "\t".join(valid) + "\t" + str((gmin, gmax))
            sys.exit()

    if (opts.b):
        print get_bowtie_option(valid[0])
    else: 
        print "\t".join(valid) + "\t" + str((gmin, gmax))


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
