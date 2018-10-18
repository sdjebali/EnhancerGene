# Calcul the Pearson correlation coefficient between the promoter and all DHS peak
# within less than 500kb.
# Usage: python Prom-other_correlation_eco.py DNase_DHS.bed
# DHS_promoter.bed > prom-other_correlation.txt

import argparse as ap
from collections import defaultdict
from collections import deque
from collections import namedtuple

from random import randint
from math import log10
import scipy.stats as stats
import numpy as np

# PEP8 suggère Camel case : nom de classe.
# http://sametmax.com/le-pep8-en-resume/


class Dhs(object):  # For each DHS.
    """
    Define each DHS peak as an object and compute firsts calculs: centered
    logcounts and standard deviation.
    """

    def __init__(self, chrom, start, end, identifier, counts):
        self._chrom = chrom
        self._start = start
        self._end = end
        self._id = identifier
        self.counts = counts  # All counts should be floats.
        self.sigma = np.std(self.log_counts)
        self.moy = np.mean(self.log_counts)
        self.log_counts_prime = np.array(
            [i - self.moy for i in self.log_counts])

    def __str__(self):
        return ("%s\t%s\t%s\t%s" % (self.chrom, self.start, self.end, self.id))

    def __repr__(self):
        return "%s:%s-%s\t%s\t%s" % (self._chrom, self._start, self._end,
                                     self._id, floats_to_str(self.counts))

    @property
    def chrom(self):
        return self._chrom

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def id(self):
        return self._id

    @property
    def log_counts(self):
        """
        Avoid false correlation and NA values
        Return: vector with the log counts.
        """
        epsilon = 0.001
        return [log10(i + epsilon) for i in self.counts]

    def is_upstream_from(self, other, dist=500000):
        """
        Returns True if self is 500kb downstream to the other Dhs
        """
        if self.chrom < other.chrom:
            return True
        elif (self.chrom == other.chrom and
              self.end < other.start - dist):
            return True
        else:
            return False

    def is_downstream_from(self, other, dist=500000):
        """
        Returns True if self is 500kb downstream to the other Dhs
        """
        if self.chrom < other.chrom:
            return True
        elif (self.chrom == other.chrom and
              self.start > other.end + dist):
            return True
        else:
            return False

    def is_a_neighbour_of(self, other, dist=500000):
        if (self.is_upstream_from(other, dist=dist) or
                self.is_downstream_from(other, dist=dist)):
            return False
        else:
            return True


def floats_to_str(counts):
    return "\t".join(["%0.2f" % x for x in counts])


def compute_pearsonr(a, b):
    return stats.pearsonr(a.log_counts, b.log_counts)[0]


# Again see http://sametmax.com/le-pep8-en-resume/
def dhs_correlation(dhs1, dhs2):
    """
    Compute the pearson correlation coefficient by hand thanks to calculs made
    in the class DHS.
    Params: take two DHS objects containing all the precomputed values.
    Return: pearson correlation coefficient.
    """
    numerator = np.dot(dhs1.log_counts_prime, dhs2.log_counts_prime)
    n = len(dhs1.log_counts_prime)
    coefP = numerator / (n * dhs1.sigma * dhs2.sigma)
    return str(coefP)   # Quelle horreur !!


def read_dhs_file(dhs_file):
    """
    Simply read the DHS file into an array of Dhs objectsfields
    """
    all_dhs = []  # Vector with all the object DHS from file_dnase.
    with open(dhs_file, "r") as f:
        next(f)   # skip first line, works only if file has always a header !!
        for line in f:
            fields = line.rstrip().split("\t")
            chrom, start, end, identifier = fields[0:4]
            start = int(start)
            end = int(end)
            counts = [float(c) for c in fields[4:]]
            all_dhs.append(Dhs(chrom, start, end, identifier, counts))
    return all_dhs


def write_correlations(correlations):
    header = "\t".join(["#prom_chr", "prom_start", "prom_end", "prom_ID",
                        "chr", "start", "end", "ID",
                        "Pearson_corr"])
    print(header)
    for pair, corr in correlations.items():
        prom, dhs = pair
        print("%s\t%s\t%3.2f" % (prom, dhs, corr))


def simulate_dhs(num_peaks=500, num_chrom=1,
                 num_tissues=10,
                 size_chrom=10000000):
    """
    Simulate num_peaks on num_chrom each of size size_chrom on num_tissues
    """
    all_dhs = []
    num_peaks_by_chrom = int(num_peaks/num_chrom)
    for chrom in range(1, num_chrom+1):
        for peak in range(0, num_peaks_by_chrom):
            start = randint(1, size_chrom)
            end = start + 500
            peak_id = "dnase_chrom%d_%d" % (chrom, peak)
            counts = 50 * np.random.random_sample(num_tissues)
            all_dhs.append(Dhs(chrom, start, end, peak_id, counts.tolist()))
    return all_dhs


# We define a named tuple
CorrelationResults = namedtuple('CorrelationResults',
                                ['corr', 'to_compute',
                                 'sorted_dhs', 'sorted_prom'])


def compute_correlations(all_dhs, promoters, dist=500000, compute=True):
    """
    Compute the correlations for each pair dhs x promoter such that
    dist(dhs, promoter) <= dist
    """
    # make sure the DHS are sorted according to chrom and pos
    all_dhs = sorted(all_dhs, key=lambda x: (x.chrom, x.start))
    promoters = sorted(promoters, key=lambda x: (x.chrom, x.start))

    correlations = defaultdict()
    to_compute = []
    # I would have liked to use a FIFO instead but FIFO doesn't allow
    # to check the first element, but we use it almost like a FIFO,
    # more precisely like an Abstract data type with the following operations
    #  - append to the right
    #  - pop from the left
    #  - check the element at the left of the queue (left formost ?)
    dq = deque()  # This queue will keep all dhs that have to be considered
    dhs = all_dhs.pop(0)
    dhs_counts = []
    for i, prom in enumerate(promoters):
        # First we remove all upstream dhs from dq, if any
        while dq and dq[0].is_upstream_from(prom, dist=dist):
            dq.popleft()
        # Second we skip all upstream dhs, if any
        while dhs.is_upstream_from(prom, dist=dist) and all_dhs:
            dhs = all_dhs.pop(0)
        # Now we insert in dq all dhs within the neighbourhood of prom
        while dhs.is_a_neighbour_of(prom, dist=dist) and all_dhs:
            dhs = all_dhs.pop(0)
            dq.append(dhs)
        # Now all the promoters in dq are within dist from dhs
        start = len(dhs_counts)
        for d in dq:
            if compute:
                r = compute_pearsonr(prom, d)
                correlations[(prom, dhs)] = r
                # print("corr( %s, %s ) = %3.2f" % (prom, d, r))
            else:
                dhs_counts.append(d.log_counts)
        end = len(dhs_counts)
        if not compute:
            to_compute.append((i, start, end))

    results = CorrelationResults(corr=correlations,
                                 to_compute=to_compute,
                                 sorted_dhs=all_dhs,
                                 sorted_prom=promoters
                                 )
    return results


if __name__ == '__main__':

    parser = ap.ArgumentParser()
    parser.add_argument("file_dnase")
    parser.add_argument("file_prom")

    args = parser.parse_args()
    file_dnase = args.file_dnase
    file_promoters = args.file_prom

    all_dhs = read_dhs_file(file_dnase)
    promoters = read_dhs_file(file_promoters)

    correlations = compute_correlations(all_dhs,
                                        promoters,
                                        dist=100000,
                                        compute=True)

    write_correlations(correlations.corr)
