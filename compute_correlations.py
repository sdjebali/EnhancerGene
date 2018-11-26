# Compute pearson and spearman correlation coefficient between
# feature signals (Dnase, ...)
# Usage:
# python compute_correlations.py [-d 500000] [-c 10] feature_signals.bdelike
#

import argparse as ap
from math import log10
import numpy as np
import sys
import os
import tempfile
import itertools
from multiprocessing import Pool, sharedctypes
import warnings
import scipy.stats as ss


def eprint(message):
    ''' A simple to sys.stderr printer wrapper '''
    print(message, file=sys.stderr)


def floats_to_str(values):
    return "\t".join(["%3.2f" % x for x in values])


class GenomeFeature(object):
    """
    An object storing a genome feature (an annotation): chrom, start, end
    and a name
    """
    def __init__(self, chrom, start, end, name=None):
        self._chrom = chrom
        self._start = start
        self._end = end
        self._name = name

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
    def name(self):
        return self._name

    @property
    def center(self):
        return (self.start + self.end - 1)/2


class FeatureSignal(GenomeFeature):
    """
    FeatureSignal signal object store a functional element (feature) and
    associated signal in different tissues or individuals
    For each feature, mean and standard deviations are computed at
    initialization and the associated standardized value
    The standard deviation can be computed for logs or for raw
    logcounts and standard deviation.

    group enables to sore in a single array signals from different sources
    """
    def __init__(self, chrom, start, end, name, values, group="A", log=False):
        super(FeatureSignal, self).__init__(chrom, start, end, name)
        self._group = group
        self.values = values
        if log:
            self.raw_x = self.log_values
        else:
            self.raw_x = self.values
        self.mean = np.mean(self.raw_x)
        self.sigma = np.std(self.raw_x, ddof=1)
        self._std_values = np.array([i-self.mean for i in self.raw_x])/self.sigma

    def __str__(self):
        return ("%s\t%s\t%s\t%s" % (self.chrom, self.start, self.end,
                                    self.name))

    def __repr__(self):
        if self._std_values is None:
            values = [float('nan')] * self.length
        else:
            values = self._std_values
        return "%s:%s-%s\t%s\t%3.2f\t%3.2f\t%s" % (self.chrom, self.start,
                                                   self.end, self.name,
                                                   self.sigma, self.mean,
                                                   floats_to_str(values))

    @property
    def group(self):
        return self._group

    @property
    def std(self):
        return self.sigma

    @property
    def std_values(self):
        return self._std_values

    @property
    def length(self):
        """
        Return the length of the counts vector.
        """
        return (len(self.values))

    @property
    def log_values(self):
        """
        Avoid false correlation and NA values
        Return: vector with the log counts.
        """
        epsilon = 0.001
        return [log10(x+epsilon) for x in self.values]

    @property
    def rangs(self):
        return ss.rankdata(self.values)

    @property
    def interval(self):
        """
        Create a vector about the position of the DHS peak.
        Return: vector [chrom,start,end]
        """
        return [self.chrom, self.start, self.end]

    def within_neighbourhood(self, right_candidate_fs, max_distance):
        """
        Return True if right_candidate_fs interval is within 500kb of self
        Param max_distance: the (half) size of the neighborhood (interger)
        Return: bolean value.
        """
        pos1 = self.center
        pos2 = right_candidate_fs.center
        if (self.chrom == right_candidate_fs.chrom and
                self.start > right_candidate_fs.start):
            eprint("The file seems to be unsorted:")
            eprint(self)
            eprint(right_candidate_fs)
            eprint("Exiting")
            exit(1)
        if (self.chrom == right_candidate_fs.chrom and
                pos2 - pos1 <= max_distance):
            return True
        else:
            return False


def ValidEntry(chrom, start, end, values, num_signals):
    """
    Returns true if the genome feature signal is elligible for correlation
    computation:
      - start < end,
      - values are not all zeros
      - at least one value different from the orhters
      - required number of values
    """
    if start < 0 or start >= end:
        return False
    if np.count_nonzero(values) == 0:
        return False
    if np.std(values, ddof=1) == 0:
        return False
    if len(values) != num_signals:
        return False
    return True


def read_signal_file(signal_file, chromosome="all", log_transform=False,
                     group="A"):
    """
    Read the signal file into an array of FeatureSignal objects fields
    File format
      chrom start end identifier s1 s2 ... sk
    """
    all_signals = []  # Vector with all the object DHS from file_dnase.
    num_signals = 0
    with open(signal_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            chrom, start, end, identifier = fields[0:4]
            if num_signals == 0:
                num_signals = len(fields[4:])
            if chromosome != "all" and chrom != chromosome:
                continue
            start = int(start)
            end = int(end)
            values = [float(c) for c in fields[4:]]
            if not ValidEntry(chrom, start, end, values, num_signals):
                eprint("Skipping \n%s\nNot a valid or a useful entry" %
                       line.rstrip())
                continue
            all_signals.append(FeatureSignal(chrom, int(start), int(end),
                                             identifier, values,
                                             log=log_transform, group=group))

    return all_signals


# def write_header(output_fh):
#     """
#     Simply write the header using the output_fh handler
#     """
#     header_line = "\t".join(["#chr", "start", "end", "ID",
#                              "chr", "start", "end", "ID",
#                              "pearson_r", "spearman_r"]) + "\n"
#     output_fh.write(header_line)


def dump_correlations(correlations, output_fh):
    """
    Dump correlations with associated genome feature coordinates and name
    """
    # all_feature_signals is a global variable
    for pair, corr in correlations:
        i, j = pair
        pea_r, spea_r = corr
        fs1 = all_feature_signals[i]
        fs2 = all_feature_signals[j]
        output_fh.write("%s\t%s\t%5.4f\t%5.4f\n" % (fs1, fs2, pea_r, spea_r))


def sort_and_write(tmp_file, output_file):
    """
    Sort, using unix sort, the signa features and append to the output file
    """
    header_line = "\t".join(["#chr", "start", "end", "ID",
                             "chr", "start", "end", "ID",
                             "pearson_r", "spearman_r"]) + "\n"
    with open(output_file, "w") as fout:
        fout.write(header_line)
    cmd = "sort -k 1,1V -k 2,2n -k 5,5V -k 6,6n " + tmp_file + " >> " + output_file
    os.system(cmd)


def write_correlations(correlations, all_values):
    """
    Write correlations to stdout
    """

    header = "\t".join(["#chr", "start", "end", "ID",
                        "chr", "start", "end", "ID",
                        "pearson_r", "spearman_r"])
    print(header)
    for pair, corr in correlations:
        i, j = pair
        pea_r, spea_r = corr
        fs1 = all_values[i]
        fs2 = all_values[j]
        print("%s\t%s\t%5.4f\t%5.4f" % (fs1, fs2, pea_r, spea_r))


def arraytonumpy(all_values):
    """
    Convert values from the class FeatureSignal to numpy arrays for
    multiprocessing.
    Param: the table with all the object.
    Return: 2 numpy arrays : logcounts,ranks
    """
    values = []
    ranks = []
    for i in all_values:
        values.append(i.std_values)
        ranks.append(i.rangs)
    return (np.array(values), np.array(ranks))


def compute_correlation(values, rank, i, j):
    pearson_r = pearsonr(values[i], values[j])
    spearman_r = spearmanr(rank[i], rank[j])
    # pearson_ss = ss.pearsonr(values[i], values[j])
    return pearson_r, spearman_r


def pearsonr(values_i, values_j):
    return np.dot(values_i, values_j)/(len(values_i)-1)


def spearmanr(rank_i, rank_j):
    n = len(rank_i)
    d = np.sum([(rank_i[k] - rank_j[k])**2 for k in range(n)])
    return 1-6*d/(n*(n**2-1))


def parallel_correlations(args):
    """
    Compute all the correlation coefficient calcul of the Dhs in start from
    start+1 to end.
    Param: to_compute with the windows border index for each Dhs.
    Return: list of the correlation coefficient.
    """
    start, end = args

    # Warning is necessary because of the use of global variables
    # in parallelization mode
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        values = np.ctypeslib.as_array(values_global)
        ranks = np.ctypeslib.as_array(ranks_global)

    correlations = []
    i = start
    for j in range(start+1, end+1):
        pear_r, spea_r = compute_correlation(values, ranks,  i, j)
        correlations.append(((i, j), (pear_r, spea_r)))
    return correlations


def parallel_correlations_pairs(pairs):
    """
    Compute all the correlations of the Dhs in start from
    start+1 to end.
    Param: to_compute with the windows border index for each Dhs.
    Return: list of the correlation coefficient.
    """

    # Warning is necessary because of the use of global variables
    # in parallelization mode
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        values = np.ctypeslib.as_array(values_global)
        ranks = np.ctypeslib.as_array(ranks_global)

    correlations = []
    for p in pairs:
        i, j = p
        pear_r, spea_r = compute_correlation(values, ranks,  i, j)
        correlations.append(((i, j), (pear_r, spea_r)))
    return correlations


def _init_parallel(values_to_populate, ranks_to_populate):
    """
    Each pool process calls this initializer.
    Load the arrays to be populated into that process's global namespace
    """
    global values_global
    global ranks_global

    values_global = values_to_populate
    ranks_global = ranks_to_populate


def grouping(n, iterable):
    """
    Group the iterable by chunks of n
    """
    args = [iter(iterable)] * n
    return ([e for e in t if e is not None] for t in itertools.zip_longest(*args))


def parallel_computation(values, ranks, to_compute, output_file, threads=1):
    """
    Create shared variables and pool for multiprocessing.
    Params: logcounts and sigma numpy array, to_compute list with the index and
    the number of processeurs to use.
    Return: all the correlations and index
    """
    # Define the corresponding shared ctype arrays
    # Create ctypes object from numpy array
    tmp_values = np.ctypeslib.as_ctypes(values)
    shared_values = sharedctypes.Array(tmp_values._type_, tmp_values,
                                       lock=False)
    #     Returns a ctypes array allocated from shared memory.
    #     tmp_dhs._type_: type of the returned array elements
    #     tmp_dhs: initialize the array
    #     lock=True (default): synchronize access to the value
    #     lock=False: access to returned object not protected

    tmp_ranks = np.ctypeslib.as_ctypes(ranks)
    shared_ranks = sharedctypes.Array(tmp_ranks._type_, tmp_ranks,
                                      lock=False)
    # results = []
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    with open(temp_file.name, "w") as temp_fh:
        with Pool(processes=threads, initializer=_init_parallel,
                  initargs=(shared_values, shared_ranks, )) as pool:
            #     controls a pool of worker processes
            #     processes: number of worker processes to use
            #     initializer is not None: each worker process will
            #     call initializer(*initargs)

            #     we split the elements to compute in batchs of size 50000
            num_pass = 0
            for group_to_compute in grouping(50000, to_compute):
                num_pass += 1
                eprint("Pass %d" % num_pass)
                nested_result = pool.map(parallel_correlations, group_to_compute, 100)
                #     supports only one iterable argument,
                #                       blocks until the result is ready
                result = [item for sublist in nested_result for item in sublist]
                dump_correlations(result, temp_fh)

    sort_and_write(temp_file.name,  output_file)

    return 0


def sort_correlations(correlations):
    """
    Sort the result from multiprocessing and add the data corresponding to each
    coefficient.
    Params: the correlations
    Result: sorted numpy array with the results to print.
    """
    results = []
    for cor in sorted(correlations,  key=lambda x: (x[0][0], x[0][1])):
        pair, r = cor
        i, j = pair
        pear_r, spea_r = r
        results.append((i, j, pear_r, spea_r))
    return results


def _within_neighborhood(all_values_A, all_values_B=None,
                         max_distance=500000):
    """
    Identify the pairs verifying the distance constraints hence the pairs
    for whixh a correlation will be computed
    """
    all_values = all_values_A
    # merge all values in a single array
    if all_values_B is not None:
        all_values.extend(all_values_B)
    # sort all values in place
    all_values.sort(key=lambda x: (x.chrom, x.start))

    to_compute = []
    for i in range(len(all_values)-1):
        right_border = i+1
        while all_values[i].within_neighbourhood(all_values[right_border],
                                                 max_distance):
            right_border += 1
            if right_border == len(all_values):
                break
        # For each peak, we compute correlation with all
        # the peaks in this window.
        to_compute.append([i, right_border-1])
    return to_compute, all_values


def _peer2peer(all_values_A, all_values_B):
    if all_values_B is None:
        eprint("Impossible line by line comparison without B values")
        exit(1)
    if len(all_values_A) != len(all_values_B):
        eprint("Impossible line by line comparison when the two files "
               "differ in the number of valid lines")
        exit(1)
    all_values = []
    to_compute = []
    for i, a, b in zip(range(0, len(all_values_A)),
                       all_values_A, all_values_B):
        all_values.extend([a, b])
        to_compute.append([2*i, 2*i+1])
    return to_compute, all_values


def get_pairs_tocompare(all_values, all_values_B=None,
                        max_distance=500000, mode="all2all"):
    """
    Identify the pairs verifying the distance constraints hence the pairs
    for whixh a correlation will be computed
    """
    to_compute = []

    if mode == "all2all" or mode == "A2B":
        to_compute = _within_neighborhood(all_values, all_values_B,
                                          max_distance=max_distance)
    elif mode == "L2L":
        to_compute = _peer2peer(all_values, all_values_B)

    return to_compute


def compute_all_correlations(all_values, to_compute, output_file, threads=1):
    (values, ranks) = arraytonumpy(all_values)
    output = parallel_computation(values, ranks, to_compute, output_file,
                                  threads=threads)
    return output


if __name__ == '__main__':
    # Input from command line.
    parser = ap.ArgumentParser(description="Compute pearson and spearman "
                               "correlations between \n"
                               "feature signals. \n"
                               "Correlations are computed bewteen features "
                               "separated \n"
                               "by at most max_distance bp "
                               "(default 500kb)\n\n"
                               "The feature signal file mus be in the "
                               "following format:\n"
                               " chr start end name val1 val2...\n",
                               formatter_class=ap.RawTextHelpFormatter)
    parser.add_argument("signal_file", help="file with the feature signals")
    parser.add_argument("-o", "--output", help="output file", required=True)
    parser.add_argument("-c", "--chrom",
                        help="a single chromosome (default:all chromosomes)",
                        default="all", type=str)
    parser.add_argument("-d", "--distance", dest='max_distance',
                        help="maximum distance between features "
                        "for correlation computation (default:500000)",
                        default=500000, type=int)
    parser.add_argument("-l", "--log-transformed", dest='log_transform',
                        help="log transform the data (default:False)",
                        action="store_true")
    parser.add_argument("-m", "--mode",
                        help="mode of comparison all to all (all2all), "
                             "file a to file b (A2B), line pairs (L2L) "
                             "(default:all2all)",
                        default="all2all", type=str)

    parser.add_argument("-t", "--threads", help="number of prcessors"
                        " (default:1)",
                        default=1, type=int)
    args = parser.parse_args()

    signal_file = args.signal_file
    output_file = args.output
    chrom = args.chrom
    max_distance = args.max_distance
    log_transform = args.log_transform
    mode = args.mode
    threads = args.threads

    global all_feature_signals

    eprint("Reading signal file")
    all_feature_signals = read_signal_file(signal_file, chrom, log_transform)
    eprint("%d peaks loaded from file" % (len(all_feature_signals)))

    eprint("Identifying correlation pairs to compute")
    to_compute, all_feature_signals = get_pairs_tocompare(all_feature_signals,
                                                          max_distance=max_distance,
                                                          mode=mode)

    eprint("Computing %d correlation intervals" % len(to_compute))
    output = compute_all_correlations(all_feature_signals, to_compute,
                                      output_file,
                                      threads=threads)
