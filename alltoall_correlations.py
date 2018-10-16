#Calcul the Pearson correlation coefficient between the promoter and
#all DHS peak within less than 500kb.
#Usage:
#python alltoall_correlations.py [-p prom_file.bed] [-d 500000] DNase_DHS.bed
#CPU = 8 by default

import argparse as ap
from math import log10
import numpy as np
import sys
from multiprocessing import Pool, sharedctypes
import warnings
import scipy.stats as ss


def eprint(*args, **kwargs):
    ''' A simple to sys.stderr printer wrapper '''
    print(*args, file=sys.stderr, **kwargs)


class Dhs(object): #For each DHS.
    """
    Define each DHS peak as an object and compute firsts calculs: centered
    logcounts and standard deviation.
    """
    def __init__(self, chrom, start, end, name, counts):
        self._chrom = chrom
        self._start = int(start)
        self._end = int(end)
        self._name = name
        self.counts = counts
        self.sigma = np.std(self.log_counts)
        self.moy = np.mean(self.log_counts)
        self.log_counts_prime = np.array([i-self.moy for i in self.log_counts])

    def __str__(self):
        return ("%s\t%s\t%s\t%s" % (self.chrom,self.start,self.end,self.name))

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
    def name(self):
        return self._name

    @property
    def length(self):
        """
        Return the length of the counts vector.
        """
        return (len(self.counts))

    @property
    def log_counts(self):
        """
        Avoid false correlation and NA values
        Return: vector with the log counts.
        """
        epsilon = 0.001
        return [log10(i+epsilon) for i in self.counts]

    @property
    def rangs(self):
        return ss.rankdata(self.counts)

    @property
    def interval(self):
        """
        Create a vector about the position of the DHS peak.
        Return: vector [chrom,start,end]
        """
        return [self.chrom,self.start,self.end]

    def within_neighbourhood(self,dhs_interval,distance, calcul_distance):
        """
        Return True if DHS interval is within 500kb of self (the promoter)
        Param dhs_interval: vector [chrom,start,end] and distance (interger)
        Return: bolean value.
        """
        if calcul_distance:
            middle1 = self.start + (self.end-self.start)/2
            middle2 = dhs_interval[1] + (dhs_interval[2]-dhs_interval[1])/2
            if self.chrom == dhs_interval[0] and (self.end < dhs_interval[1] and middle1+distance >= middle2):
                return True
            else:
                return False
        else:
            if self.chrom == dhs_interval[0] and (self.end < dhs_interval[1] and self.end+distance >= dhs_interval[1]):
                return True
            else:
                return False


def read_dhs_file(dhs_file, chromosome="all"):
    """
    Simply read the DHS file into an array of Dhs objectsfields
    """
    all_dhs = []  # Vector with all the object DHS from file_dnase.
    with open(dhs_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            chrom, start, end, identifier = fields[0:4]
            if chromosome != "all" and chrom != chromosome:
                continue
            start = int(start)
            end = int(end)
            counts = [float(c) for c in fields[4:]]
            all_dhs.append(Dhs(chrom, start, end, identifier, counts))

    return all_dhs


def write_correlations(output, all_dhs):

    header = "\t".join(["#chr", "start", "end", "ID",
                        "chr", "start", "end", "ID",
                        "corr"])
    print(header)
    for i, j, corr in output:
        dhs1 = all_dhs[i]
        dhs2 = all_dhs[j]
        print("%s\t%s\t%5.4f" % (dhs1, dhs2, corr))


def arraytonumpy(all_dhs, method="pearson"):
    """
    Convert values from the class Dhs (logcounts and std) to numpy arrays for
    multiprocessing.
    Param: the table with all the object.
    Return: 3 numpy arrays : logcounts, std, ranks
    """
    logcounts = []
    stdeviation = []
    ranks = []
    for i in all_dhs:
        logcounts.append(i.log_counts_prime)
        stdeviation.append(i.sigma)
        if method == "spearman":
            ranks.append(i.rangs)
    return (np.array(logcounts), np.array(stdeviation), np.array(ranks))


def calcul_correlation(logcount1, logcount2, sigma1, sigma2, rangs, bis):
    """
    Compute the pearson correlation coefficient by hand thanks to calculs made
    in the class Dhs.
    Params: take two Dhs objects containing all the precomputed values.
    Return: pearson correlation coefficient.
    """
    if len(rangs) == 0:
        numerator = np.dot(logcount1, logcount2)
        coefP = numerator/(10*sigma1*sigma2)
        return coefP
    else:
        numerator = np.dot(logcount1, logcount2)
        coefP = numerator/(10*sigma1*sigma2)
        diff = 0
        for i in range(len(rangs)):
            diff += (rangs[i] - bis[i])**2
        coefS = 1 - (6*diff)/(10*((10**2)-1))
        return coefP, coefS


def compute_correlation(logcounts, sigma, rank, i, j, mode="pearson"):
        if mode == "pearson":
            return pearsonr(logcounts[i], logcounts[j],
                            sigma[i], sigma[j])
        elif mode == "spearman":
            return spearmanr(rank[i], rank[j])
        else:
            eprint("Unknown correlation method")
            exit(1)


def pearsonr(logcount_i, logcount_j, sigma_i, sigma_j):
    numerator = np.dot(logcount_i, logcount_j)
    coefP = numerator/(logcount_i.size*sigma_i*sigma_j)
    return coefP


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
    start, end, mode = args

    # Warning is necessary because of the use of global variables
    # in parallelization mode
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        logcounts = np.ctypeslib.as_array(logcounts_global)
        sigma = np.ctypeslib.as_array(sigma_global)
        ranks = np.ctypeslib.as_array(ranks_global)

    correlations = []
    i = start
    for j in range(start+1, end+1):
        r = compute_correlation(logcounts, sigma, ranks,  i, j, mode)
        correlations.append(((i, j), r))

    return correlations


def _init_parallel(logcounts_to_populate, sigma_to_populate,
                   ranks_to_populate):
    """
    Each pool process calls this initializer.
    Load the arrays to be populated into that process's global namespace
    """
    global logcounts_global
    global sigma_global
    global ranks_global
    logcounts_global = logcounts_to_populate
    sigma_global = sigma_to_populate
    ranks_global = ranks_to_populate


def parallel_implementation_corr(logcounts, sigma, ranks, to_compute,
                                 threads=1):
    """
    Create shared variables and pool for multiprocessing.
    Params: logcounts and sigma numpy array, to_compute list with the index and
    the number of processeurs to use.
    Return: all the correlations and index
    """
    # Define the corresponding shared ctype arrays
    # Create ctypes object from numpy array
    tmp_logcounts = np.ctypeslib.as_ctypes(logcounts)
    shared_logcounts = sharedctypes.Array(tmp_logcounts._type_, tmp_logcounts,
                                          lock=False)
    #     Returns a ctypes array allocated from shared memory.
    #     tmp_dhs._type_: type of the returned array elements
    #     tmp_dhs: initialize the array
    #     lock=True (default): synchronize access to the value
    #     lock=False: access to returned object not protected

    tmp_sigma = np.ctypeslib.as_ctypes(sigma)
    shared_sigma = sharedctypes.Array(tmp_sigma._type_, tmp_sigma, lock=False)

    tmp_ranks = np.ctypeslib.as_ctypes(ranks)
    shared_ranks = sharedctypes.Array(tmp_ranks._type_, tmp_ranks,
                                          lock=False)

    pool = Pool(processes=threads, initializer=_init_parallel,
                initargs=(shared_logcounts, shared_sigma, shared_ranks, ))
    #                       controls a pool of worker processes
    #                       processes: number of worker processes to use
    #                       initializer is not None: each worker process will
    #                       call initializer(*initargs)
    result = pool.map(parallel_correlations, to_compute)
    #                       supports only one iterable argument,
    #                       blocks until the result is ready

    results = [item for sublist in result for item in sublist]
    return results


def sort_correlations(correlations, all_dhs):
    """
    Sort the result from multiprocessing and add the data corresponding to each
    coefficient.
    Params: list with index and correlation coefficient and list of all the
    object Dhs
    Result: sorted numpy array with the results to print.
    """
    results = []
    for cor in sorted(correlations,  key=lambda x: (x[0][0], x[0][1])):
        pair, r = cor
        i, j = pair
        results.append((i, j, r))
    return results


def get_pairs_tocompare(all_dhs, distance, calcul_distance, method="pearson"):
    """
    Identify the paris verifying the distance constraints hence the pairs
    for whixh a correlation will be computed
    """
    to_compute = []
    for i in range(len(all_dhs)-1):
        right_border = i+1
        while all_dhs[i].within_neighbourhood(all_dhs[right_border].interval,
                                              distance, calcul_distance):
            right_border += 1
            if right_border == len(all_dhs):
                break
        # For each peak, we compute correlation with all
        # the peaks in this window.
        to_compute.append([i, right_border-1, method])
    return to_compute


def compute_all_correlations(all_dhs, to_compute, calcul_distance, distance=500000,
                             method="pearson", threads=1):
    #to_compute = get_pairs_tocompare(all_dhs, distance, calcul_distance, method)
    (logcounts, sigma, ranks) = arraytonumpy(all_dhs, method=method)
    output = parallel_implementation_corr(logcounts, sigma, ranks,
                                          to_compute, threads=threads)
    return output


if __name__ == '__main__':
    # Input from command line.
    parser = ap.ArgumentParser(description="Compute Dnase correlations from a "
                                            "file in bed-like format: \n "
                                            " chr start end name val1 val2 ...\n"
                                            "\nWarning, the file has to be sorted (using sort -k1,1 -k2,2n)",
                               formatter_class=ap.RawTextHelpFormatter)
    parser.add_argument("file_dnase", help="file with the read count per peak")
    parser.add_argument("-c", "--chrom", help="if you want the calcul only on \n"
                        "one chromosome.", default="all", type=str)
    parser.add_argument("-d", "--distance", help="distance maximum between peaks \n"
                        "to compute the correlation. (default:500000)",
                        default=500000, type=int)
    parser.add_argument("-l", "--length", help="how to calcul the distance between \n"
                        "two peaks, if actif distance compute between middle of \n"
                        "each peaks. (default: start-end)", action="store_true") #nom a revoir
    parser.add_argument("-m", "--method", help="correlation ncoefficient,\n"
                        "pearson or spearman. (default:pearson)",
                        default="pearson", type=str)
    parser.add_argument("-t", "--threads", help="number of prcessors\
                         (default:1)",
                        default=1, type=int)
    args = parser.parse_args()

    dnase_file = args.file_dnase
    distance = args.distance
    threads = args.threads
    chrom = args.chrom
    method = args.method
    calcul_distance = args.length

    eprint("Reading dhs file")
    all_dhs = read_dhs_file(dnase_file, chrom)
    eprint("%d peaks loaded" % (len(all_dhs)))

    eprint("Identifying correlation pairs to compute")
    to_compute = get_pairs_tocompare(all_dhs, distance, calcul_distance)

    eprint("Computing the correlations")
    output = compute_all_correlations(all_dhs, to_compute, calcul_distance, distance=distance,
                                      method=method, threads=threads)

    eprint("Sorting the correlations")
    sorted_correlations = sort_correlations(output, all_dhs)

    eprint("Writing the correlations")
    write_correlations(sorted_correlations, all_dhs)
