"""
Script based post_stidsim.pl from
simulation-scripts.zip downloaded from here:
https://sites.google.com/eng.ucsd.edu/datasets/astral-ii

Modified from the script written by EKM (molloy.erin.k@gmail.com) in Spring 2018.
"""
import argparse
import numpy
import sys


def build_freqs(freq, size):
    nfrq = len(freq)
    fmat = numpy.zeros((size, nfrq))
    for j in range(nfrq):
            fmat[:, j] = numpy.random.gamma(freq[j], 1.0, size)
    rowsum = numpy.sum(fmat, axis=1)
    for j in range(nfrq):
        fmat[:, j] = fmat[:, j] / rowsum
    return fmat


def build_rate_matrix(rate, size):
    nrat = len(rate)
    rmat = numpy.zeros((size, nrat))
    for j in range(nrat):
            rmat[:, j] = numpy.random.gamma(rate[j], 1.0, size)
    rowsum = numpy.sum(rmat, axis=1)
    for j in range(nrat):
        rmat[:, j] = rmat[:, j] / rowsum
    for j in range(nrat):
        rmat[:, j] = rmat[:, j] / rmat[:, 1]  # Divide by AG for indelible
    return rmat


def build_alpha_vector(alpha, size):
    avec = numpy.log(numpy.random.uniform(low=0.0,
                                          high=1.0,
                                          size=size))
    avec = -1.0 * avec / alpha
    if numpy.where(avec == numpy.inf)[0].shape[0] != 0:
        sys.exit("Alpha was set to 0!")
    return avec


def build_length_vector(size, low=1000, high=1000):
    mus = numpy.random.uniform(low=numpy.log(low),
                               high=numpy.log(high),
                               size=size)
    sigmas = numpy.random.uniform(low=0.0,
                                  high=0.3,
                                  size=size)
    seqlen = []
    for mu, sigma in zip(mus, sigmas):
        seqlen.append(int(numpy.random.lognormal(mean=mu,
                                                 sigma=sigma)))
    return seqlen


def set_indelible_params(output, ntype, size):
    """
    """
    dtype = ["default"]
    freqs = [[36, 26, 28, 32]]
    rates = [[16, 3, 5, 5, 6, 15]]
    alphs = [1.199133]

    with open(output, 'w') as f:
        gene = 1
        f.write("GENE,TYPE,ALPH,SQLN,P_A,P_C,P_G,P_T,P_CA,P_AG,P_TA,P_CG,P_TC,P_TG\n")
        for i in range(len(dtype)):
            data = dtype[i]
            print(size)
            fmat = build_freqs(freqs[i], size)
            rmat = build_rate_matrix(rates[i], size)
            avec = build_alpha_vector(alphs[i], size)
            lvec = [ntype[i]]*size
            for j in range(size):
                fj = fmat[j, :]
                rj = rmat[j, :]
                f.write("%d,%s,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n"
                        % (gene, data, avec[j], lvec[j],
                           fj[0], fj[1], fj[2], fj[3],
                           rj[0], rj[1], rj[2], rj[3], rj[4], rj[5]))
                gene = gene + 1


def main(args):
    set_indelible_params(args.output, [args.seqlen], args.genenum)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--seqlen", type=int,
                        help="Sequence length", required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output file", required=True)
    parser.add_argument("-n", "--genenum", type=int,
                        help="Number of genes", required=True)

    main(parser.parse_args())
