"""
Runs INDELible using the file from set_indelible_params.py

Written by EKM (molloy.erin.k@gmail.com) in Spring 2018.
"""
import argparse
import numpy
import pandas
import os
import sys


def run_indelible(trees, params, outdir, tmpdir):
    os.makedirs(tmpdir)
    os.chdir(tmpdir)

    #ipath, iexec = indelible.rsplit('/', 1)
    #print(ipath, iexec)
    #print(os.getcwd())
    #os.system("cp " + indelible + " ./")

    pad = len(str(len(trees)))
    for i, r in params.iterrows():
        gene = str(r["GENE"]).zfill(pad)
        with open("control.txt", "w") as f:
            f.write("[TYPE] NUCLEOTIDE 1\n")
            f.write("[MODEL] modelname\n")
            f.write("[submodel] GTR %f %f %f %f %f\n"
                    % (r["P_TC"], r["P_TA"],
                       r["P_TG"], r["P_CA"], r["P_CG"]))  # TC TA TG CA CG AG
            f.write("[statefreq] %f %f %f %f\n"
                    % (r["P_T"], r["P_C"],
                       r["P_A"], r["P_G"]))  # T C A G
            f.write("[rates] 0 %f 0 \n" % r["ALPH"])
            f.write("[TREE] treename  " + trees[i] + "\n")
            f.write("[PARTITIONS] partitionname\n")
            f.write("[treename modelname %d]\n" % r["SQLN"])
            f.write("[EVOLVE] partitionname 1 %s\n" % gene)
        #os.system("./" + iexec)
        os.system("indelible")
        os.rename(gene + "_TRUE.phy", gene + ".phy")
        os.remove(gene + ".fas")

    #os.remove(iexec)
    os.remove("trees.txt")
    os.remove("control.txt")
    os.remove("LOG.txt")
    os.system("mv " + tmpdir + "/*.phy " + outdir)
    os.system("rm -rf " + tmpdir)
    #os.system("cat *phy | sed '/^ *$/d' > all-genes.phylip")


def main(args):
    with open(args.trees, 'r') as f:
        trees = [l for l in f]

    params = pandas.read_csv(args.params)
    params = params[(params["GENE"] >= args.start) &
                    (params["GENE"] <= args.end)]

    tmpdir = args.output + "/" + "tmp-" + str(args.start) + "-" + str(args.end)

    if not os.path.exists(tmpdir):
        run_indelible(trees, params,
                      args.output, tmpdir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    #parser.add_argument("-x", "--indelible", type=str,
    #                    help="Path to indelible executable", required=True)
    parser.add_argument("-s", "--start", type=int,
                        help="Start index", required=True)
    parser.add_argument("-e", "--end", type=int,
                        help="End index", required=True)
    parser.add_argument("-t", "--trees", type=str,
                        help="Tree list file", required=True)
    parser.add_argument("-p", "--params", type=str,
                        help="Parameter list file", required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output directory", required=True)

    main(parser.parse_args())
