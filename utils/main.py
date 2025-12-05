import argparse
import os
from utils import jac_calculator, plot_jac
from utils.message import Message
from utils.__version__ import __version__


def get_opts():
    group = argparse.ArgumentParser()
    group.add_argument("-g", "--genome", help="Input genome file", required=True)
    group.add_argument(
        "-l",
        "--list",
        help="Haplotype group file, contain two columns: Chromosome\tHaplotype",
        required=True,
    )
    group.add_argument(
        "-w",
        "--window",
        help="Window size, could be scientific notation, default=1e6",
        type=float,
        default=1e6,
    )
    group.add_argument(
        "-s",
        "--step",
        help="Step size, could be scientific notation, default=5e5",
        type=float,
        default=5e5,
    )
    group.add_argument("-k", help="k size of kmer, default=21", type=int, default=21)
    group.add_argument(
        "--method",
        help="Method for calculating jaccard similarity, could be exact or sample, default=exact",
        choices=["exact", "sample"],
        default="exact",
    )
    group.add_argument(
        "--sample_ratio",
        help="Sample ratio for calculating, default=0.1",
        type=float,
        default=0.1,
    )
    group.add_argument(
        "--sample_seed", help="Random seed, default=None", type=int, default=None
    )
    group.add_argument("-o", "--output", help="Output directory", required=True)
    group.add_argument(
        "--cmap", help='CMAP for drawing heatmap, default="RdBu_r"', default="RdBu_r"
    )
    group.add_argument(
        "--fmt", help='Heatmap file format, default="pdf"', default="pdf"
    )
    group.add_argument(
        "--log_scale",
        help="If set, heatmap would be scaled with log10(jaccard similarity * 100+1)",
        action="store_true",
    )
    group.add_argument(
        "-t", "--threads", help="Threads, default=10", type=int, default=10
    )
    group.add_argument(
        "--verbose", help="Print detail information", action="store_true"
    )
    group.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )
    return group.parse_args()


def main():
    opts = get_opts()
    genome = opts.genome
    group_list = opts.list
    wsize = int(opts.window)
    ssize = int(opts.step)
    k = opts.k
    calc_method = opts.method
    sample_ratio = opts.sample_ratio
    if opts.sample_seed:
        sample_seed = opts.sample_seed
    else:
        sample_seed = None
    out_dir = opts.output
    cmap = opts.cmap
    pic_fmt = opts.fmt
    log_scale = opts.log_scale
    threads = opts.threads
    verbose = opts.verbose

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    Message.info("Calculating Jaccard similarity")
    out_jac = os.path.join(out_dir, "final.jac")
    if os.path.exists(out_jac):
        Message.info("Jac file found, skipping...")
    else:
        jac_calculator.win_kmer_jac_similarity(
            genome,
            group_list,
            wsize,
            ssize,
            k,
            sample_ratio,
            sample_seed,
            out_dir,
            threads,
            calc_method,
            verbose,
        )

    out_pic_dir = os.path.join(out_dir, "pic")
    if not os.path.exists(out_pic_dir):
        os.makedirs(out_pic_dir)

    Message.info("Plotting heatmap for each chromosome")
    plot_jac.plot_jac(out_jac, group_list, out_pic_dir, cmap, pic_fmt, log_scale)
    Message.info("Finished")
