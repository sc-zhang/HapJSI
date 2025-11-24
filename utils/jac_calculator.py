import multiprocessing
import os
import random
from utils.message import Message


GENOME_DB = {}


def load_genome(in_genome, group_db):
    seq_db = {}
    with open(in_genome, "r") as fin:
        for line in fin:
            if line[0] == ">":
                sid = line.strip()[1:]
                seq_db[sid] = []
            else:
                seq_db[sid].append(line.strip().upper())
    fa_db = {}
    for sid in seq_db:
        if sid in group_db:
            fa_db[sid] = "".join(seq_db[sid])
    return fa_db


def load_group(group_file):
    group_db = {}
    with open(group_file, "r") as fin:
        for line in fin:
            data = line.strip().split()
            group_db[data[1]] = data[0]
    return group_db


def reverse_seq(seq):
    base_db = {"A": "T", "T": "A", "G": "C", "C": "G"}
    rev_seq = [base_db[_] if _ in base_db else _ for _ in seq[::-1]]
    return "".join(rev_seq)


def gen_kmer(seq, k, sample_ratio=1.0, sample_seed=None):
    if not sample_seed:
        random.seed()
    else:
        random.seed(sample_seed)
    kmer_set = set()
    for i in range(len(seq) - k + 1):
        if random.random() < sample_ratio:
            kmer_set.add(seq[i : i + k])
    return kmer_set


def calc_jaccard(hap1, sp1, ep1, hap2, sp2, ep2, k, sample_ratio, sample_seed, verbose):
    if verbose:
        Message.info(
            "\tComparing %s: %d-%d, %s: %d-%d"
            % (hap1, sp1 + 1, ep1, hap2, sp2 + 1, ep2)
        )
    kmer_set1 = gen_kmer(GENOME_DB.get(hap1)[sp1:ep1], k, sample_ratio, sample_seed)
    kmer_rev_set1 = gen_kmer(reverse_seq(GENOME_DB.get(hap1)[sp1:ep1]), k, sample_ratio, sample_seed)
    kmer_set2 = gen_kmer(GENOME_DB.get(hap2)[sp2:ep2], k, sample_ratio, sample_seed)

    try:
        jac = len(kmer_set1 & kmer_set2) * 1.0 / (len(kmer_set1 | kmer_set2))
    except ZeroDivisionError:
        jac = 0.0
    try:
        jac_rev = (
            len(kmer_rev_set1 & kmer_set2) * 1.0 / (len(kmer_rev_set1 | kmer_set2))
        )
    except ZeroDivisionError:
        jac_rev = 0.0

    # positive first
    return jac if jac >= jac_rev else -jac_rev


def win_kmer_jac_similarity(
    in_genome,
    group_file,
    wsize,
    ssize,
    k,
    sample_ratio,
    sample_seed,
    out_dir,
    threads,
    method="exact",
    verbose=False,
):
    Message.info("Loading group file")
    group_db = load_group(group_file)

    try:
        multiprocessing.set_start_method("fork")
    except RuntimeError:
        pass

    Message.info("Loading genome")
    global GENOME_DB
    GENOME_DB = load_genome(in_genome, group_db)

    Message.info("Pairwise comparing")
    pool = multiprocessing.Pool(processes=threads)
    compare_region_db = {}
    for hap in group_db:
        chrn = group_db[hap]
        if chrn not in compare_region_db:
            compare_region_db[chrn] = []
        for _ in range(0, len(GENOME_DB[hap]) - ssize + 1, ssize):
            compare_region_db[chrn].append(
                [hap, _, min(_ + wsize, len(GENOME_DB[hap]))]
            )

    res = []
    for chrn in compare_region_db:
        for i in range(len(compare_region_db[chrn]) - 1):
            hap1, sp1, ep1 = compare_region_db[chrn][i]
            for j in range(i + 1, len(compare_region_db[chrn])):
                hap2, sp2, ep2 = compare_region_db[chrn][j]
                r = pool.apply_async(
                    calc_jaccard,
                    (
                        hap1,
                        sp1,
                        ep1,
                        hap2,
                        sp2,
                        ep2,
                        k,
                        1.0 if method == "exact" else sample_ratio,
                        sample_seed,
                        verbose,
                    ),
                )
                res.append([chrn, hap1, sp1, hap2, sp2, r])

    pool.close()
    pool.join()

    Message.info("Writing Jac")
    out_jac = os.path.join(out_dir, "final.jac")
    with open(out_jac, "w") as fout:
        for chrn, hap1, sp1, hap2, sp2, r in sorted(res):
            try:
                jac = r.get()
                fout.write(
                    "%s\t%s\t%d\t%s\t%d\t%f\n" % (chrn, hap1, sp1, hap2, sp2, jac)
                )
            except Exception as e:
                print(
                    "Error found while comparing {} {} {} and {} {} with {}".format(
                        chrn, hap1, sp1, hap2, sp2, e
                    )
                )

    Message.info("Finished")
