import multiprocessing
import os
from utils.message import Message


def rev_kmer(kmer):
    base_db = {"A": "T", "T": "A", "G": "C", "C": "G"}
    rev_kmer = [base_db[_] if _ in base_db else _ for _ in kmer[::-1]]
    return "".join(rev_kmer)


def kmer_jac(kmer_file1, rev_kmer_file1, kmer_file2):
    kmer_set1 = set()
    kmer_rev_set1 = set()
    kmer_set2 = set()
    with open(kmer_file1, "r") as fin:
        for line in fin:
            kmer_set1.add(line.strip())

    with open(rev_kmer_file1, "r") as fin:
        for line in fin:
            kmer_rev_set1.add(line.strip())
    with open(kmer_file2, "r") as fin:
        for line in fin:
            kmer_set2.add(line.strip())

    jac = (
        len(kmer_set1.intersection(kmer_set2)) * 1.0 / (len(kmer_set1.union(kmer_set2)))
    )
    jac_rev = (
        len(kmer_rev_set1.intersection(kmer_set2))
        * 1.0
        / (len(kmer_rev_set1.union(kmer_set2)))
    )

    # positive first
    return jac if jac >= jac_rev else -jac_rev


def gen_kmer_file(sid, seq, wsize, ssize, k, out_dir):
    for _ in range(0, len(seq) - ssize + 1, ssize):
        sub_seq = seq[_ : _ + wsize]
        fn = "%s::%d::%d.for" % (sid, _ + 1, min(_ + wsize, len(seq)))
        full_fn = os.path.join(out_dir, fn)
        rev_fn = "%s::%d::%d.rev" % (sid, _ + 1, min(_ + wsize, len(seq)))
        rev_full_fn = os.path.join(out_dir, rev_fn)
        kmer_set = set()
        rev_kmer_set = set()
        for idx in range(len(sub_seq) - k + 1):
            kmer_set.add(sub_seq[idx : idx + k])
            rev_kmer_set.add(rev_kmer(sub_seq[idx : idx + k]))
        with open(full_fn, "w") as fout:
            fout.write("%s\n" % ("\n".join(sorted(kmer_set))))
        with open(rev_full_fn, "w") as fout:
            fout.write("%s\n" % ("\n".join(sorted(rev_kmer_set))))


def win_kmer_jac_similarity(in_genome, group_file, wsize, ssize, k, out_dir, threads):
    Message.info("Loading group file")
    group_db = {}
    with open(group_file, "r") as fin:
        for line in fin:
            data = line.strip().split()
            group_db[data[1]] = data[0]
    Message.info("Loading genome")
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

    Message.info("Creating kmers")
    kmer_dir = os.path.join(out_dir, "kmers")
    if not os.path.exists(kmer_dir):
        os.makedirs(kmer_dir)
        pool = multiprocessing.Pool(processes=threads)
        for sid in fa_db:
            pool.apply_async(
                gen_kmer_file,
                (
                    sid,
                    fa_db[sid],
                    wsize,
                    ssize,
                    k,
                    kmer_dir,
                ),
            )
        pool.close()
        pool.join()
    else:
        Message.info("Kmers found, skipping...")

    Message.info("Pairwise comparing")

    hap_kmer_file_db = {}

    for fn in os.listdir(kmer_dir):
        if fn.endswith(".for"):
            rev_fn = fn.replace(".for", ".rev")
            hap, sp, _ = fn.split(".")[0].split("::")
            chrn = group_db[hap]
            full_fn = os.path.join(kmer_dir, fn)
            ref_full_fn = os.path.join(kmer_dir, rev_fn)
            if chrn not in hap_kmer_file_db:
                hap_kmer_file_db[chrn] = []
            hap_kmer_file_db[chrn].append([hap, int(sp), full_fn, ref_full_fn])

    res = []
    pool = multiprocessing.Pool(processes=threads)
    for chrn in hap_kmer_file_db:
        Message.info("\tDealing %s" % (chrn))
        for i in range(len(hap_kmer_file_db[chrn]) - 1):
            for j in range(i + 1, len(hap_kmer_file_db[chrn])):
                r = pool.apply_async(
                    kmer_jac,
                    (
                        hap_kmer_file_db[chrn][i][2],
                        hap_kmer_file_db[chrn][i][3],
                        hap_kmer_file_db[chrn][j][2],
                    ),
                )
                res.append(
                    [
                        chrn,
                        hap_kmer_file_db[chrn][i][0],
                        hap_kmer_file_db[chrn][i][1],
                        hap_kmer_file_db[chrn][j][0],
                        hap_kmer_file_db[chrn][j][1],
                        r,
                    ]
                )

    pool.close()
    pool.join()

    Message.info("Writing Jac")
    out_jac = os.path.join(out_dir, "final.jac")
    with open(out_jac, "w") as fout:
        for chrn, hap1, sp1, hap2, sp2, r in sorted(res):
            jac = r.get()
            fout.write("%s\t%s\t%d\t%s\t%d\t%f\n" % (chrn, hap1, sp1, hap2, sp2, jac))

    Message.info("Finished")
