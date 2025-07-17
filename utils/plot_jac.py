import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from utils.message import Message

mpl.use("Agg")


def plot_jac(
    in_jac, group_file, out_dir, cmap="YlOrRd", pic_fmt="pdf", log_scale=False
):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    hap_order_db = {}
    with open(group_file, "r") as fin:
        for line in fin:
            chrn, hap = line.strip().split()
            if chrn not in hap_order_db:
                hap_order_db[chrn] = []
            hap_order_db[chrn].append(hap)

    Message.info("Loading Jac")
    mat_size_db = {}
    jac_db = {}
    pos_db = {}
    with open(in_jac, "r") as fin:
        for line in fin:
            data = line.strip().split()
            chrn = data[0]
            hap1 = data[1]
            sp1 = int(data[2])
            hap2 = data[3]
            sp2 = int(data[4])
            jac = float(data[5])
            if chrn not in mat_size_db:
                mat_size_db[chrn] = {}
            if hap1 not in mat_size_db[chrn]:
                mat_size_db[chrn][hap1] = []
            mat_size_db[chrn][hap1].append(sp1)
            if chrn not in jac_db:
                jac_db[chrn] = []
            jac_db[chrn].append([hap1, hap2, sp1, sp2, jac])
            if chrn not in pos_db:
                pos_db[chrn] = {}
            if hap1 not in pos_db[chrn]:
                pos_db[chrn][hap1] = set()
            pos_db[chrn][hap1].add(sp1)
            if hap2 not in pos_db[chrn]:
                pos_db[chrn][hap2] = set()
            pos_db[chrn][hap2].add(sp2)

    Message.info("Plotting jac")
    for chrn in mat_size_db:
        Message.info("\tInit %s..." % chrn)
        idx_db = {}
        idx = 0
        for hap in hap_order_db[chrn]:
            if hap not in idx_db:
                idx_db[hap] = {}
            for sp in sorted(pos_db[chrn][hap]):
                idx_db[hap][sp] = idx
                idx += 1

        mat = [[0 for _ in range(idx)] for _ in range(idx)]
        for _ in range(idx):
            mat[_][_] = 1
        for hap1, hap2, sp1, sp2, jac in jac_db[chrn]:
            idx1 = idx_db[hap1][sp1]
            idx2 = idx_db[hap2][sp2]
            mat[idx1][idx2] = jac
            mat[idx2][idx1] = jac

        Message.info("\tplotting")
        plt.figure(figsize=(10, 10), dpi=300)

        X = []
        haps = []
        last_border = 0
        for hap in hap_order_db[chrn]:
            haps.append(hap)
            border = idx_db[hap][sorted(idx_db[hap], reverse=True)[0]]
            if len(haps) < len(idx_db):
                plt.plot((border, border), (0, idx), color="grey", linestyle=":", lw=1)
                plt.plot((border, border), (idx, 0), color="grey", linestyle=":", lw=1)
                plt.plot((0, idx), (border, border), color="grey", linestyle=":", lw=1)
                plt.plot((idx, 0), (border, border), color="grey", linestyle=":", lw=1)
            X.append((border - last_border) / 2.0 + last_border)
            last_border = border

        if log_scale:
            for i in range(len(mat)):
                for j in range(len(mat[i])):
                    if mat[i][j] >= 0:
                        mat[i][j] = np.log10(mat[i][j] * 100.0 + 1)
                    elif mat[i][j] < 0:
                        mat[i][j] = -np.log10(-mat[i][j] * 100.0 + 1)
            vmin = -np.log10(101)
            vmax = np.log10(101)
        else:
            vmin = -1
            vmax = 1
        ax = plt.gca()
        hmap = ax.imshow(
            mat,
            interpolation="nearest",
            origin="lower",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            aspect="equal",
        )

        cbar = plt.colorbar(
            mappable=hmap,
            cax=None,
            ax=None,
            shrink=0.5,
            label=(
                "Jaccard Similarity (scaled)" if log_scale else "Jaccard Similarity"
            ),
        )

        cbar.set_ticks(
            ticks=(
                [
                    (-1 if _ < 0 else 1) * np.log10(np.abs(_) + 1)
                    for _ in range(-100, 101, 25)
                ]
                if log_scale
                else [_ / 100.0 for _ in range(-100, 101, 25)]
            ),
            labels=["%d%%" % _ for _ in range(-100, 101, 25)],
        )
        plt.xticks(X, haps)
        plt.yticks(X, haps)
        plt.tick_params(labelsize=12)
        for ticks in ax.get_xticklabels():
            ticks.set_rotation(90)
        for ticks in ax.get_yticklabels():
            ticks.set_rotation(0)
        plt.title(chrn, y=1.01, fontsize=20)
        pic_fn = os.path.join(out_dir, "%s.%s" % (chrn, pic_fmt))
        plt.savefig(pic_fn, bbox_inches="tight", dpi=300)
        plt.close("all")

    Message.info("Plotted")
