#!/usr/bin/env python3
import argparse
import re
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Polygon
import matplotlib as mpl
# Set Type 42 font embedding
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

#Author : Pushpa Itagi
## Usage : python bafplots_3_colord_Ideo.py --infile titamCuratedBinFile/sampleName_cluterNumber.titan.ichor.cna.txt --cytoband cytoBand_hg38.txt --chr 6 --label_bands --ideogram --point_size 11 --min_depth 45
#Or directly enter the params in the script

# -------------------------
# TITAN state colors (use YOUR final mapping)
# -------------------------
TITAN_STATE_COLORS = {
    "HOMD":  "#00FF00",
    "DLOH":  "#006400",
    "NLOH":  "#006400",
    "GAIN":  "#BEBEBE",
    "ALOH":  "#006400",
    "HET":   "#BEBEBE",
    "ASCNA": "#BEBEBE",
    "BCNA":  "#BEBEBE",
    "UBCNA": "#BEBEBE",}

# -------------------------
# Helpers
# -------------------------
def norm_chr(c):
    c = str(c)
    c = re.sub(r"^chr", "", c, flags=re.IGNORECASE)
    return c

def chr_order_key(c):
    c = norm_chr(c)
    if c.isdigit():
        return (0, int(c))
    if c.upper() == "X":
        return (1, 23)
    if c.upper() == "Y":
        return (1, 24)
    if c.upper() in ("M", "MT"):
        return (2, 25)
    return (3, 999)

# -------------------------
# Load TITAN SNP table
# -------------------------
def load_titan_het_like(path, min_depth, state_col="TITANstate"):
    df = pd.read_csv(path, sep="\t", low_memory=False)

    required = ["Chr", "Position", "Depth", "AllelicRatio"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}\nFound: {list(df.columns)}")

    df = df.rename(columns={
        "Chr": "chr",
        "Position": "pos",
        "Depth": "depth",
        "AllelicRatio": "baf"
    }).copy()

    df["chr"] = df["chr"].map(norm_chr)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["depth"] = pd.to_numeric(df["depth"], errors="coerce")
    df["baf"] = pd.to_numeric(df["baf"], errors="coerce")

    df = df.dropna(subset=["chr", "pos", "depth", "baf"])
    df = df[(df["depth"] >= min_depth) & (df["baf"].between(0, 1))].copy()

    # mirrored allele fractions
    df["altFrac"] = df["baf"]
    df["refFrac"] = 1.0 - df["baf"]

    # normalize state column if present
    if state_col in df.columns:
        df[state_col] = df[state_col].astype(str).str.strip()
        df["LOH"] = df[state_col].map(TITAN_STATE_COLORS).map(
            lambda color: "LOH" if color in ["#006400", "#00FF00"] else "NotLOH"
        ).fillna("NotLOH")
    return df

# -------------------------
# Genome x-axis (concatenate) - your original
# -------------------------
def add_genome_x(df, pos_col="pos"):
    chrs = sorted(df["chr"].unique(), key=chr_order_key)
    maxpos = df.groupby("chr")[pos_col].max().reindex(chrs)
    offsets = maxpos.cumsum().shift(fill_value=0)

    df = df.copy()
    df["x"] = df[pos_col] + df["chr"].map(offsets)
    centers = (offsets + maxpos / 2.0).to_dict()
    return df, chrs, centers

# -------------------------
# Cytobands + ideogram drawing (hg38 UCSC cytoBand format)
# -------------------------
def load_cytobands(cytoband_path):
    cb = pd.read_csv(
        cytoband_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "band", "stain"],
        dtype={"chrom": str, "start": int, "end": int, "band": str, "stain": str},
    ).copy()
    cb["chr"] = cb["chrom"].map(norm_chr)
    cb["stain"] = cb["stain"].fillna("gneg").astype(str)
    cb["band"] = cb["band"].astype(str)
    return cb

def stain_to_facecolor(stain):
    s = str(stain).lower()
    if s == "gneg":
        return "white"
    if s.startswith("gpos"):
        m = re.findall(r"gpos(\d+)", s)
        p = int(m[0]) if m else 50
        g = max(0.1, 1.0 - (p / 100.0) * 0.9)  # 1.0..0.1 grayscale
        return (g, g, g)
    if s in ("gvar", "stalk"):
        return (0.85, 0.85, 0.85)
    if s == "acen":
        return None
    return (0.9, 0.9, 0.9)

def draw_ideogram(ax_ideo, cb, chr_focus, height=1.0, edge_lw=2.0,
                  label_bands=True, label_min_bp=1_000_000,
                  label_rotation=90, label_fontsize=6):
    """
    Draw a single-chromosome ideogram in the SAME x-coordinates as your chr plot.
    IMPORTANT: This assumes your x-axis is basepair position on that chromosome.
    """
    chr_focus = norm_chr(chr_focus)
    sub = cb[cb["chr"] == chr_focus].copy()
    if sub.empty:
        return

    # headroom so labels don't clip
    ax_ideo.set_ylim(0, 1.30 * height)
    ax_ideo.set_yticks([])
    ax_ideo.set_xticks([])
    ax_ideo.set_frame_on(False)

    # Draw rectangles and labels
    for _, r in sub.iterrows():
        start = float(r["start"])
        end   = float(r["end"])
        stain = str(r["stain"]).lower()

        if stain == "acen":
            continue

        fc = stain_to_facecolor(stain)
        ax_ideo.add_patch(
            Rectangle(
                (start, 0.15 * height),
                end - start,
                0.70 * height,
                facecolor=fc,
                edgecolor="black",
                linewidth=edge_lw,
                zorder=2,
            )
        )

        if label_bands and (end - start) >= label_min_bp:
            mid = (start + end) / 2.0
            ax_ideo.text(
                mid,
                1.05 * height,
                str(r["band"]),
                ha="center",
                va="bottom",
                fontsize=label_fontsize,
                rotation=label_rotation,
                clip_on=False,
                zorder=10,
            )

    # Centromere triangles
    acen = sub[sub["stain"].str.lower() == "acen"].sort_values("start")
    if len(acen) >= 2:
        left, right = acen.iloc[0], acen.iloc[1]
        l0 = float(left["start"]); l1 = float(left["end"])
        r0 = float(right["start"]); r1 = float(right["end"])
        mid = (l1 + r0) / 2.0

        for tri in (
            [[l0, 0.15*height], [l1, 0.15*height], [mid, 0.50*height]],
            [[l0, 0.85*height], [l1, 0.85*height], [mid, 0.50*height]],
            [[r0, 0.15*height], [r1, 0.15*height], [mid, 0.50*height]],
            [[r0, 0.85*height], [r1, 0.85*height], [mid, 0.50*height]],
        ):
            ax_ideo.add_patch(
                Polygon(
                    tri,
                    closed=True,
                    facecolor=(0.8, 0.2, 0.2),
                    edgecolor="black",
                    linewidth=edge_lw,
                    zorder=3,
                )
            )

# -------------------------
# Plot (color by TITAN state)
# -------------------------
def plot_allelic_ratios(df, mode="genome", chr_focus=None,
                        point_size=1, alpha=0.7,
                        downsample=None,
                        outfile="baf_plot.pdf",
                        sample_name=None,
                        state_col="TITANstate",
                        show_legend=True,
                        ideogram=False,
                        cytoband_path=None,
                        label_bands=True,
                        label_min_bp=1_000_000):
    """
    Colors points by TITAN state (your mapping). Adds an ideogram below (chr mode only).
    """

    if chr_focus is not None:
        chr_focus = norm_chr(chr_focus)
        df = df[df["chr"] == chr_focus].copy()
        mode = "chr"

    if downsample is not None and len(df) > downsample:
        df = df.sample(downsample, random_state=1).sort_values(["chr", "pos"])

    # Per-point colors
    if state_col in df.columns:
        states = df[state_col].astype(str)
        df = df.copy()
        df["_color"] = states.map(TITAN_STATE_COLORS).fillna("black")
    else:
        df = df.copy()
        df["_color"] = "black"

    # Make figure
    ###fig = plt.figure(figsize=(14, 4.5) if mode == "genome" else (12, 4.5))
    fig = plt.figure(figsize=(7, 7) if mode == "genome" else (5, 5))
    ax = fig.add_subplot(111)

    if mode == "genome":
        df2, chrs, centers = add_genome_x(df, pos_col="pos")

        ax.scatter(df2["x"], df2["refFrac"], s=point_size, alpha=alpha,
                   c=df2["_color"], edgecolors="none")
        ax.scatter(df2["x"], df2["altFrac"], s=point_size, alpha=alpha,
                   c=df2["_color"], edgecolors="none")

        ax.set_xticks([])
        ax.set_xlabel("Genome position (chromosomes concatenated)", fontsize=16)

        if ideogram:
            print("[note] ideogram is only supported in --mode chr (to avoid changing your genome axis).")

    else:
        ax.scatter(df["pos"], df["refFrac"], s=point_size, alpha=alpha,
                   c=df["_color"], edgecolors="none")
        ax.scatter(df["pos"], df["altFrac"], s=point_size, alpha=alpha,
                   c=df["_color"], edgecolors="none")
        ##ax.set_xlabel(f"Position on chr{chr_focus}", fontsize=16)
        
        # Remove x-axis ticks
        ax.set_xticks([])

        # ---- Ideogram below (minimal change)
        if ideogram:
            if cytoband_path is None:
                raise ValueError("To use --ideogram, provide --cytoband (hg38 cytoBand file).")
            cb = load_cytobands(cytoband_path)

            # Put ideogram in an inset axis below x-axis
            ax_ideo = ax.inset_axes([0.0, -0.22, 1.0, 0.12], transform=ax.transAxes)
            ax_ideo.set_xlim(ax.get_xlim())
            draw_ideogram(
                ax_ideo, cb, chr_focus,
                height=1.0, edge_lw=1.0,
                label_bands=label_bands,
                label_min_bp=label_min_bp,
                label_rotation=90,
                label_fontsize=6,
            )

            # Make room for it
            plt.subplots_adjust(bottom=0.28)

    ax.set_ylim(-0.02, 1.02)
    ax.set_ylabel("Allele fraction", fontsize=16)
    ax.tick_params(axis="both", which="major", labelsize=14)
    
    # Make top and bottom spines thicker
    ax.spines['top'].set_linewidth(2.0)
    ax.spines['bottom'].set_linewidth(2.0)
    ax.spines['left'].set_linewidth(2.0)
    ax.spines['right'].set_linewidth(2.0)

    title = "Allelic ratios on Chr"
    if sample_name:
        title = f"{sample_name} - {title} -{chr_focus if chr_focus else 'genome'}"
    ax.set_title(title, fontsize=16)

    # Legend (optional)
    if show_legend and state_col in df.columns:
        present_states = [s for s in TITAN_STATE_COLORS.keys() if (df[state_col] == s).any()]
        legend_elems = [
            Line2D([0], [0], marker='o', color='w', label=s,
                   markerfacecolor=TITAN_STATE_COLORS[s], markersize=6)
            for s in present_states
        ]
        if legend_elems:
            ax.legend(handles=legend_elems, frameon=False, ncol=min(5, len(legend_elems)),
                      loc="upper right", fontsize=9)

    plt.tight_layout()
    plt.savefig(outfile, bbox_inches="tight")
    # Also a png
    plt.savefig(outfile.replace(".pdf", ".png"), bbox_inches="tight", dpi=150)
    plt.close()
    # Save only the columns needed to reproduce the plot source data.
    keep_cols = ["chr", "pos", "depth", "baf", "altFrac", "refFrac"]
    if state_col in df.columns:
        keep_cols.append(state_col)
    if "LOH" in df.columns:
        keep_cols.append("LOH")
    if mode == "genome":
        keep_cols.append("x")

    keep_cols = [col for col in keep_cols if col in df.columns]
    plot_df = df[keep_cols].copy()
 
    data_outfile = outfile.replace(".pdf", "_data.csv").replace(".png", "_data.csv")
    plot_df.drop(columns=["TITANcall"], inplace=True, errors="ignore")
    #Also add the sample name to the output data file
    plot_df.insert(0, "sample", sample_name)    
    plot_df.to_csv(data_outfile, index=False)
    return 0;

# -------------------------
# Main
# -------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True, help="TITAN SNP table")
    ap.add_argument("--min_depth", type=int, default=10)
    ap.add_argument("--mode", choices=["genome", "chr"], default="genome")
    ap.add_argument("--chr", dest="chr_focus", default=None)
    ap.add_argument("--point_size", type=float, default=1)
    ap.add_argument("--alpha", type=float, default=1.0)
    ap.add_argument("--downsample", type=int, default=200000)

    ap.add_argument("--state_col", default="TITANcall",
                    help="Column name holding TITAN state calls (default TITANcall).")

    # ideogram
    ap.add_argument("--ideogram", action="store_true",
                    help="Add cytoband ideogram BELOW the x-axis (chr mode only).")
    ap.add_argument("--cytoband", default=None,
                    help="hg38 cytoBand file (UCSC format). Example: /mnt/data/cytoBand_hg38.txt")
    ap.add_argument("--label_bands", action="store_true",
                    help="Label cytobands on ideogram (p21.2 etc).")
    ap.add_argument("--label_min_bp", type=int, default=1_000_000,
                    help="Only label bands >= this size (bp). For full chr, 1e6 is ok; increase if too crowded.")

    # output
    ap.add_argument("--outfile", default="output_dir/baf_plot_color_ideo.pdf")
    ap.add_argument("--legend", action="store_true", help="Show legend for TITAN states (can be busy).")

    args = ap.parse_args()

    sample_name = os.path.basename(args.infile).split("_cluster")[0]
    if args.outfile == "output_dir/baf_plot_color_ideo.pdf":
        if args.chr_focus:
            args.outfile = f"output_dir/baf_plot_{sample_name}_chr{norm_chr(args.chr_focus)}_color_ideo.pdf"
        else:
            args.outfile = f"output_dir/baf_plot_{sample_name}_genome_color_ideo.pdf"

    df = load_titan_het_like(args.infile, min_depth=args.min_depth, state_col=args.state_col)

    plot_allelic_ratios(
        df,
        mode=args.mode,
        chr_focus=args.chr_focus,
        point_size=args.point_size,
        alpha=args.alpha,
        downsample=args.downsample,
        outfile=args.outfile,
        sample_name=sample_name,
        state_col=args.state_col,
        show_legend=args.legend,
        ideogram=args.ideogram,
        cytoband_path=args.cytoband,
        label_bands=args.label_bands,
        label_min_bp=args.label_min_bp,
    )

    print(f"Wrote: {args.outfile}")

if __name__ == "__main__":
    # Batch processing parameters
    batch_mode = True  # Set to False to use command line args
    
    if batch_mode:
        # Your batch parameters
        in_dir = "/path/to/CuratedSoultions/"
        starts_with = "00-010"  # Add the common prefix of your files here, e.g. "00-010" to match all patients starting with 00-010
        
        point_size = 11
        min_depth = 20 #   min_depth = 35 for 17-026 chr6, min_depth = 20 for all FF samples,  min_depth = 25 for chr 19 17-020
        chr_focus = "6"
        
        # Recursively find files
        for root, dirs, files in os.walk(in_dir):
            for file in files:
                if file.startswith(starts_with) and file.endswith(".titan.ichor.cna.txt"):
                    infile = os.path.join(root, file)
                    print(f"Processing {infile}...")
                    
                    # Generate sample name and output file
                    sample_name = os.path.basename(file).split("_cluster")[0]
                    outfile = f"output_dir/baf_plot_{sample_name}_chr{chr_focus}_color_ideo.pdf"
                    
                    # Create args namespace with all parameters
                    args = argparse.Namespace(
                        infile=infile,
                        min_depth=min_depth,
                        mode="chr",
                        chr_focus=chr_focus,
                        point_size=point_size,
                        alpha=1.0,
                        downsample=None,
                        state_col="TITANcall",
                        ideogram=True,
                        cytoband="cytoBand_hg38.txt",
                        label_bands=False,
                        label_min_bp=1_000_000,
                        outfile=outfile,
                        legend=False,
                    )
                    
                    # Process file
                    df = load_titan_het_like(args.infile, min_depth=args.min_depth, state_col=args.state_col)
                    plot_allelic_ratios(
                        df,
                        mode=args.mode,
                        chr_focus=args.chr_focus,
                        point_size=args.point_size,
                        alpha=args.alpha,
                        downsample=args.downsample,
                        outfile=args.outfile,
                        sample_name=sample_name,
                        state_col=args.state_col,
                        show_legend=args.legend,
                        ideogram=args.ideogram,
                        cytoband_path=args.cytoband,
                        label_bands=args.label_bands,
                        label_min_bp=args.label_min_bp,
                    )
                    print(f"Wrote: {args.outfile}")
    else:
        # Use command line arguments
        main()