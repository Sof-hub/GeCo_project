# -*- coding: utf-8 -*-
"""
Full Comparative Genomic Analysis (v2)
Author: Sofia
Date: 2025-11-02

Description:
This script performs an integrated comparative genomic analysis across multiple bacterial genomes.
All analyses are based on complete genome sequences and coding sequences (CDS) extracted from
GenBank files (.gbff), ensuring nucleotide-level biological validity.

Analyses performed:

1. General genome features analysis
   - Genome size
   - GC content
   - CDS count
   - Basic compositional statistics

2. GC Skew and cumulative GC Skew analysis
   - Strand compositional asymmetry
   - Detection of replication-associated biases

3. Replication origin and terminus inference
   - Identification of putative oriC and ter sites
   - Symmetry and replication architecture analysis

4. Functional annotation comparison
   - COG functional category distribution
   - Cross-organism functional profile comparison

5. Codon usage analysis (CDS nucleotide sequences only)
   - Global codon usage frequency
   - Relative Synonymous Codon Usage (RSCU)

6. Codon Adaptation Index (CAI) analysis
   - CAI computation using a reference organism
   - Assessment of translational optimization

7. Dinucleotide frequency analysis
   - Genome-wide dinucleotide composition
   - Detection of compositional and mutational biases

8. Comparative visualization and summary outputs
   - Heatmaps (RSCU, dinucleotide frequencies)
   - Tabulated comparative results
   - Integrated multi-panel summary dashboard

"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# ============================================================
# CONFIGURATION
# ============================================================

# Changer le working directory
os.chdir(r"C:\Users\sofia\Desktop\GeCoproject")

OUTPUT_DIR = "results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

genomes = {
    "Apilactobacillus kunkeei": {
        "fna": r"C:\Users\sofia\Desktop\GeCoproject\data\1_GCF_019575995.1_ASM1957599v1_genomic.fna",
        "faa": r"C:\Users\sofia\Desktop\GeCoproject\data\1_GCF_019575995.1_ASM1957599v1_protein.faa",
        "cusp": r"C:\Users\sofia\Desktop\GeCoproject\data\1_GCF_019575995.1_ASM1957599v1_genomic_cusp.txt",
        "emapper": r"C:\Users\sofia\Desktop\GeCoproject\data\1_GCF_019575995.1_out.emapper.annotations.xlsx",
        "genbank": r"C:\Users\sofia\Desktop\GeCoproject\data\1_GCF_019575995.1_ASM1957599v1.gbff"
    },
    "Apilactobacillus bombintestini": {
        "fna": r"C:\Users\sofia\Desktop\GeCoproject\data\2_GCF_003627035.1_ASM362703v1_genomic.fna",
        "faa": r"C:\Users\sofia\Desktop\GeCoproject\data\2_GCF_003627035.1_ASM362703v1_protein.faa",
        "cusp": r"C:\Users\sofia\Desktop\GeCoproject\data\2_GCF_003627035.1_ASM362703v1_genomic_cusp.txt",
        "emapper": r"C:\Users\sofia\Desktop\GeCoproject\data\2_GCF_003627035.1_out.emapper.annotations.xlsx",
        "genbank": r"C:\Users\sofia\Desktop\GeCoproject\data\2_GCF_003627035.1_ASM362703v1.gbff"
    },
    "Apilactobacillus apinorum": {
        "fna": r"C:\Users\sofia\Desktop\GeCoproject\data\3_GCF_946888465.1_Fhon13_genomic.fna",
        "faa": r"C:\Users\sofia\Desktop\GeCoproject\data\3_GCF_946888465.1_Fhon13_protein.faa",
        "cusp": r"C:\Users\sofia\Desktop\GeCoproject\data\3_GCF_946888465.1_Fhon13_genomic_cusp.txt",
        "emapper": r"C:\Users\sofia\Desktop\GeCoproject\data\3_GCF_946888465.1_out.emapper.annotations.xlsx",
        "genbank": r"C:\Users\sofia\Desktop\GeCoproject\data\3_GCF_946888465.1_Fhon13.gbff"
    },
    "Bombilactobacillus folatiphilus": {
        "fna": r"C:\Users\sofia\Desktop\GeCoproject\data\4_GCF_023380265.1_ASM2338026v1_genomic.fna",
        "faa": r"C:\Users\sofia\Desktop\GeCoproject\data\4_GCF_023380265.1_ASM2338026v1_protein.faa",
        "cusp": r"C:\Users\sofia\Desktop\GeCoproject\data\4_GCF_023380265.1_ASM2338026v1_genomic_cusp.txt",
        "emapper": r"C:\Users\sofia\Desktop\GeCoproject\data\4_GCF_023380265.1_out.emapper.annotations.xlsx",
        "genbank": r"C:\Users\sofia\Desktop\GeCoproject\data\4_GCF_023380265.1_ASM2338026v1.gbff"
    },
    "Bombilactobacillus thymidiniphilus": {
        "fna": r"C:\Users\sofia\Desktop\GeCoproject\data\5_GCF_023380245.1_ASM2338024v1_genomic.fna",
        "faa": r"C:\Users\sofia\Desktop\GeCoproject\data\5_GCF_023380245.1_ASM2338024v1_protein.faa",
        "cusp": r"C:\Users\sofia\Desktop\GeCoproject\data\5_GCF_023380245.1_ASM2338024v1_genomic_cusp.txt",
        "emapper": r"C:\Users\sofia\Desktop\GeCoproject\data\5_GCF_023380245.1_out.emapper.annotations.xlsx",
        "genbank": r"C:\Users\sofia\Desktop\GeCoproject\data\5_GCF_023380245.1_ASM2338024v1.gbff"
    }
}

# ============================================================
# FUNCTIONS
# ============================================================

def gc_skew(seq, window=None, step=None, cumulative=False):
    """Compute GC skew and optionally cumulative GC skew."""
    n = len(seq)
    if window is None:
        window = max(10000, n // 1000)
    if step is None:
        step = max(1000, window // 10)
    skews, positions = [], []
    for i in range(0, n - window, step):
        frag = seq[i:i + window]
        g, c = frag.count("G"), frag.count("C")
        skew = (g - c) / (g + c) if (g + c) > 0 else 0
        skews.append(skew)
        positions.append(i + window // 2)
    if cumulative:
        skews = np.cumsum(skews)
    return positions, skews

import itertools
from collections import Counter

def compute_dinucleotide_freq(fasta_file):
    dinucs = ["".join(p) for p in itertools.product("ACGT", repeat=2)]
    counts = Counter()
    total = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        for i in range(len(seq) - 1):
            dinuc = seq[i:i+2]
            if all(b in "ACGT" for b in dinuc):
                counts[dinuc] += 1
                total += 1
    return {d: counts[d] / total if total > 0 else 0 for d in dinucs}


def compute_codon_usage(fasta_file):
    codons = ["".join(p) for p in itertools.product("ACGT", repeat=3)]
    counts = Counter()
    total = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        for i in range(0, len(seq)-2, 3):
            codon = seq[i:i+3]
            if all(b in "ACGT" for b in codon):
                counts[codon] += 1
                total += 1
    return {c: counts[c] / total if total > 0 else 0 for c in codons}


def compute_RSCU(codon_freqs):
    aa_table = {
        "F":["TTT","TTC"], "L":["TTA","TTG","CTT","CTC","CTA","CTG"],
        "I":["ATT","ATC","ATA"], "M":["ATG"], "V":["GTT","GTC","GTA","GTG"],
        "S":["TCT","TCC","TCA","TCG","AGT","AGC"], "P":["CCT","CCC","CCA","CCG"],
        "T":["ACT","ACC","ACA","ACG"], "A":["GCT","GCC","GCA","GCG"],
        "Y":["TAT","TAC"], "H":["CAT","CAC"], "Q":["CAA","CAG"],
        "N":["AAT","AAC"], "K":["AAA","AAG"], "D":["GAT","GAC"],
        "E":["GAA","GAG"], "C":["TGT","TGC"], "W":["TGG"],
        "R":["CGT","CGC","CGA","CGG","AGA","AGG"],
        "G":["GGT","GGC","GGA","GGG"], "*":["TAA","TAG","TGA"]
    }
    rscu = {}
    for aa, codons in aa_table.items():
        total = sum(codon_freqs.get(c,0) for c in codons)
        n = len(codons)
        for c in codons:
            exp = total/n if n > 0 else 0
            rscu[c] = codon_freqs.get(c,0)/exp if exp > 0 else 0
    return rscu


def compute_CAI(fasta_file, ref_usage):
    usage = compute_codon_usage(fasta_file)
    w = [usage[c]/ref_usage[c] for c in usage if c in ref_usage and ref_usage[c] > 0 and usage[c] > 0]
    return float(np.exp(np.mean(np.log(w)))) if w else 0

def is_nucleotide_fasta(fasta_file, max_records=10):
    """
    Vérifie si un FASTA contient des séquences nucléotidiques (ACGTN uniquement).
    Analyse les max_records premières séquences.
    """
    allowed = set("ACGTN")
    for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
        if i >= max_records:
            break
        seq = set(str(record.seq).upper())
        if not seq.issubset(allowed):
            return False
    return True

def extract_cds_from_gbff(gbff_file, output_fasta):
    cds_records = []
    cds_count = 0

    for record in SeqIO.parse(gbff_file, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue

            if "translation" not in feature.qualifiers:
                continue

            try:
                seq = feature.extract(record.seq).upper()
            except Exception:
                continue

            if len(seq) % 3 != 0:
                continue

            if not set(seq).issubset({"A", "T", "G", "C"}):
                continue

            locus = feature.qualifiers.get("locus_tag", ["unknown"])[0]
            protein_id = feature.qualifiers.get("protein_id", ["unknown"])[0]

            cds_records.append(
                SeqRecord(
                    Seq(str(seq)),
                    id=protein_id,
                    description=f"locus_tag={locus}"
                )
            )
            cds_count += 1

    SeqIO.write(cds_records, output_fasta, "fasta")
    print(f"✅ {cds_count} CDS extraites → {output_fasta}")


# ============================================================
# 1. GC SKEW COMPARATIVE PLOT
# ============================================================

plt.figure(figsize=(10, 6))
colors = ["red", "blue", "green", "orange", "purple"]

for (org, paths), color in zip(genomes.items(), colors):
    record = next(SeqIO.parse(paths["fna"], "fasta"))
    seq = str(record.seq).upper()
    _, skews = gc_skew(seq)
    plt.plot(skews, label=org, color=color)

plt.title("GC Skew comparison (window=10,000 bp)")
plt.xlabel("Genome position (x1000 bp)")
plt.ylabel("GC skew")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "GC_skew_all.png"))
plt.close()
print("✅ GC skew comparison plot saved.")

fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(10, 12), sharex=True)

colors = ["red", "blue", "green", "orange", "purple"]

for (ax, ((org, paths), color)) in zip(axes, zip(genomes.items(), colors)):
    record = next(SeqIO.parse(paths["fna"], "fasta"))
    seq = str(record.seq).upper()
    _, skews = gc_skew(seq)

    ax.plot(skews, color=color)
    ax.set_title(org)
    ax.set_ylabel("GC skew")

axes[-1].set_xlabel("Genome position (x1000 bp)")
fig.suptitle("GC Skew comparison (window=10,000 bp)")
fig.tight_layout(rect=[0, 0, 1, 0.96])

plt.savefig(os.path.join(OUTPUT_DIR, "GC_skew_all_subplots.png"))
plt.close()
print("✅ GC skew comparison subplot figure saved.")


# ============================================================
# 2. MULTI-PANEL GC SKEW WITH ORI/TER & SYMMETRY ANALYSIS
# ============================================================

fig, axes = plt.subplots(3, 2, figsize=(12, 10))
axes = axes.flatten()
replication_sites = []

for i, (org, paths) in enumerate(genomes.items()):
    try:
        seq_record = next(SeqIO.parse(paths["fna"], "fasta"))
        seq = str(seq_record.seq).upper()
        genome_size = len(seq)

        # Coding GC
        coding_gc = None
        with open(paths["cusp"]) as f:
            for line in f:
                if line.startswith("#Coding GC"):
                    coding_gc = float(line.strip().split()[2].replace("%", ""))
                    break

        # Cumulative GC skew
        positions, cskews = gc_skew(seq, cumulative=True)
        ori_pos = positions[np.argmin(cskews)]
        ter_pos = positions[np.argmax(cskews)]

        # Plot cumulative skew + annotations
        axes[i].plot(positions, cskews, color="teal")
        axes[i].axvline(ori_pos, color="red", linestyle="--", label="ori")
        axes[i].axvline(ter_pos, color="blue", linestyle="--", label="ter")
        y_min, y_max = np.min(cskews), np.max(cskews)
        axes[i].text(ori_pos, y_min, f"ori\n{ori_pos:,} bp", color="red", fontsize=8,
                     ha="right", va="bottom", rotation=90)
        axes[i].text(ter_pos, y_max, f"ter\n{ter_pos:,} bp", color="blue", fontsize=8,
                     ha="left", va="bottom", rotation=90)

        # Symmetry metrics
        ori_ter_dist = abs(ter_pos - ori_pos)
        ori_ter_dist_pct = (ori_ter_dist / genome_size) * 100
        midpoint = (ori_pos + ter_pos) / 2
        symmetry_index = 100 * (1 - abs(midpoint - genome_size / 2) / (genome_size / 2))

        replication_sites.append({
            "Organism": org,
            "Genome_size(bp)": genome_size,
            "Origin_position(bp)": ori_pos,
            "Terminus_position(bp)": ter_pos,
            "Ori–Ter_distance(bp)": ori_ter_dist,
            "Ori–Ter_distance(%)": ori_ter_dist_pct,
            "Symmetry_index(%)": symmetry_index,
            "Coding_GC(%)": coding_gc
        })

        title = (f"{org}\n{genome_size/1e6:.2f} Mb | GC={coding_gc or 0:.2f}% | "
                 f"Sym={symmetry_index:.1f}%")
        axes[i].set_title(title, fontsize=9)
        axes[i].set_xlabel("Genome position (bp)")
        axes[i].set_ylabel("Cumulative GC skew")
        axes[i].grid(True, linestyle="--", alpha=0.3)
        axes[i].legend(fontsize=8)

    except Exception as e:
        print(f"⚠️ Error processing {org}: {e}")

for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.suptitle("Cumulative GC skew with replication origin (ori) and terminus (ter)", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])
cskew_fig = os.path.join(OUTPUT_DIR, "Cumulative_GC_skew_with_labels.png")
plt.savefig(cskew_fig, dpi=300)
plt.close()
print(f"✅ Labeled cumulative GC skew saved: {cskew_fig}")

# Save replication data
rep_df = pd.DataFrame(replication_sites)
rep_table = os.path.join(OUTPUT_DIR, "Replication_sites.tsv")
rep_df.to_csv(rep_table, sep="\t", index=False)
print(f"✅ Replication metrics saved: {rep_table}")

# Plot replication symmetry summary
plt.figure(figsize=(8, 5))
plt.bar(rep_df["Organism"], rep_df["Symmetry_index(%)"], color="mediumseagreen")
plt.ylabel("Replication symmetry index (%)")
plt.title("Replication symmetry across genomes")
plt.ylim(0, 100)
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
sym_plot = os.path.join(OUTPUT_DIR, "Replication_symmetry_summary.png")
plt.savefig(sym_plot, dpi=300)
plt.close()
print(f"✅ Replication symmetry summary plot saved: {sym_plot}")


# ============================================================
# 3. GENERAL FEATURES TABLE
# ============================================================

features = []

for org, paths in genomes.items():
    # --- CDS depuis .faa ---
    faa = list(SeqIO.parse(paths["faa"], "fasta"))
    cds_count = len(faa)
    cds_lengths = [len(p.seq) for p in faa]
    total_cds_len = sum(cds_lengths)
    mean_cds_len = total_cds_len / cds_count
    max_cds_len = max(cds_lengths)

    # --- Génome + tRNA/rRNA depuis GBFF (multi-contigs possible) ---
    records = list(SeqIO.parse(paths["genbank"], "genbank"))
    genome_size = sum(len(rec.seq) for rec in records)
    trna_count = sum(1 for rec in records for f in rec.features if f.type == "tRNA")
    rrna_count = sum(1 for rec in records for f in rec.features if f.type == "rRNA")
    coding_density = 100 * total_cds_len / genome_size


    # --- GC content (coding) depuis cusp ---
    gc_coding = None
    with open(paths["cusp"]) as f:
        for line in f:
            if line.startswith("#Coding GC"):
                gc_coding = float(line.strip().split()[2].replace("%", ""))
                break

    # --- Annotation eggNOG/COG ---
    df = pd.read_excel(paths["emapper"], comment="#", engine="openpyxl")
    df_clean = df[~df["Description"].str.contains("hypothetical|unknown",
                                                  case=False, na=False)]
    proteins_pred_fun = len(df_clean)
    proteins_with_cog = df["COG_category"].notna().sum()

    features.append({
        "Organism": org,
        "Genome_size(bp)": genome_size,
        "GC_coding(%)": gc_coding,
        "Predicted_CDS": cds_count,
        "Coding_density(%)": coding_density,
        "Mean_CDS_length(bp)": mean_cds_len,
        "Max_CDS_length(bp)": max_cds_len,
        "tRNA": trna_count,
        "rRNA": rrna_count,
        "Proteins_predicted_function": proteins_pred_fun,
        "Proteins_with_COG": proteins_with_cog,
    })

features_df = pd.DataFrame(features)
features_out = os.path.join(OUTPUT_DIR, "General_features.tsv")
features_df.to_csv(features_out, sep="\t", index=False)
print("✅ General features table saved.")

# ============================================================
# Codon usage, CAI and RSCU (NUCLEOTIDE CDS ONLY)
# ============================================================


CDS_DIR = r"C:\Users\sofia\Desktop\GeCoproject\data\CDS_ffn"
os.makedirs(CDS_DIR, exist_ok=True)

for org, paths in genomes.items():
    gbff = paths["genbank"]

    if not os.path.exists(gbff):
        print(f" {org}: fichier gbff introuvable")
        continue

    out_fasta = os.path.join(
        CDS_DIR,
        org.replace(" ", "_") + "_CDS.ffn"
    )

    extract_cds_from_gbff(
        gbff_file=gbff,
        output_fasta=out_fasta
    )

    # On enrichit le dictionnaire existant
    genomes[org]["cds_fasta"] = out_fasta


ref_fasta = genomes["Apilactobacillus kunkeei"]["cds_fasta"]

if not is_nucleotide_fasta(ref_fasta):
    raise ValueError(
        "FASTA de référence CAI non nucléotidique : "
        "impossible de calculer CAI/RSCU scientifiquement."
    )

ref_usage = compute_codon_usage(ref_fasta)


dinuc_data = {}
cai_data = []

for org, paths in genomes.items():
    fasta = paths.get("cds_fasta")


    if not is_nucleotide_fasta(fasta):
        print(f" {org}: FASTA protéique → CAI & dinucléotides ignorés")
        continue

    cai = compute_CAI(fasta, ref_usage)
    dinuc = compute_dinucleotide_freq(fasta)

    dinuc_data[org] = dinuc
    cai_data.append({"Organism": org, "CAI": round(cai, 3)})


cai_df = pd.DataFrame(cai_data)
cai_df.to_csv(os.path.join(OUTPUT_DIR, "CAI_values.tsv"),
              sep="\t", index=False)

dinuc_df = pd.DataFrame(dinuc_data).T
dinuc_df.to_csv(os.path.join(OUTPUT_DIR, "Dinucleotide_frequency.tsv"),
                sep="\t")

plt.figure(figsize=(10,6))
plt.imshow(dinuc_df, aspect="auto", cmap="viridis")
plt.colorbar(label="Frequency")
plt.xticks(range(len(dinuc_df.columns)), dinuc_df.columns, rotation=90)
plt.yticks(range(len(dinuc_df.index)), dinuc_df.index)
plt.title("Dinucleotide frequencies")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "Dinucleotide_frequencies_heatmap.png"), dpi=300)
plt.close()

for org in genomes:
    print(org, genomes[org]["cds_fasta"])


rscu_data = {}

for org, paths in genomes.items():
    fasta = paths.get("cds_fasta")


    if not is_nucleotide_fasta(fasta):
        print(f" {org}: FASTA protéique → RSCU ignoré")
        continue

    codon_usage = compute_codon_usage(fasta)
    rscu_data[org] = compute_RSCU(codon_usage)


rscu_df = pd.DataFrame(rscu_data)
rscu_out = os.path.join(OUTPUT_DIR, "RSCU_all_organisms.tsv")
rscu_df.to_csv(rscu_out, sep="\t")
print(f"✅ RSCU table saved: {rscu_out}")

plt.figure(figsize=(18,12))
plt.imshow(rscu_df, aspect="auto", cmap="coolwarm")
plt.colorbar(label="RSCU")
plt.xticks(range(len(rscu_df.columns)), rscu_df.columns, rotation=45, ha="right")
plt.yticks(range(len(rscu_df.index)), rscu_df.index)
plt.title("Relative Synonymous Codon Usage (RSCU)")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "RSCU_heatmap.png"), dpi=300)
plt.close()



# ============================================================
# 4. FUNCTIONAL COG ANALYSIS
# ============================================================

cog_all = {}
for org, paths in genomes.items():
    df = pd.read_excel(paths["emapper"], comment="#", engine="openpyxl")
    df = df.dropna(subset=["COG_category"])
    grouped = df.groupby("COG_category").size()
    percent = 100 * grouped / grouped.sum()
    cog_all[org] = percent

cog_all_df = pd.DataFrame(cog_all).fillna(0)
cog_file = os.path.join(OUTPUT_DIR, "COG_all_organisms.tsv")
cog_all_df.to_csv(cog_file, sep="\t")

cog_all_df.plot(kind="bar", figsize=(21, 14))
plt.ylabel("Percentage of proteins")
plt.title("COG category distribution across the 5 organisms")
plt.legend(title="Organisms")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "COG_all_comparison.png"))
plt.close()
print("✅ COG comparison completed.")

import numpy as np

fig, ax = plt.subplots(figsize=(22, 12))

categories = cog_all_df.index
organisms = cog_all_df.columns

x = np.arange(len(categories))
bar_width = 0.15  # 5 organismes → ~0.15

for i, org in enumerate(organisms):
    ax.bar(
        x + i * bar_width,
        cog_all_df[org],
        width=bar_width,
        label=org
    )

ax.set_xticks(x + bar_width * (len(organisms) - 1) / 2)
ax.set_xticklabels(categories, rotation=90)
ax.set_ylabel("Percentage of proteins")
ax.set_title("COG category distribution across the 5 organisms")
ax.legend(title="Organisms")

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "COG_all_comparison_grouped.png"), dpi=300)
plt.close()


# ============================================================
# 4bis. COG NCBI vs EggNOG (A. kunkeei)
# ============================================================

# Chemins vers les fichiers de A. kunkeei dans le sous-dossier data
ncbi_cog_file = os.path.join("data", "1_GCF_019575995.1_ASM1957599v1_cog_result_table.tsv")
eggnog_file   = genomes["Apilactobacillus kunkeei"]["emapper"]

# --- 4bis.1. Profil COG (NCBI) à partir du tableau téléchargé sur le site COG ---
ncbi_df = pd.read_csv(ncbi_cog_file, sep="\t")

# La colonne de catégories fonctionnelles est "Cat"
ncbi_counts = ncbi_df["Cat"].value_counts()
ncbi_percent = 100 * ncbi_counts / ncbi_counts.sum()

# --- 4bis.2. Profil COG (EggNOG) à partir du fichier emapper ---
emap_df = pd.read_excel(eggnog_file, comment="#", engine="openpyxl")

# Colonne des catégories COG dans eggNOG-mapper (adapter si nom différent)
# Ici on suppose une colonne "COG_category" qui contient les lettres (A, C, "E H", etc.)
emap_df_cog = emap_df.dropna(subset=["COG_category"])
egg_counts = emap_df_cog["COG_category"].value_counts()
egg_percent = 100 * egg_counts / egg_counts.sum()

# --- 4bis.3. Mettre NCBI et EggNOG dans la même table ---
all_cats = sorted(set(ncbi_percent.index).union(egg_percent.index))
comp_df = pd.DataFrame(index=all_cats, columns=["NCBI", "EggNOG"]).fillna(0.0)
comp_df.loc[ncbi_percent.index, "NCBI"] = ncbi_percent.values
comp_df.loc[egg_percent.index, "EggNOG"] = egg_percent.values

# Sauvegarde TSV
cog_comp_out = os.path.join(OUTPUT_DIR, "COG_comparison_Apilactobacillus_kunkeei.tsv")
comp_df.to_csv(cog_comp_out, sep="\t")
print(f"✅ COG NCBI vs EggNOG table saved: {cog_comp_out}")

# --- 4bis.4. Figure de comparaison NCBI vs EggNOG ---
plt.figure(figsize=(18, 10))      
x = np.arange(len(comp_df))      
width = 0.8                      

plt.bar(x - width/4, comp_df["NCBI"], width=width/2, label="NCBI")
plt.bar(x + width/4, comp_df["EggNOG"], width=width/2, label="EggNOG")


plt.xticks(ticks=x, labels=comp_df.index, rotation=90)
plt.ylabel("Percentage")
plt.title("COG Category Comparison (NCBI vs EggNOG) - A. kunkeei")
plt.legend()
plt.tight_layout()

cog_comp_fig = os.path.join(OUTPUT_DIR, "COG_Apilactobacillus_kunkeei_comparison.png")
plt.savefig(cog_comp_fig, dpi=300)
plt.close()
print(f"✅ COG NCBI vs EggNOG plot saved: {cog_comp_fig}")

# ============================================================
# 4ter. COG CATEGORY COMPARISON PLOT (5 ORGANISMS)
# ============================================================

plt.figure(figsize=(18, 10))

x = range(len(cog_all_df.index))
width = 0.8  

organisms_list = list(cog_all_df.columns)

for i, org in enumerate(organisms_list):
    plt.bar([v + (i - len(organisms_list)/2) * width for v in x],
            cog_all_df[org],
            width=width,
            label=org)

plt.xticks(ticks=x, labels=cog_all_df.index, rotation=90)
plt.ylabel("Percentage")
plt.xlabel("COG_category")
plt.title("COG Category Comparison across 5 Organisms")
plt.legend(title="Organisms")
plt.tight_layout()

cog_all_fig = os.path.join(OUTPUT_DIR, "COG_all_comparison.png")
plt.savefig(cog_all_fig, dpi=300)
plt.close()
print(f"✅ COG all-organisms comparison plot saved: {cog_all_fig}")



# ============================================================
# 5. INTEGRATED SUMMARY FIGURE
# ============================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

# (A) GC Skew
for (org, paths), color in zip(genomes.items(), colors):
    record = next(SeqIO.parse(paths["fna"], "fasta"))
    seq = str(record.seq).upper()
    _, skews = gc_skew(seq)
    axes[0].plot(skews, label=org, color=color)
axes[0].set_title("A. GC Skew (structural bias)")
axes[0].legend()

# (B) GC content
axes[1].bar(features_df["Organism"], features_df["GC_coding(%)"], color="skyblue")
axes[1].set_title("B. Coding GC content (%)")
axes[1].set_ylabel("GC%")
axes[1].tick_params(axis='x', rotation=45)

# (C) Functional proteins
axes[2].bar(features_df["Organism"], features_df["Proteins_predicted_function"], color="orange")
axes[2].set_title("C. Proteins with predicted function")
axes[2].tick_params(axis='x', rotation=45)

# (D) Average COG profile
cog_all_df.mean(axis=1).plot(kind="bar", ax=axes[3], color="green")
axes[3].set_title("D. Average COG distribution")
axes[3].set_ylabel("Average %")
axes[3].tick_params(axis='x', rotation=90)

plt.tight_layout()
summary_path = os.path.join(OUTPUT_DIR, "Summary_comparative_dashboard.png")
plt.savefig(summary_path, dpi=300)
plt.close()
print(f"✅ Summary dashboard generated.\nAll results saved in: {OUTPUT_DIR}")


