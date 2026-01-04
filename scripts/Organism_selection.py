# -*- coding: utf-8 -*-
"""
Created on Sun Sep 28 21:36:00 2025

@author: sofia
"""

import os
import pandas as pd
import requests
from io import StringIO

# Changer le working directory
os.chdir(r"C:\Users\sofia\Desktop\GeCoproject")  # adapte le chemin

# ----------------------------
# 1. Télécharger assembly_summary (RefSeq Bacteria)
# ----------------------------
URL_ASM = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"

print("Téléchargement assembly_summary...")
r = requests.get(URL_ASM)
r.raise_for_status()

# Lecture avec noms de colonnes NCBI (officiels)
colnames = [
    "assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category",
    "taxid", "species_taxid", "organism_name", "infraspecific_name", "isolate",
    "version_status", "assembly_level", "release_type", "genome_rep",
    "seq_rel_date", "asm_name", "asm_submitter", "gbrs_paired_asm", "paired_asm_comp",
    "ftp_path", "excluded_from_refseq", "relation_to_type_material",
    "asm_not_live_date", "assembly_type", "group", "genome_size", "genome_size_ungapped", "gc_percent", "replicon_count", "scaffold_count",
    "contig_count", "annotation_provider", "annotation_name", "annotation_date",
    "total_gene_count", "protein_coding_gene_count", "non_coding_gene_count", "pubmed_id"
]

asm_df = pd.read_csv(StringIO(r.text), sep="\t", comment="#", names=colnames, dtype=str)
print(f"Assemblies chargées : {len(asm_df)}")

# Convertir colonnes numériques
for col in ["total_gene_count", "protein_coding_gene_count"]:
    asm_df[col] = pd.to_numeric(asm_df[col], errors="coerce")

# ----------------------------
# 2. Charger le fichier COG local
# ----------------------------
COG_FILE = r"C:\Users\sofia\Desktop\GeCoproject\data\cog-24.cog.csv"  

cog_columns = [
    "locus_tag", "assembly_accession", "protein_id", "protein_length",
    "protein_coords", "protein_length_alt", "cog_id", "cog_id_alt",
    "unknown1", "score", "evalue", "alignment_length", "alignment_coords"
]

print("Chargement du fichier COG local...")
cog_df = pd.read_csv(COG_FILE, sep=",", header=None, names=cog_columns)
print(f"Lignes dans COG : {len(cog_df)}")

# ----------------------------
# 3. Compter le nombre de CDS par assembly
# ----------------------------
cds_counts = cog_df.groupby('assembly_accession').size()
print(f"Nombre d'assemblies uniques dans COG : {len(cds_counts)}")

# ----------------------------
# 4. Filtrer assembly_summary
# ----------------------------
asm_df_cog = asm_df[asm_df['assembly_accession'].isin(cds_counts.index)].copy()

filtered = asm_df_cog[
    (asm_df_cog["genome_rep"] == 'Full') &
    (asm_df_cog['assembly_accession'].isin(cds_counts[(cds_counts >= 1000) & (cds_counts <= 2000)].index))
].copy()

print(f"Génomes candidats après filtres : {len(filtered)}")

# ----------------------------
# 5. Supprimer les organismes déjà utilisés
# ----------------------------
organisms_used = [
    "Pediococcus claussenii ATCC BAA-344", "Streptococcus pyogènes M1 GAS",
    "Latilactobacillus sakei", "Streptococcus pneumoniae TIGR4",
    "Borrelia hermsii HS1", "Pediococcus pentosaceus", "Pediococcus acidilactici",
    "Limosilactobacillus mucosae", "Lactobacillus paragasseri", "Streptococcus mutans",
    "Streptococcus thermophilus", "Lactococcus lactis subsp. lactis", "Lactococcus garvieae",
    "Latilactobacillus graminis", "Latilactobacillus curvatus", "Lactobacillus acidophilus",
    "Fructobacillus americanaquae", "Streptococcus sanguinis", "Streptococcus salivarius",
    "Floricoccus penangensis", "Borrelia turcica", "Borrelia miyamotoi",
    "Borreliella mayonii", "Borreliella burgdorferi"
]
filtered = filtered[~filtered['organism_name'].isin(organisms_used)].copy()
filtered.reset_index(drop=True, inplace=True)

# ----------------------------
# 6. Choisir un organisme valide
# ----------------------------
organism_of_interest = None
same_genus_pool, other_genera_pool = None, None

for idx, row in filtered.sample(frac=1).iterrows():
    genus_interest = row['organism_name'].split()[0]

    # Autres du même genre
    same_genus_pool = filtered[
        (filtered['organism_name'].str.startswith(genus_interest)) &
        (filtered['organism_name'] != row['organism_name'])
    ]
    if len(same_genus_pool) < 2:
        continue

    # Genres différents
    other_genera_pool = filtered[~filtered['organism_name'].str.startswith(genus_interest)].copy()
    other_genera_pool['genus'] = other_genera_pool['organism_name'].apply(lambda x: x.split()[0])
    if other_genera_pool['genus'].nunique() < 2:
        continue

    organism_of_interest = row
    break

if organism_of_interest is None:
    raise ValueError("Aucun organisme ne remplit toutes les conditions.")

print(f"\n✅ Organisme d’intérêt choisi : {organism_of_interest['organism_name']}")

# ----------------------------
# 7. Sélection finale
# ----------------------------
same_genus = same_genus_pool.sample(n=2)
close_genera = other_genera_pool.groupby('genus').first().sample(n=2)
other_organisms = other_genera_pool[other_genera_pool['genus'].isin(close_genera.index)].sample(n=2)

final_selection = pd.concat([pd.DataFrame([organism_of_interest]), same_genus, other_organisms])
final_selection = final_selection[['organism_name', 'assembly_accession', 'ftp_path']]

print("\n✅ Sélection finale des 5 génomes :")
print(final_selection)

# ----------------------------
# 8. Sauvegarder
# ----------------------------
os.makedirs("results", exist_ok=True)  # crée results/ si besoin
out_path = os.path.join("results", "genomes_selected.tsv")

final_selection.to_csv(out_path, sep="\t", index=False)
print(f"\nTable sauvegardée dans '{out_path}'")
