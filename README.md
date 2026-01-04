# Comparative Genomics of Bee-Associated Lactobacillaceae

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![BioPython](https://img.shields.io/badge/BioPython-1.79+-green.svg)](https://biopython.org/)

> **Genomic signatures, replication organization, and functional annotation in five bee-associated _Lactobacillaceae_ species**

## 📋 Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Data](#data)
- [Usage](#usage)
- [Analyses Performed](#analyses-performed)
- [Results Interpretation](#results-interpretation)
- [Methodology](#methodology)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

---

## 🔬 Overview

This repository contains the complete analysis pipeline, data, and results for a comparative genomic study of five _Lactobacillaceae_ species associated with honeybees:
- **_Apilactobacillus kunkeei_** (GCF_019575995.1)
- **_Apilactobacillus bombintestini_** (GCF_003627035.1)
- **_Apilactobacillus apinorum_** (GCF_946888465.1)
- **_Bombilactobacillus folatiphilus_** (GCF_023380265.1)
- **_Bombilactobacillus thymidiniphilus_** (GCF_023380245.1)

### Key Findings

✅ **Replication Architecture**: Clear bidirectional replication with _ori_/_ter_ identified by GC skew  
✅ **Genomic Signatures**: AT-rich composition (GC 34-38%), strong CG depletion  
✅ **Codon Usage**: AT-rich codon preference, CAI variation between genera  
✅ **Synteny**: Conserved within _Apilactobacillus_, rearranged between genera  
✅ **Pan-genome**: Open structure with core + accessory + species-specific genes  
✅ **Phylogeny**: Three _Apilactobacillus_ form tight clade; _Bombilactobacillus_ more divergent  

---

## 📁 Repository Structure

```
Lactobacillaceae-Comparative-Genomics/
├── README.md # This file
├── LICENSE # MIT License
├── requirements.txt # Python dependencies
├── CHANGELOG.md # Version history
│
├── data/ # Input genomic data
│ ├── 1_GCF_019575995.1/ # A. kunkeei
│ ├── 2_GCF_003627035.1/ # A. bombintestini
│ ├── 3_GCF_946888465.1/ # A. apinorum
│ ├── 4_GCF_023380265.1/ # B. folatiphilus
│ └── 5_GCF_023380245.1/ # B. thymidiniphilus
│ ├── *.fna # Genomic sequence
│ ├── *.faa # Proteins
│ ├── *.gbff # GenBank
│ ├── *_cusp.txt # Codon usage
│ └── *.emapper.annotations.xlsx # COG annotations
│
├── scripts/ # Analysis scripts
│ ├── Organism_selection.py # Genome selection
│ ├── GeCo_analysis_v2.py # Main pipeline
│ └── Mauve_Newick_tree.py # Phylogenetic tree
│
└── results/ # Output files
├── tables/ # TSV results
└── figures/ # PNG figures (300 dpi)

```

---

## 🛠️ Installation

### Prerequisites

- Python 3.8+
- BioPython 1.79+
- pandas, numpy, matplotlib, seaborn

### Quick Setup

```bash
# 1. Clone repository
git clone https://github.com/yourusername/Lactobacillaceae-Comparative-Genomics.git
cd Lactobacillaceae-Comparative-Genomics

# 2. Create virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt
```
requirements.txt:
```
biopython>=1.79
pandas>=1.3.0
numpy>=1.21.0
matplotlib>=3.4.0
seaborn>=0.11.0
openpyxl>=3.0.0
requests>=2.26.0
```
📊 Data
Genome Overview

| Species            | Accession       | Size (Mb) | GC%  | CDS  |
| ------------------ | --------------- | --------- | ---- | ---- |
| A. kunkeei         | GCF_019575995.1 | 1.54      | 36.6 | 1343 |
| A. bombintestini   | GCF_003627035.1 | 1.29      | 34.1 | 1159 |
| A. apinorum        | GCF_946888465.1 | 1.46      | 34.7 | 1329 |
| B. folatiphilus    | GCF_023380265.1 | 1.64      | 38.4 | 1506 |
| B. thymidiniphilus | GCF_023380245.1 | 1.49      | 36.4 | 1412 |

Required Files per Genome
.fna - Genomic DNA sequence (FASTA)

.faa - Protein sequences (FASTA)

.gbff - GenBank file (for CDS extraction)

*_cusp.txt - EMBOSS Cusp output

*.emapper.annotations.xlsx - eggNOG-mapper annotations

External Tools for Data Preparation
```bash
# Download genome from NCBI
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/[...]/genome.fna.gz
gunzip genome.fna.gz

# Functional annotation with eggNOG-mapper
emapper.py -i protein.faa --output_dir annotations -o genome_name --excel

# Codon usage with EMBOSS
cusp -sequence genome.fna -outfile genome_cusp.txt
```
🚀 Usage
Quick Start
Run complete analysis pipeline:

```bash
cd scripts/
python GeCo_analysis_v2.py
```
Runtime: ~5-10 minutes
Output: All results in results/ directory

Configuration
Edit file paths in GeCo_analysis_v2.py:
```python
genomes = {
    "Apilactobacillus kunkeei": {
        "fna": r"path/to/genome.fna",
        "faa": r"path/to/protein.faa",
        "cusp": r"path/to/cusp.txt",
        "emapper": r"path/to/annotations.xlsx",
        "genbank": r"path/to/genome.gbff"
    },
    # ... (repeat for 5 species)
}
```

Individual Scripts
1. Genome Selection
```bash
python Organism_selection.py
```
Output: genomes_selected.tsv

2. Main Analysis Pipeline
```bash
python GeCo_analysis_v2.py
```
Output: TSV tables + PNG figures

3. Phylogenetic Tree
```bash
python Mauve_Newick_tree.py
```
Output: Phylogenetic_tree_Mauve.png

🧬 Analyses Performed
1. General Genome Features
Genome size, GC content, CDS count

Coding density, tRNA/rRNA counts

2. GC Skew Analysis
Local GC skew: (G-C)/(G+C)

Cumulative GC skew → ori/ter detection

Replication symmetry index

3. Codon Usage
RSCU: Relative Synonymous Codon Usage

CAI: Codon Adaptation Index (reference: A. kunkeei)

4. Dinucleotide Frequencies
Observed vs. expected frequencies

Obs/Exp ratios for compositional bias

5. Functional Annotation
COG category distribution

NCBI COG vs eggNOG-mapper comparison

6. Synteny Analysis
D-GENIES dot-plots (pairwise alignments)

Mauve LCBs (Locally Collinear Blocks)

7. Pan-genome
Core vs accessory vs species-specific genes

Venn diagram + UpSet plot

8. Phylogenetic Analysis
Mauve-based tree (genome alignment)

OrthoVenn-based tree (protein orthologs)

📈 Results Interpretation
Key Output Files
Tables (TSV)

| File                       | Description                              |
| -------------------------- | ---------------------------------------- |
| General_features.tsv       | Genome stats (size, GC%, CDS)            |
| Replication_sites.tsv      | Ori/ter positions, symmetry              |
| CAI_values.tsv             | Codon adaptation index                   |
| Dinucleotide_frequency.tsv | All 16 dinucleotide frequencies          |
| RSCU_all_organisms.tsv     | Codon usage bias (64 codons × 5 species) |
| COG_all_organisms.tsv      | COG category percentages                 |

Figures (PNG, 300 dpi)
Replication:

Cumulative_GC_skew_with_labels.png

Replication_symmetry_summary.png

Signatures:

Dinucleotide_frequencies_heatmap.png

RSCU_heatmap.png

Synteny:

DGENIES_1vs2-5.png (4 dot-plots)

Alignement_1vs2-5.jpg (4 Mauve alignments)

Functional:

COG_all_comparison_grouped.png

COG_Apilactobacillus_kunkeei_comparison.png

Pan-genome:

jVenn_chart.png

UpSetJS.png

Phylogeny:

Phylogenetic_tree_Mauve.png

orthovenn_fasta_tree.png

Summary:

Summary_comparative_dashboard.png

How to Interpret Results
1. GC Skew Profiles
V-shaped cumulative curve = Normal replication

Minimum → Replication origin (ori)

Maximum → Replication terminus (ter)

Symmetry ~50% = Balanced replichores

Asymmetry > 50% = Chromosomal rearrangements

Example:
```
A. bombintestini: 66.2% symmetry
→ Suggests unequal replichore lengths
→ May indicate deletion/insertion events
```
2. Dinucleotide Obs/Exp Ratios

| Dinucleotide | Typical Ratio | Meaning                     |
| ------------ | ------------- | --------------------------- |
| CG           | < 0.5         | CpG depletion (methylation) |
| TA           | > 1.0         | AT-rich signature           |
| AA, TT       | > 1.0         | AT enrichment               |

3. RSCU (Codon Usage)
```
RSCU > 1.5 → Strong preference
RSCU ≈ 1.0 → Neutral
RSCU < 0.5 → Strong avoidance
```
AT-rich genomes:

TTA (Leucine): RSCU often > 2.0 ✓

CTG (Leucine): RSCU often < 0.2 ✓

4. CAI Values
```
Reference (A. kunkeei) = 1.000
CAI > 1.05 → Better translational optimization
CAI < 0.95 → Lower optimization
```
Results:

Bombilactobacillus (1.17-1.20) > Apilactobacillus (0.94-1.00)

Suggests better codon adaptation in Bombilactobacillus

5. Synteny (Dot-plots)

| Pattern             | Interpretation     |
| ------------------- | ------------------ |
| Continuous diagonal | Conserved synteny  |
| Diagonal break      | Inversion          |
| Offset diagonal     | Translocation      |
| Gaps                | Insertion/deletion |

Observed:

Apilactobacillus intra-genus: Conserved

Inter-genus: Multiple rearrangements

6. Pan-genome Structure
```
Core genome: Shared by all 5 species (essential functions)
Accessory: Present in 2-4 species (adaptation)
Species-specific: Unique to one species (niche-specific)
```
COG enrichment:

Core: J (translation), K (transcription), L (replication)

Accessory: G (carbohydrate), P (transport), V (defense)

7. Phylogenetic Trees
Branch lengths:

Apilactobacillus: 0.015-0.247 (recent radiation)

Bombilactobacillus: 0.309-0.327 (greater divergence)

Validation:

Mauve tree ≈ OrthoVenn tree → Robust signal ✓

🔬 Methodology
GC Skew Calculation
Formula:
```python
GC_skew = (G_count - C_count) / (G_count + C_count)
```
Parameters:

Window: 10,000 bp

Step: 5,000 bp

Cumulative sum identifies ori (min) and ter (max)

Dinucleotide Analysis
```python
# Observed
freq_obs[XY] = count[XY] / total_dinucleotides

# Expected
freq_exp[XY] = freq[X] × freq[Y]

# Ratio
ratio = freq_obs / freq_exp
```

RSCU Calculation
```python
for amino_acid in genetic_code:
    synonymous_codons = get_synonymous(amino_acid)
    n = len(synonymous_codons)
    total = sum(count[codon] for codon in synonymous_codons)
    
    for codon in synonymous_codons:
        expected = total / n
        RSCU[codon] = count[codon] / expected
```
(suite du README.md)

### CAI Calculation

**Reference set:** _A. kunkeei_ CDS sequences

```python
# Step 1: Relative adaptiveness
for amino_acid in genetic_code:
    synonymous_codons = get_synonymous(amino_acid)
    max_freq = max(ref_freq[codon] for codon in synonymous_codons)
    
    for codon in synonymous_codons:
        w[codon] = ref_freq[codon] / max_freq

# Step 2: CAI for each genome
CAI = exp( (1/L) × sum(ln(w[codon_i])) )
# where L = total codons in all CDS
```
Replication Symmetry
```python
ori_ter_distance = abs(ter_position - ori_position)
midpoint = (ori_position + ter_position) / 2
genome_center = genome_size / 2

symmetry_index = 100 × (1 - abs(midpoint - genome_center) / genome_center)
```
Interpretation:

50% = Perfect symmetry

50% = One replichore longer than the other

🐛 Troubleshooting
Common Issues
1. File Path Errors
Problem: FileNotFoundError: [Errno 2] No such file or directory

Solution:
```python
# Use raw strings for Windows paths
path = r"C:\data\genome.fna"

# Or use forward slashes (works on all OS)
path = "C:/data/genome.fna"

# Or use pathlib (recommended)
from pathlib import Path
path = Path("data") / "genome.fna"
```

2. Memory Errors
Problem: MemoryError with large genomes

Solutions:

Reduce window size: window=5000 instead of 10000

Process genomes one at a time

Increase system RAM

Use chunked processing

3. Missing Dependencies
Problem: ModuleNotFoundError: No module named 'Bio'

Solution:
```bash
pip install biopython pandas matplotlib seaborn openpyxl
```

4. EMBOSS Cusp Format
Problem: Cannot extract coding GC% from cusp file

Solution: Ensure cusp output contains:
```
Coding GC 36.60%
```
Run cusp correctly:
```bash
cusp -sequence genome.fna -outfile genome_cusp.txt
```
5. eggNOG File Issues
Problem: KeyError: 'COG_category'

Solution:

Verify Excel file has column COG_category

Re-run eggNOG-mapper with --excel flag

Check file encoding (UTF-8)

6. BioPython GenBank Parsing
Problem: No CDS extracted from .gbff

Solution:

Verify GenBank file is not corrupted

Check that translation qualifier exists in CDS features

Ensure CDS sequences are multiples of 3

7. Matplotlib Display Issues
Problem: Figures not displaying

Solution:
```python
# For non-GUI environments
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
```
Getting Help
GitHub Issues: Open an issue

Email: Contact author (see below)

Error Reports: Include:

Python version (python --version)

Error message (full traceback)

Operating system

Input file formats

📖 Advanced Usage
Custom Genome Analysis
To analyze your own genomes:
```python
# Edit GeCo_analysis_v2.py
genomes = {
    "Your_Species_Name": {
        "fna": "path/to/genome.fna",
        "faa": "path/to/proteins.faa",
        "cusp": "path/to/cusp_output.txt",
        "emapper": "path/to/emapper.xlsx",
        "genbank": "path/to/genome.gbff"
    },
    # Add more species...
}
```
Batch Processing
Process multiple genome sets:
```bash
#!/bin/bash
# batch_analysis.sh

for dataset in dataset1 dataset2 dataset3; do
    echo "Processing $dataset..."
    python scripts/GeCo_analysis_v2.py --config configs/${dataset}.yaml
    mv results/ results_${dataset}/
done
```
Integration with Other Tools
Export to R for statistical analysis
```python
import pandas as pd

# Load results
df = pd.read_csv('results/tables/General_features.tsv', sep='\t')

# Export for R
df.to_csv('for_R_analysis.csv', index=False)
```
```r
# In R
data <- read.csv("for_R_analysis.csv")
cor.test(data$Genome_size_bp, data$GC_coding)
```

🔄 Reproducibility Checklist
 All dependencies installed (see requirements.txt)

 Python version 3.8+ verified

 Input files in correct format (.fna, .faa, .gbff, etc.)

 File paths updated in scripts

 EMBOSS cusp run with correct parameters

 eggNOG-mapper run with --excel flag

 Results directory created

 Scripts executed in order (selection → analysis → phylogeny)

 Output files generated without errors

 Figures display correctly

🎓 Educational Use
This repository is ideal for teaching:

Graduate Courses:
Comparative genomics

Computational biology

Microbial genomics

Bioinformatics methods

Learning Objectives:
Understand replication-associated compositional biases

Interpret codon usage patterns

Analyze synteny and chromosomal rearrangements

Construct and validate phylogenetic trees

Explore pan-genome structure

Hands-on Exercises:
Modify window sizes for GC skew analysis

Compare different CAI reference organisms

Add additional genomes to the analysis

Implement new visualization methods

Statistical testing of observed patterns

📚 Citation
If you use this repository in your research or teaching, please cite:

This Work:
```
@mastersthesis{Quinteros2026,
  author       = {Sofia Quinteros},
  title        = {Genomic Signatures, Replication Organization, and Functional 
                  Annotation in Five Bee-Associated Lactobacillaceae Species},
  school       = {Aix-Marseille Université},
  year         = {2026},
  type         = {Master's Thesis},
  note         = {Master 2 Structural Biology and Genomics, 
                  Comparative Genomics Course (GeCo 2025)},
  url          = {https://github.com/yourusername/Lactobacillaceae-Comparative-Genomics}
}
```
Key References:
Bee Microbiome:

Ellegaard & Engel (2019). Genomic diversity landscape of the honey bee gut microbiota. Nat Commun 10:446.

Bradford et al. (2022). Comparative genomics of lactobacillaceae from honeybees. Environ Microbiol 24:5841-5854.

Genome Analysis Methods:

Lobry & Sueoka (2002). Asymmetric directional mutation pressures in bacteria. Genome Biol 3:research0058.

Sharp & Li (1987). The codon adaptation index. Nucleic Acids Res 15:1281-1295.

Software Tools:

Cantalapiedra et al. (2021). eggNOG-mapper v2. Mol Biol Evol 38:5825-5829.

Cabanettes & Klopp (2018). D-GENIES: dot plot large genomes. PeerJ 6:e4958.

Darling et al. (2010). progressiveMauve. PLoS ONE 5:e11147.

Xu et al. (2019). OrthoVenn2. Nucleic Acids Res 47:W52-W58.

🤝 Contributing
Contributions are welcome! We especially appreciate:

Types of Contributions:
🐛 Bug reports

💡 Feature requests

📝 Documentation improvements

🧪 Additional test cases

🎨 Visualization enhancements

🌍 Translations

How to Contribute:
1. Fork the repository
2. Create a feature branch

```bash
git checkout -b feature/YourAmazingFeature
```
3. Make your changes
4. Add tests if applicable
5. Commit with clear message

```bash
git commit -m "Add feature: YourAmazingFeature"
```
6. Push to branch

```bash
git push origin feature/YourAmazingFeature
```
7. Open a Pull Request

Code Style:
- Follow PEP 8 for Python
- Use descriptive variable names
- Add docstrings to functions
- Comment complex logic
- Keep functions < 50 lines when possible

Testing:
```bash
# Run tests before submitting PR
pytest tests/

# Check code style
flake8 scripts/

# Format code
black scripts/
```

👤 Contact & Support
Author
Sofia Quinteros
Master 2 Structural Biology and Genomics
Parcours: Génomique et Analyse Bio-informatique de Données
Aix-Marseille Université, France

📧 Email: sofia.quinteros@etu.univ-amu.fr
🔗 GitHub: @sofia-quinteros
🔗 LinkedIn: Sofia Quinteros

Supervisor
Dr. Emmanuel Talla
Course: Comparative Genomics (GeCo 2025)
Aix-Marseille Université

Support
💬 Questions: Open a Discussion

🐛 Bug Reports: Open an Issue

📧 Private Inquiries: Email author

📖 Documentation: See docs/ folder

🙏 Acknowledgments
This work was conducted as part of the Master 2 Structural Biology and Genomics program, Comparative Genomics Course (GeCo 2025), under the supervision of Dr. Emmanuel Talla.

Special Thanks:
NCBI RefSeq - High-quality genome assemblies

eggNOG Database - Functional annotations

EMBOSS Team - Codon usage analysis tools

Mauve Developers - Genome alignment software

D-GENIES Team - Interactive dot-plot visualization

OrthoVenn Team - Ortholog clustering

BioPython Community - Sequence analysis framework

Aix-Marseille Université - Academic support

Data Sources:
NCBI RefSeq: https://www.ncbi.nlm.nih.gov/refseq/

eggNOG 5.0: http://eggnog5.embl.de/

COG Database: https://www.ncbi.nlm.nih.gov/research/cog

🗺️ Roadmap
Version 1.1 (Q2 2026)
 Add interactive Plotly visualizations

 Implement Snakemake workflow

 Add unit tests (pytest)

 Create Docker container

 Add more example datasets

Version 2.0 (Q3 2026)
 Web interface (Flask/Dash)

 Real-time analysis monitoring

 Metabolic pathway reconstruction

 Machine learning predictions

 Comparative metagenomics module

Future Enhancements
 Support for >5 genomes

 Automated statistical testing

 Integration with Galaxy

 Cloud computing support (AWS, GCP)

 Mobile-responsive visualization

🔗 Related Resources
Similar Projects:
Bee Microbiome Database

Lactobacillus Genomics

[Comparative Genomics Pipeline](https://github.com/









