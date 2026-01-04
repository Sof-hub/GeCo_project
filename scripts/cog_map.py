#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Colorisation COG pour Apilactobacillus kunkeei - VERSION FINALE
Compatible avec le format: COG_category=['J'] et COG: ['J']
"""

from Bio import SeqIO
import os

COG_COLORS = {
    'J': '#4169E1',  # Bleu royal - Translation, ribosomal structure
    'K': '#87CEEB',  # Bleu ciel - Transcription
    'L': '#32CD32',  # Vert lime - Replication, recombination, repair
    'D': '#9370DB',  # Violet moyen - Cell cycle control
    'O': '#FF69B4',  # Rose chaud - Post-translational modification
    'M': '#20B2AA',  # Turquoise - Cell wall/membrane biogenesis
    'N': '#DDA0DD',  # Prune - Cell motility
    'P': '#00CED1',  # Cyan foncé - Inorganic ion transport
    'T': '#FF00FF',  # Magenta - Signal transduction
    'C': '#FF4500',  # Rouge-orange - Energy production
    'G': '#FF1493',  # Rose profond - Carbohydrate metabolism
    'E': '#FFA500',  # Orange - Amino acid metabolism
    'F': '#FFD700',  # Or - Nucleotide metabolism
    'H': '#9370DB',  # Violet - Coenzyme metabolism
    'I': '#FF6347',  # Tomate - Lipid metabolism
    'Q': '#BA55D3',  # Orchidée - Secondary metabolites
    'R': '#808080',  # Gris - General function prediction
    'S': '#C0C0C0',  # Argent - Function unknown
    'U': '#696969',  # Gris foncé - Intracellular trafficking
    'V': '#8B0000',  # Rouge foncé - Defense mechanisms
    'W': '#8B4513',  # Brun - Extracellular structures
    'Y': '#2F4F4F',  # Gris ardoise - Nuclear structure
    'Z': '#000080',  # Bleu marine - Cytoskeleton
}

def blend_colors(color1, color2):
    """Mélange deux couleurs hex"""
    r1, g1, b1 = int(color1[1:3], 16), int(color1[3:5], 16), int(color1[5:7], 16)
    r2, g2, b2 = int(color2[1:3], 16), int(color2[3:5], 16), int(color2[5:7], 16)
    r, g, b = (r1+r2)//2, (g1+g2)//2, (b1+b2)//2
    return f"#{r:02x}{g:02x}{b:02x}"

def get_cog_color(cog_category):
    """Retourne la couleur (simple ou mélangée)"""
    if not cog_category:
        return '#808080'
    if len(cog_category) == 1:
        return COG_COLORS.get(cog_category, '#808080')
    # Multiple : mélanger les 2 premières couleurs
    colors = [COG_COLORS.get(c, '#808080') for c in cog_category if c in COG_COLORS]
    if len(colors) >= 2:
        return blend_colors(colors[0], colors[1])
    elif len(colors) == 1:
        return colors[0]
    return '#808080'

def extract_cog_category(feature):
    """
    Extrait la catégorie COG depuis les qualifiers
    Compatible avec COG_category=['J'] et COG: ['J']
    Gère aussi les multiples: ['EH'], ['JK'], etc.
    """
    # Méthode 1 : Chercher dans 'COG_category'
    cog_cat = feature.qualifiers.get('COG_category', [])
    if cog_cat:
        # Extraire les lettres de la première valeur
        cat_str = str(cog_cat[0]) if isinstance(cog_cat, list) else str(cog_cat)
        letters = ''.join([c for c in cat_str if c.isalpha()])
        if letters:
            return letters
    
    # Méthode 2 : Chercher dans 'COG'
    cog = feature.qualifiers.get('COG', [])
    if cog:
        cat_str = str(cog[0]) if isinstance(cog, list) else str(cog)
        letters = ''.join([c for c in cat_str if c.isalpha()])
        if letters:
            return letters
    
    return None

def colorize_apilactobacillus_kunkeei(input_gbff, output_gbff=None):
    """
    Colorisation des COG pour A. kunkeei
    """
    if output_gbff is None:
        base, ext = os.path.splitext(input_gbff)
        output_gbff = f"{base}_colored{ext}"
    
    print("="*70)
    print("🧬 COLORISATION COG - Apilactobacillus kunkeei")
    print("="*70)
    print(f"\n📂 Lecture : {input_gbff}")
    
    try:
        record = SeqIO.read(input_gbff, "genbank")
    except:
        records = list(SeqIO.parse(input_gbff, "genbank"))
        record = records[0]
        print("⚠️  Plusieurs records, utilisation du premier")
    
    colored_count = 0
    cog_stats = {}
    total_cds = 0
    
    # Traiter chaque CDS
    for feature in record.features:
        if feature.type == "CDS":
            total_cds += 1
            
            # Extraire la catégorie COG
            cog_cat = extract_cog_category(feature)
            
            if cog_cat:
                # Statistiques
                cog_stats[cog_cat] = cog_stats.get(cog_cat, 0) + 1
                
                # Obtenir la couleur
                color = get_cog_color(cog_cat)
                
                # Ajouter les qualifiers pour SnapGene
                feature.qualifiers["ApEinfo_fwdcolor"] = [color]
                feature.qualifiers["ApEinfo_revcolor"] = [color]
                feature.qualifiers["ApEinfo_graphicformat"] = ["arrow_data {{0 1 2 0 0 -1} {} 0}"]
                
                colored_count += 1
    
    # Sauvegarder
    SeqIO.write(record, output_gbff, "genbank")
    
    print(f"\n{'='*70}")
    print("📊 RÉSULTATS :")
    print(f"{'='*70}")
    print(f"  Total CDS               : {total_cds}")
    print(f"  CDS avec COG colorées   : {colored_count}")
    print(f"  CDS sans COG            : {total_cds - colored_count}")
    print(f"  Taux de couverture COG  : {100*colored_count/total_cds:.1f}%")
    print(f"\n💾 Fichier sauvegardé : {output_gbff}")
    
    if cog_stats:
        print(f"\n{'='*70}")
        print("🎨 CATÉGORIES COG COLORÉES :")
        print(f"{'='*70}")
        
        simple = {k: v for k, v in cog_stats.items() if len(k) == 1}
        multiple = {k: v for k, v in cog_stats.items() if len(k) > 1}
        
        if simple:
            print("\n🔹 Catégories simples :")
            for cat in sorted(simple.keys()):
                color = COG_COLORS.get(cat, '#808080')
                print(f"   {cat:3s} → {color} ({simple[cat]:4d} gènes)")
        
        if multiple:
            print("\n🔸 Catégories multiples (couleurs mélangées) :")
            for cat in sorted(multiple.keys()):
                color = get_cog_color(cat)
                components = ' + '.join([c for c in cat])
                print(f"   {cat:3s} ({components}) → {color} ({multiple[cat]:4d} gènes)")
    
    print(f"\n{'='*70}")
    print("✨ TERMINÉ !")
    print("   Ouvrez le fichier *_colored.gbff dans SnapGene")
    print(f"{'='*70}")
    
    return output_gbff

if __name__ == "__main__":
    # Chemin de votre fichier
    input_file = r"C:\Users\sofia\Desktop\GeCoproject\Apilactobacillus_kunkeei_COG.gbff"
    
    if os.path.exists(input_file):
        colorize_apilactobacillus_kunkeei(input_file)
    else:
        print(f"❌ Fichier non trouvé : {input_file}")

"""
Génération de la palette COG avec légende visuelle
Crée une image PNG montrant chaque catégorie COG avec sa couleur
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import numpy as np

# Palette de couleurs COG avec descriptions complètes
COG_LEGEND = {
    'J': ('#4169E1', 'Translation, ribosomal structure and biogenesis'),
    'K': ('#87CEEB', 'Transcription'),
    'L': ('#32CD32', 'Replication, recombination and repair'),
    'D': ('#9370DB', 'Cell cycle control, cell division'),
    'O': ('#FF69B4', 'Post-translational modification, chaperones'),
    'M': ('#20B2AA', 'Cell wall/membrane/envelope biogenesis'),
    'N': ('#DDA0DD', 'Cell motility'),
    'P': ('#00CED1', 'Inorganic ion transport and metabolism'),
    'T': ('#FF00FF', 'Signal transduction mechanisms'),
    'C': ('#FF4500', 'Energy production and conversion'),
    'G': ('#FF1493', 'Carbohydrate transport and metabolism'),
    'E': ('#FFA500', 'Amino acid transport and metabolism'),
    'F': ('#FFD700', 'Nucleotide transport and metabolism'),
    'H': ('#9370DB', 'Coenzyme transport and metabolism'),
    'I': ('#FF6347', 'Lipid transport and metabolism'),
    'Q': ('#BA55D3', 'Secondary metabolites biosynthesis'),
    'R': ('#808080', 'General function prediction only'),
    'S': ('#C0C0C0', 'Function unknown'),
    'U': ('#696969', 'Intracellular trafficking and secretion'),
    'V': ('#8B0000', 'Defense mechanisms'),
}

def create_cog_legend_horizontal(output_file='COG_palette_legend.png', dpi=300):
    """
    Crée une légende horizontale compacte avec carrés de couleur
    """
    fig, ax = plt.subplots(figsize=(16, 10))
    ax.axis('off')
    
    # Titre
    fig.suptitle('Palette de couleurs COG - Catégories fonctionnelles', 
                 fontsize=18, fontweight='bold', y=0.98)
    
    # Organisation en 2 colonnes
    categories = list(COG_LEGEND.keys())
    n_per_col = len(categories) // 2 + len(categories) % 2
    
    y_start = 0.92
    y_step = 0.042
    col_width = 0.45
    
    for i, cat in enumerate(categories):
        color, description = COG_LEGEND[cat]
        
        # Déterminer la colonne et position Y
        col = i // n_per_col
        row = i % n_per_col
        x_pos = 0.05 + col * (col_width + 0.05)
        y_pos = y_start - row * y_step
        
        # Carré de couleur
        rect = Rectangle((x_pos, y_pos), 0.025, 0.03, 
                        facecolor=color, edgecolor='black', linewidth=1.5,
                        transform=fig.transFigure)
        fig.patches.append(rect)
        
        # Texte: Catégorie + Description
        text = f"{cat}  -  {description}"
        plt.text(x_pos + 0.035, y_pos + 0.015, text,
                transform=fig.transFigure,
                fontsize=11, verticalalignment='center',
                fontweight='bold' if cat in ['J', 'K', 'L', 'E', 'C'] else 'normal')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✅ Légende horizontale sauvegardée : {output_file}")

def create_cog_legend_vertical(output_file='COG_palette_vertical.png', dpi=300):
    """
    Crée une légende verticale avec barres de couleur larges
    """
    categories = list(COG_LEGEND.keys())
    n_cats = len(categories)
    
    fig, ax = plt.subplots(figsize=(10, 12))
    
    # Titre
    ax.text(0.5, 0.98, 'Palette COG - Catégories fonctionnelles',
            ha='center', fontsize=16, fontweight='bold',
            transform=ax.transAxes)
    
    # Créer les barres
    y_positions = np.arange(n_cats)[::-1]  # Inverser pour avoir J en haut
    
    for i, (cat, (color, desc)) in enumerate(COG_LEGEND.items()):
        y = y_positions[i]
        
        # Barre de couleur
        ax.barh(y, 0.5, height=0.8, left=0, color=color, 
                edgecolor='black', linewidth=1.5)
        
        # Lettre COG (sur la barre)
        ax.text(0.25, y, cat, ha='center', va='center',
                fontsize=14, fontweight='bold', color='white')
        
        # Description (à droite)
        ax.text(0.6, y, desc, ha='left', va='center',
                fontsize=10, fontweight='normal')
    
    ax.set_xlim(0, 5)
    ax.set_ylim(-1, n_cats)
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✅ Légende verticale sauvegardée : {output_file}")

def create_cog_palette_grid(output_file='COG_palette_grid.png', dpi=300):
    """
    Crée une grille compacte style mosaïque avec code couleur
    """
    fig, axes = plt.subplots(5, 4, figsize=(12, 10))
    fig.suptitle('Palette de couleurs COG', fontsize=16, fontweight='bold')
    
    axes = axes.flatten()
    
    for i, (cat, (color, desc)) in enumerate(COG_LEGEND.items()):
        ax = axes[i]
        
        # Rectangle coloré
        ax.add_patch(Rectangle((0, 0), 1, 1, facecolor=color, edgecolor='black', linewidth=2))
        
        # Catégorie (grande lettre)
        ax.text(0.5, 0.65, cat, ha='center', va='center',
                fontsize=32, fontweight='bold', color='white')
        
        # Couleur hex
        ax.text(0.5, 0.35, color, ha='center', va='center',
                fontsize=9, color='white', fontweight='bold')
        
        # Description (en bas)
        ax.text(0.5, 0.05, desc, ha='center', va='bottom',
                fontsize=7, color='white', wrap=True)
        
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
    
    # Masquer les axes non utilisés
    for j in range(i+1, len(axes)):
        axes[j].axis('off')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✅ Grille palette sauvegardée : {output_file}")

def create_cog_color_table(output_file='COG_color_table.txt'):
    """
    Crée un fichier texte avec le tableau des couleurs
    Format compatible pour LaTeX, Markdown, etc.
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("PALETTE DE COULEURS COG - TABLE DE RÉFÉRENCE\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"{'Cat':<5} {'Couleur':<10} {'Description':<60}\n")
        f.write("-"*80 + "\n")
        
        for cat, (color, desc) in COG_LEGEND.items():
            f.write(f"{cat:<5} {color:<10} {desc:<60}\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("\nFORMAT LATEX (pour tableau) :\n")
        f.write("-"*80 + "\n")
        
        for cat, (color, desc) in COG_LEGEND.items():
            # Convertir hex en format LaTeX RGB
            r = int(color[1:3], 16) / 255
            g = int(color[3:5], 16) / 255
            b = int(color[5:7], 16) / 255
            
            f.write(f"{cat} & \\cellcolor[RGB]{{{int(r*255)},{int(g*255)},{int(b*255)}}} {color} & {desc} \\\\\n")
    
    print(f"✅ Table texte sauvegardée : {output_file}")

def create_all_palettes():
    """
    Génère toutes les versions de la palette
    """
    print("\n" + "="*70)
    print("🎨 GÉNÉRATION DES PALETTES COG")
    print("="*70 + "\n")
    
    create_cog_legend_horizontal('COG_palette_horizontal.png', dpi=300)
    create_cog_legend_vertical('COG_palette_vertical.png', dpi=300)
    create_cog_palette_grid('COG_palette_grid.png', dpi=300)
    create_cog_color_table('COG_color_table.txt')
    
    print("\n" + "="*70)
    print("✨ TOUTES LES PALETTES GÉNÉRÉES !")
    print("="*70)
    print("\n📁 Fichiers créés :")
    print("  - COG_palette_horizontal.png  (légende compacte)")
    print("  - COG_palette_vertical.png    (barres colorées)")
    print("  - COG_palette_grid.png        (grille mosaïque)")
    print("  - COG_color_table.txt         (table texte/LaTeX)")
    print("\n💡 Vous pouvez insérer ces images dans votre rapport ou présentation")
    print("="*70 + "\n")

if __name__ == "__main__":
    import os
    
    # Changer vers le dossier de travail
    os.chdir(r"C:\Users\sofia\Desktop\GeCoproject")
    
    # Générer toutes les palettes
    create_all_palettes()

