# -*- coding: utf-8 -*-
"""
Created on Sun Jan  4 10:37:00 2026

@author: sofia

Generation d'un arbre phylogénétique à partir du guide tree de Mauve
"""

import matplotlib.pyplot as plt
from Bio import Phylo
from io import StringIO
import matplotlib.patches as mpatches
import os

# Changer le working directory
os.chdir(r"C:\Users\sofia\Desktop\GeCoproject\results")

# Données de l'arbre (format Newick du fichier Mauve)
newick_string = "(seq5:0.30925,(seq4:0.326779:0.098448,(seq2:0.247818,seq3:0.218232,seq1:0.234038):0.0149703):0.0149703);"

# Mapping des noms de séquences vers les noms réels
organism_names = {
    'seq1': 'A. kunkeei',
    'seq2': 'A. bombintestini',
    'seq3': 'A. apinorum',
    'seq4': 'B. folatiphilus',
    'seq5': 'B. thymidiniphilus'
}

# Couleurs pour les genres
genus_colors = {
    'A.': '#1f77b4',  # Bleu pour Apilactobacillus
    'B.': '#ff7f0e'   # Orange pour Bombilactobacillus
}

def create_phylogenetic_tree():
    """Crée et visualise l'arbre phylogénétique"""
    
    # Parser l'arbre Newick
    tree = Phylo.read(StringIO(newick_string), "newick")
    
    # Renommer les terminaux avec les noms complets
    for terminal in tree.get_terminals():
        if terminal.name in organism_names:
            terminal.name = organism_names[terminal.name]
    
    # Créer la figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Dessiner l'arbre
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)
    
    # Personnaliser l'apparence
    ax.set_xlabel('Distance évolutive', fontsize=12, fontweight='bold')
    ax.set_title('Arbre phylogénétique basé sur les alignements génomiques (Mauve)\nCinq espèces de Lactobacillaceae associées aux abeilles', 
                fontsize=14, fontweight='bold', pad=20)
    
    # Ajouter la légende des genres
    legend_elements = [
        mpatches.Patch(facecolor=genus_colors['A.'], label='Apilactobacillus', alpha=0.7),
        mpatches.Patch(facecolor=genus_colors['B.'], label='Bombilactobacillus', alpha=0.7)
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=11, framealpha=0.9)
    
    # Ajouter une note sur la méthode
    note_text = "Arbre construit par Mauve progressiveMauve à partir des alignements de génomes complets.\nLes longueurs de branches représentent les distances évolutives estimées."
    ax.text(0.02, 0.02, note_text, transform=ax.transAxes, 
           fontsize=9, verticalalignment='bottom',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig('Phylogenetic_tree_Mauve.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✓ Arbre phylogénétique sauvegardé: Phylogenetic_tree_Mauve.png")


def create_simple_tree_matplotlib():
    """Crée un arbre phylogénétique simple avec matplotlib pur (alternative)"""
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Positions verticales pour chaque organisme
    y_positions = {
        'A. kunkeei': 5,
        'A. bombintestini': 4,
        'A. apinorum': 3,
        'B. folatiphilus': 2,
        'B. thymidiniphilus': 1
    }
    
    # Dessiner les branches principales
    # Groupe Apilactobacillus
    ax.plot([0, 0.234], [5, 5], 'k-', linewidth=2)  # A. kunkeei
    ax.plot([0, 0.247], [4, 4], 'k-', linewidth=2)  # A. bombintestini
    ax.plot([0, 0.218], [3, 3], 'k-', linewidth=2)  # A. apinorum
    
    # Nœud interne Apilactobacillus
    ax.plot([0.015, 0.015], [3, 5], 'k-', linewidth=2)
    
    # Groupe Bombilactobacillus
    ax.plot([0, 0.327], [2, 2], 'k-', linewidth=2)  # B. folatiphilus
    ax.plot([0.098, 0.098], [2, 4], 'k-', linewidth=2)
    ax.plot([0.098, 0.015], [4, 4], 'k-', linewidth=2)
    
    # B. thymidiniphilus
    ax.plot([0, 0.309], [1, 1], 'k-', linewidth=2)
    ax.plot([0.015, 0.015], [1, 2], 'k-', linewidth=2)
    
    # Ajouter les noms des organismes
    for org, y in y_positions.items():
        genus = org.split('.')[0]
        color = genus_colors[genus + '.']
        ax.text(0.35, y, f'  {org}', va='center', fontsize=11, 
               fontweight='bold', color=color)
        # Ajouter des points aux terminaisons
        distance = {'A. kunkeei': 0.234, 'A. bombintestini': 0.247, 
                   'A. apinorum': 0.218, 'B. folatiphilus': 0.327, 
                   'B. thymidiniphilus': 0.309}
        ax.plot(distance[org], y, 'o', color=color, markersize=8)
    
    ax.set_xlim(-0.02, 0.6)
    ax.set_ylim(0, 6)
    ax.set_xlabel('Distance évolutive', fontsize=12, fontweight='bold')
    ax.set_yticks([])
    ax.set_title('Arbre phylogénétique basé sur les alignements génomiques complets\n(Guide tree Mauve progressiveMauve)', 
                fontsize=14, fontweight='bold', pad=20)
    
    # Légende
    legend_elements = [
        mpatches.Patch(facecolor=genus_colors['A.'], label='Genre Apilactobacillus', alpha=0.7),
        mpatches.Patch(facecolor=genus_colors['B.'], label='Genre Bombilactobacillus', alpha=0.7)
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    # Note méthodologique
    note_text = ("Arbre phylogénétique construit par Mauve progressiveMauve\n"
                "à partir d'alignements de génomes complets.\n"
                "Les distances de branches représentent les substitutions\n"
                "par site aligné.")
    ax.text(0.02, 0.98, note_text, transform=ax.transAxes, 
           fontsize=9, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig('Phylogenetic_tree_Mauve_simple.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("✓ Arbre phylogénétique simple sauvegardé: Phylogenetic_tree_Mauve_simple.png")


def print_tree_info():
    """Affiche les informations sur l'arbre"""
    
    print("\n" + "="*60)
    print("INFORMATIONS SUR L'ARBRE PHYLOGÉNÉTIQUE")
    print("="*60 + "\n")
    
    print("Source: Guide tree de Mauve progressiveMauve")
    print("Méthode: Alignement de génomes complets\n")
    
    print("Distances évolutives (longueurs de branches):")
    print("  - A. kunkeei: 0.234")
    print("  - A. bombintestini: 0.247")
    print("  - A. apinorum: 0.218")
    print("  - B. folatiphilus: 0.327")
    print("  - B. thymidiniphilus: 0.309\n")
    
    print("Observations:")
    print("  1. Les trois Apilactobacillus forment un clade monophylétique")
    print("  2. Les deux Bombilactobacillus sont plus divergents")
    print("  3. B. thymidiniphilus apparaît comme groupe externe")
    print("  4. Les distances intra-Apilactobacillus sont courtes (0.015)")
    print("     suggérant une radiation récente\n")


if __name__ == "__main__":
    print("Génération de l'arbre phylogénétique à partir des données Mauve...\n")
    
    # Méthode 1: Utiliser BioPython (nécessite: pip install biopython)
    try:
        create_phylogenetic_tree()
    except ImportError:
        print("⚠ BioPython non installé. Utilisation de la méthode alternative...")
    
    # Méthode 2: Arbre simple avec matplotlib pur
    create_simple_tree_matplotlib()
    
    # Afficher les informations
    print_tree_info()
    
    print("\n" + "="*60)
    print("TERMINÉ !")
    print("="*60)
