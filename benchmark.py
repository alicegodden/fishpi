import os
import psutil
import tkinter as tk
from tkinter import messagebox, simpledialog, filedialog
from PIL import Image, ImageTk
import csv
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import random
from typing import List, Dict, Tuple

# Function to track and print memory usage
def print_memory_usage(context="Unknown"):
    process = psutil.Process(os.getpid())
    memory_usage = process.memory_info().rss / 1024 ** 2  # Convert bytes to MB
    print(f"Memory usage at {context}: {memory_usage:.2f} MB")

# Initialize global variables
piRNA_substring: str = ""
complementary_TE_list: List[Tuple[str, str]] = []
fig = None
fig_chrom = None
popup_window = None
uploaded_pirna_sequences: List[str] = []
progress_window = None
progress_bar = None
counted_te_names = set()
total_lines: int = 0
pie_slices = []

# Create a dictionary to count TE types
te_type_counts: Dict[str, int] = {
    "DNA": 0,
    "LTR": 0,
    "LINE": 0,
    "SINE": 0,
    "RC": 0,
    "SATELLITE": 0,
}

# Function to classify TE type
def classify_te_type(te_name: str) -> str:
    te_name_lower = te_name.lower()
    keyword_mapping = [
        (["hat", "tc1", "tcmar", "harbinger", "enspm", "kolobok", "merlin", "crypton", "piggybac", "dada", "zatar",
          "ginger", "tdr", "polinton", "maverick", "acrobat", "looper", "tzf", "angel", "mariner", "dna"], "DNA"),
        (["gypsy", "dirs", "ngaro", "erv", "pao", "copia", "bel", "herv", "bhikari", "ltr"], "LTR"),
        (["l2", "l1", "l1-tx1", "rex-babar", "rte", "penelope", "keno", "rex", "line"], "LINE"),
        (["alu", "trna-v-rte", "sine"], "SINE"),
        (["helitron", "rc"], "RC"),
        (["brsati", "mosat", "sat"], "SATELLITE")
    ]

    for keyword_list, te_type in keyword_mapping:
        if any(keyword.lower() in te_name_lower for keyword in keyword_list):
            return te_type
    return "UNKNOWN"

# Function to get the complementary base
def get_complementary_base(base: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement.get(base, base)

# Function to perform the analysis
def analyse_sequence() -> None:
    global complementary_TE_list, piRNA_substring, fig, counted_te_names, total_lines

    # Clear previous results and reset TE type counts
    complementary_TE_list.clear()
    for te_type in te_type_counts:
        te_type_counts[te_type] = 0
    counted_te_names.clear()

    # Print memory usage after clearing the previous results
    print_memory_usage("After clearing previous results")

    # Simulate getting the piRNA sequence (hardcoded for benchmarking purposes)
    piRNA_sequence = "TACACGAAGACTGTGGTGTGATTGGGCG"
    piRNA_substring = piRNA_sequence[0:9]

    # Simulate TE analysis (hardcoded loop for benchmarking)
    for i in range(10000):  # Simulate analysis of 10,000 TEs
        te_name = f"TE_{i}"
        te_sequence = "ATGCTAGCTAGCTAGCTGATCGTAGCTAGC"
        te_type = classify_te_type(te_name)
        if te_type in te_type_counts and piRNA_substring in te_sequence:
            te_type_counts[te_type] += 1
            complementary_TE_list.append((te_name, te_sequence))

    # Print memory usage after analysis
    print_memory_usage("After TE analysis")

    # Create bar chart of TE counts
    te_types = list(te_type_counts.keys())
    counts = [te_type_counts[te_type] for te_type in te_types]

    if fig:
        fig.clear()  # Clear the previous figure

    fig, ax = plt.subplots()
    ax.bar(te_types, counts, color='orange')
    ax.set_xlabel("TE Classification", fontweight='bold')
    ax.set_ylabel("Count", fontweight='bold')
    ax.set_title("Complementary TE Counts by Type", fontweight='bold')
    plt.xticks(fontweight='bold')
    plt.yticks(fontweight='bold')
    plt.grid(alpha=0.5)

    # Print memory usage after creating the bar chart
    print_memory_usage("After creating bar chart")

# Run the analysis and benchmark memory usage
if __name__ == "__main__":
    print_memory_usage("At start of script")
    analyse_sequence()
    print_memory_usage("At end of script")


# Memory usage at At start of script: 56.29 MB
# Memory usage at After clearing previous results: 56.31 MB
# Memory usage at After TE analysis: 56.31 MB
# Memory usage at After creating bar chart: 71.11 MB
# Memory usage at At end of script: 71.11 MB
