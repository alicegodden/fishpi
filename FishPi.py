# Title: FishPi
# First Author : Dr. Alice M. Godden

import tkinter as tk
from tkinter import messagebox, simpledialog, filedialog
from PIL import Image, ImageTk
import csv
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.lines import Line2D


# test piRNA sequence should give, sequence: TACACGAAGACTGTGGTGTGATTGGGCG
# result: dna 1223, ltr 369, line 136, sine 42, rc 15, satellite 8 counts

# Initialize global variables
piRNA_substring = ""
complementary_TE_list = []
fig = None
fig_chrom = None
popup_window = None
uploaded_pirna_sequences = []

# Define a function to check for piRNA TE complementarity
def is_complementary(piRNA_substring, TE_sequence):
    complementary_sequence = "".join(get_complementary_base(base) for base in piRNA_substring)
    return complementary_sequence in TE_sequence

# Define a function to get the complementary base
def get_complementary_base(base):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement.get(base, base)

# Create a dictionary to count TE types
te_type_counts = {
    "DNA": 0,
    "LTR": 0,
    "LINE": 0,
    "SINE": 0,
    "RC": 0,
    "SATELLITE": 0,
}

# Define a function to classify TE type based on keywords
counted_te_names = set()

def classify_te_type(TE_name):
    TE_name_lower = TE_name.lower()

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
        if any(keyword.lower() in TE_name_lower for keyword in keyword_list):
            if TE_name not in counted_te_names:
                counted_te_names.add(TE_name)
                return te_type
            else:
                return None

    return None

# Define a function to perform the analysis
def analyse_sequence():
    global complementary_TE_list, piRNA_substring, fig, popup_window, counted_te_names, fig_chrom

    # Clear the previous results and reset te_type_counts
    complementary_TE_list.clear()
    for te_type in te_type_counts:
        te_type_counts[te_type] = 0

    # Reset the set of counted TE names
    counted_te_names = set()

    # Get the new piRNA sequence
    piRNA_sequence = piRNA_entry.get()
    if not piRNA_sequence and uploaded_pirna_sequences:
        piRNA_sequence = uploaded_pirna_sequences.pop(0)  # Get the first sequence from the uploaded file
    piRNA_substring = piRNA_sequence[0:9]

    # Open the TE sequence file here
    with open("GRCz11.teseqs.fasta", "r") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                TE_name = line.strip().split()[0]
                TE_sequence = next(fasta_file).strip()

                te_type = classify_te_type(TE_name)
                if te_type in te_type_counts and is_complementary(piRNA_substring, TE_sequence):
                    te_type_counts[te_type] += 1
                    complementary_TE_list.append((TE_name, TE_sequence))

    # Create a bar chart of TE counts for different types with orange bars
    te_types = list(te_type_counts.keys())
    counts = [te_type_counts[te_type] for te_type in te_types]
    colors = ['orange' for _ in te_types]

    if fig:
        fig.clear()  # Clear the previous figure

    fig, ax = plt.subplots()
    ax.bar(te_types, counts, color=colors)
    ax.set_xlabel("TE Classification", fontweight='bold')
    ax.set_ylabel("Count", fontweight='bold')
    ax.set_title("Complementary TE Counts by Type", fontweight='bold')
    plt.xticks(fontweight='bold')
    plt.yticks(fontweight='bold')
    plt.grid(alpha=0.5)

    if popup_window:
        popup_window.destroy()  # Destroy the previous popup window if it exists

    # Create a new popup window for the bar chart
    popup_window = tk.Toplevel(app)
    popup_window.title("Complementary TE Counts")
    popup_window.geometry("1600x800")

    # Layout for the bar chart and buttons
    bar_chart_frame = tk.Frame(popup_window)
    bar_chart_frame.pack(side=tk.LEFT, padx=5, fill=tk.Y)

    # Display the bar chart
    canvas = FigureCanvasTkAgg(fig, master=bar_chart_frame)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    # Create a frame for the export buttons under the bar chart
    export_buttons_frame = tk.Frame(popup_window)
    export_buttons_frame.pack(side=tk.TOP, padx=50, pady=10, fill=tk.Y)

    # Add export buttons for the bar chart
    export_chart_button = tk.Button(export_buttons_frame, text="Export Bar Chart", command=export_bar_chart, bg="orange",
                                    fg="navy", borderwidth="0")
    export_chart_button.pack(pady=5)

    # Add an export button for the complementary TE list
    export_te_list_button = tk.Button(export_buttons_frame, text="Export Complementary piRNA:TEs",
                                      command=export_complementary_TE_list, bg="orange",
                                      fg="navy", borderwidth="0")
    export_te_list_button.pack(pady=5)

    # Plot TE chromosomal locations
    plot_te_chromosomal_locations()
# te colors

def classify_te_type(TE_name):
    TE_name_lower = TE_name.lower()

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
        if any(keyword.lower() in TE_name_lower for keyword in keyword_list):
            return te_type

    return "UNKNOWN"
# Define a function to plot TE chromosomal locations
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

def plot_te_chromosomal_locations():
    global fig_chrom
    # Read chromosome lengths from chr_length.txt
    chromosome_lengths = {}
    with open("chrom_end.txt", "r") as length_file:
        for line in length_file:
            parts = line.strip().split()
            chrom = parts[0]
            length = int(parts[1])
            chromosome_lengths[chrom] = length

    # Read centromere positions from chrcen.txt
    centromere_positions = {}
    with open("chrcen.txt", "r") as cen_file:
        for line in cen_file:
            parts = line.strip().split()
            chrom = parts[0]
            position = int(parts[1])
            centromere_positions[chrom] = position

    # Define TE color mapping with Inferno colormap
    # Create an Inferno colormap
    colormap = plt.get_cmap('PuOr')
    # Extract distinct colors from the colormap
    num_colors = len(te_type_counts) + 1  # +1 for unknown
    colors = [to_hex(colormap(i / num_colors)) for i in range(num_colors)]

    # Create a mapping from TE types to colors
    te_types = list(te_type_counts.keys()) + ['UNKNOWN']
    te_colors = dict(zip(te_types, colors))

    # Initialize the plot
    fig_chrom, ax_chrom = plt.subplots(figsize=(10, 6))

    # Plot chromosomes
    for chrom, length in chromosome_lengths.items():
        ax_chrom.plot([chrom, chrom], [0, length], color='black', lw=2)  # Draw chromosome line

    # Plot centromeres as grey dots
    for chrom, position in centromere_positions.items():
        if chrom in chromosome_lengths:
            length = chromosome_lengths[chrom]
            ax_chrom.plot(chrom, position, 'D', color='grey', markersize=5, label='Centromere', zorder=20)

    # Plot TE locations
    te_types_plotted = set()  # To keep track of TE types for the legend

    for te_name, _ in complementary_TE_list:
        parts = te_name.split("::")
        if len(parts) < 2:
            continue

        chrom = parts[1].split(":")[0]
        positions = parts[1].split(":")[1].split("-")
        start_pos = int(positions[0])
        end_pos = int(positions[1].split("(")[0])

        # Determine TE type
        te_type = classify_te_type(te_name)
        color = te_colors.get(te_type, te_colors['UNKNOWN'])  # Default to 'UNKNOWN' color if type is not in te_colors

        # Plot TE start and end positions
        ax_chrom.plot(chrom, start_pos, 'o', color=color, markersize=5, alpha=0.75, label=f"{te_type}")
        te_types_plotted.add(te_type)

    # Set labels and title
    ax_chrom.set_xticks(list(chromosome_lengths.keys()))
    ax_chrom.set_xticklabels(list(chromosome_lengths.keys()), rotation=0)
    ax_chrom.set_xlabel("Chromosome", fontweight='bold')
    ax_chrom.set_ylabel("Position", fontweight='bold')
    ax_chrom.set_title("Complementary TE genomic locations ", fontweight='bold')
    plt.grid(alpha=0.5)


    # Create custom legend excluding 'UNKNOWN'
    legend_elements = [Line2D([0], [0], marker='D', color='w', markerfacecolor='grey', markersize=8, label='Centromere')]
    for te_type in te_types:
        if te_type != 'UNKNOWN':  # Exclude 'UNKNOWN' from the legend
            color = te_colors.get(te_type, te_colors['UNKNOWN'])
            legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=8, label=te_type))

    ax_chrom.legend(handles=legend_elements)

#####LEAVE
    # Create a frame for the chromosomal plot and export button
    chrom_plot_frame = tk.Frame(popup_window)
    chrom_plot_frame.pack(side=tk.LEFT, padx=5, fill=tk.Y)

    # Display the chromosomal plot
    chrom_canvas = FigureCanvasTkAgg(fig_chrom, master=chrom_plot_frame)
    chrom_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    # Add an export button for the chromosomal plot
    export_chrom_plot_button = tk.Button(chrom_plot_frame, text="Export Chromosomal Plot",
                                         command=lambda: export_chrom_plot(fig_chrom), bg="orange",
                                         fg="navy", borderwidth="0")
    export_chrom_plot_button.pack(pady=5)

# Define a function to export the chromosomal plot as an image
def export_chrom_plot(fig_chrom):
    filename = simpledialog.askstring("Input",
                                      "Enter the file name for the chromosomal plot (e.g., te_chromosomal_plot.png):")
    if filename:
        fig_chrom.savefig(filename, dpi=600, bbox_inches='tight')
        messagebox.showinfo("File Saved", f"Chromosomal plot saved as {filename}")

# Define a function to export the complementary TE list as a CSV file
def export_complementary_TE_list():
    filename = simpledialog.askstring("Input", "Enter the file name (e.g., complementary_TE_list.csv):")
    if filename:
        with open(filename, "w", newline="") as csvfile:
            fieldnames = ["TE Name", "TE Sequence"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for te_name, te_sequence in complementary_TE_list:
                writer.writerow({"TE Name": te_name, "TE Sequence": te_sequence})
        messagebox.showinfo("File Saved", f"Complementary TE List saved to {filename}")

# Define a function to export the bar chart as an image
def export_bar_chart():
    global fig
    filename = simpledialog.askstring("Input", "Enter the file name for the chart (e.g., te_counts_chart.png):")
    if filename:
        # Save the chart with a DPI of 600
        fig.savefig(filename, dpi=600, bbox_inches='tight')
        messagebox.showinfo("File Saved", f"Bar chart saved as {filename}")

# Define a function to upload and process a file containing piRNA sequences
def upload_pirna_file():
    global uploaded_pirna_sequences
    file_path = filedialog.askopenfilename(filetypes=[("FASTA Files", "*.fasta *.fa"), ("All Files", "*.*")])
    if file_path:
        with open(file_path, "r") as file:
            uploaded_pirna_sequences.clear()
            piRNA_sequences = file.read().strip().split('>')
            for sequence in piRNA_sequences:
                if sequence:
                    lines = sequence.splitlines()
                    header = lines[0]
                    piRNA_sequence = ''.join(lines[1:])
                    uploaded_pirna_sequences.append(piRNA_sequence)
        messagebox.showinfo("File Uploaded",
                            "piRNA file uploaded successfully. Click 'Analyse Uploaded piRNA' to proceed.")

# Define a function to analyse the next piRNA sequence from the uploaded file
def analyse_uploaded_pirna():
    if uploaded_pirna_sequences:
        analyse_sequence()
    else:
        messagebox.showwarning("No piRNA Sequences", "No piRNA sequences to analyse. Please upload a file first.")

# Create the main application window FishPi
app = tk.Tk()
app.title("FishPi: piRNA:TE Complementarity Analyser")
app.geometry("600x600")

load2 = Image.open("fishpi.png")
load2 = load2.resize((400, 400))
render2 = ImageTk.PhotoImage(load2)
img_label2 = tk.Label(image=render2, bg="orange")
img_label2.image = render2
img_label2.place(x=100, y=200)

app.configure(bg="orange")

piRNA_label = tk.Label(app, text="Enter piRNA Sequence:", bg="orange", fg="navy", font=("bold"))
piRNA_label.pack()
piRNA_label.place(x=200, y=30)

piRNA_entry = tk.Entry(app, bg="white", fg="navy", borderwidth="0")
piRNA_entry.pack()
piRNA_entry.place(x=200, y=60, width=200)

analyse_button = tk.Button(app, text="Analyse piRNA sequence", command=analyse_sequence, bg="orange", fg="navy",
                           borderwidth="0")
analyse_button.pack()
analyse_button.place(x=200, y=90, width=200)

# Add an upload button for piRNA file
upload_button = tk.Button(app, text="Upload piRNA fasta file", command=upload_pirna_file, bg="orange", fg="navy",
                          borderwidth="0")
upload_button.pack()
upload_button.place(x=200, y=140, width=200)

# Add a button to analyse the uploaded piRNA file
analyse_uploaded_button = tk.Button(app, text="Analyse Uploaded piRNAs", command=analyse_uploaded_pirna, bg="orange",
                                    fg="navy", borderwidth="0")
analyse_uploaded_button.pack()
analyse_uploaded_button.place(x=200, y=170, width=200)

app.mainloop()
