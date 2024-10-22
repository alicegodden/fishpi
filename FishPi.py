# Title: FishPi
# Authors: Dr. Alice M. Godden & Dr. Benjamin T. Rix
# Test piRNA sequence should give, sequence: TACACGAAGACTGTGGTGTGATTGGGCG
# Result: dna 1223, ltr 369, line 136, sine 42, rc 15, satellite 8 counts

from tkinter import messagebox, simpledialog, filedialog, ttk
from PIL import Image, ImageTk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.lines import Line2D
from matplotlib.colors import to_hex
from scipy.stats import binom
import csv
from typing import List, Dict, Tuple
import random
import os


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
pie_slices = []  # List to keep track of pie slice labels

try:
    import Tkinter as tk  # Python 2
except ImportError:
    import tkinter as tk  # Python 3

    # same as 2> /dev/null in the shell
    f = open("/dev/null", "w")
    os.dup2(f.fileno(), 2)
    f.close()  # to avoid error message about mac compatibility
    # Secure coding is not enabled for restorable state!

# Create a dictionary to count TE types
te_type_counts: Dict[str, int] = {
    "DNA": 0,
    "LTR": 0,
    "LINE": 0,
    "SINE": 0,
    "RC": 0,
    "SATELLITE": 0,
}

# Define file paths for each species
species_files: Dict[str, Dict[str, str]] = {
    "zebrafish": {
        "te_fasta": "GRCz11_ensembl_teseqs.fishpi.fasta",
        "chrom_end": "zebrafish_chrom_end.txt",
        "chrcen": "zebrafish_chrcen.txt",
    },
    "tilapia": {
        "te_fasta": "Onil_1.2_ensembl_teseqs.fishpi.fasta",
        "chrom_end": "tilapia_chrom_end.txt",
    },
    "medaka": {
        "te_fasta": "oryLat2_ensembl_teseqs.fishpi.fasta",
        "chrom_end": "medaka_chrom_end.txt",
    },
}


# Function to classify TE type
def classify_te_type(te_name: str) -> str:
    te_name_lower = te_name.lower()
    keyword_mapping = [
        (["hat", "tc1", "tcmar", "harbinger", "enspm", "kolobok", "merlin",
          "crypton", "piggybac", "dada", "zatar", "ginger", "tdr", "polinton",
          "maverick", "acrobat", "looper", "tzf", "angel", "mariner", "dna"],
         "DNA"),
        (["gypsy", "dirs", "ngaro", "erv", "pao", "copia", "bel", "herv",
          "bhikari", "ltr"], "LTR"),
        (["l2", "l1", "l1-tx1", "rex-babar", "rte", "penelope", "keno", "rex",
          "line"], "LINE"),
        (["alu", "trna-v-rte", "sine"], "SINE"),
        (["helitron", "rc"], "RC"),
        (["brsati", "mosat", "sat"], "SATELLITE"),
    ]

    for keyword_list, te_type in keyword_mapping:
        if any(keyword.lower() in te_name_lower for keyword in keyword_list):
            return te_type
    return "UNKNOWN"


# Define a function to introduce mismatches in the seed sequence
def introduce_mismatches(seed: str, num_mismatches: int) -> str:
    """Introduce mismatches into the piRNA seed region."""
    if num_mismatches == 0:
        return seed  # No mismatches, return the seed unchanged
    seed_list = list(seed)
    indices_to_mutate = random.sample(range(len(seed)), num_mismatches)
    bases = ['A', 'T', 'C', 'G']
    for i in indices_to_mutate:
        current_base = seed_list[i]
        possible_bases = [b for b in bases if b != current_base]
        seed_list[i] = random.choice(possible_bases)
    return ''.join(seed_list)


# Define a function to check for piRNA TE complementarity
def is_complementary(piRNA_substring: str, te_sequence: str,
                     reverse_complement: bool = False) -> bool:
    complementary_sequence = "".join(
        get_complementary_base(base) for base in piRNA_substring)
    if reverse_complement:
        complementary_sequence = complementary_sequence[::-1]
    return complementary_sequence in te_sequence


# Define a function to get the complementary base
def get_complementary_base(base: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement.get(base, base)


# Function to perform the analysis
def analyse_sequence() -> None:
    global complementary_TE_list, piRNA_substring, fig, counted_te_names
    global total_lines

    # Clear previous results and reset TE type counts
    complementary_TE_list.clear()
    for te_type in te_type_counts:
        te_type_counts[te_type] = 0

    counted_te_names.clear()

    # Get the new piRNA sequence
    piRNA_sequence = piRNA_entry.get()
    if not piRNA_sequence and uploaded_pirna_sequences:
        piRNA_sequence = uploaded_pirna_sequences.pop(0)

    # Get the range for the piRNA substring from the user input
    range_input = search_range_entry.get()
    start, end = 0, 9  # Default range to match `0:9`
    try:
        if range_input:
            start, end = map(int, range_input.split('-'))
            start -= 1
            end -= 1  # Adjust to make inclusive in a one-based system
        piRNA_substring = piRNA_sequence[start:end]

    except ValueError:
        messagebox.showerror("Invalid", "Please enter"
                                        " a valid range like '1-10'.")
        return

    # Get the number of mismatches from user input
    mismatches = 0
    mismatch_input = mismatch_entry.get()
    try:
        if mismatch_input:
            mismatches = int(mismatch_input)
            if mismatches < 0 or mismatches > 15:
                messagebox.showerror("Invalid Input",
                                     "Please enter a number between 0 and 15.")
                return
    except ValueError:
        messagebox.showerror("Invalid Input", "Please enter a valid number.")
        return

    # Introduce mismatches if applicable
    piRNA_substring_with_mismatches = introduce_mismatches(
        piRNA_substring, mismatches)

    # Determine if reverse complement should be searched
    reverse_complement_search = search_reverse_complement.get()

    # Determine the selected species for TE reference file
    species = species_var.get()
    te_file = species_files[species]["te_fasta"]

    # Create progress window
    create_progress_window()

    # Open the TE sequence file and analyze it
    try:
        with open(te_file, "r") as fasta_file:
            total_lines = sum(1 for line in fasta_file)  # Count total lines
            fasta_file.seek(0)  # Reset file pointer

            for line_number, line in enumerate(fasta_file):
                if line.startswith(">"):
                    te_name = line.strip().split()[0]
                    te_sequence = next(fasta_file).strip()

                    te_type = classify_te_type(te_name)
                    if te_type in te_type_counts and (
                            is_complementary(piRNA_substring_with_mismatches,
                                             te_sequence, False) or
                            (reverse_complement_search and is_complementary
                                (piRNA_substring_with_mismatches,
                                 te_sequence, True))):
                        te_type_counts[te_type] += 1
                        complementary_TE_list.append((te_name, te_sequence))

                # Update progress if necessary
                if line_number % 100 == 0 or line_number == total_lines - 1:
                    progress_value = (line_number + 1) / total_lines * 100
                    update_progress_bar(progress_value)
                    move_fish(progress_value)

        # Finalize progress bar update and create the result
        update_progress_bar(100)  # Ensure progress bar reaches 100%
        move_fish(100)
        progress_window.after(500, progress_window.destroy)
        create_results()

    except FileNotFoundError:
        messagebox.showerror("File Not Found",
                             f"The specified TE sequence file for {species} "
                             f"was not found.")
        progress_window.destroy()
        return
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {str(e)}")
        progress_window.destroy()
        return


# Define a function to create the progress window
def create_progress_window() -> None:
    global progress_window, fish_label, progress_bar, fish_open_photo
    global fish_closed_photo, pie_slices
    progress_window = tk.Toplevel(app)
    progress_window.title("Analyzing...")
    progress_window.geometry("300x150")

    # Progress bar style
    style = ttk.Style()
    style.theme_use('clam')
    style.configure(
        "Custom.Horizontal.TProgressbar",
        troughcolor="white",
        background="orange",
        thickness=20,
    )

    progress_bar = ttk.Progressbar(progress_window,
                                   style="Custom.Horizontal.TProgressbar",
                                   length=250, mode='determinate')
    progress_bar.pack(pady=20)

    # Create pie slices at every 10% mark, except for 100%
    pie_slices.clear()  # Clear previous slices
    for i in range(1, 10):  # Create pie slices
        pie_slice_image = Image.open("pieslice.png")
        pie_slice_image = pie_slice_image.resize((25, 25), Image.LANCZOS)
        pie_slice_photo = ImageTk.PhotoImage(pie_slice_image)
        pie_slice_label = tk.Label(progress_window, image=pie_slice_photo)
        pie_slice_label.image = pie_slice_photo
        pie_slice_label.place(x=(i * 25), y=70)  # Adjust the x position
        pie_slices.append(pie_slice_label)

    # Load fish images
    fish_open_image = (
        Image.open("fish_open.png").resize((50, 30), Image.LANCZOS))
    global fish_open_photo
    fish_open_photo = ImageTk.PhotoImage(fish_open_image)

    fish_closed_image = Image.open("fish_closed.png")
    fish_closed_image = fish_closed_image.resize((50, 30), Image.LANCZOS)
    global fish_closed_photo
    fish_closed_photo = ImageTk.PhotoImage(fish_closed_image)

    fish_label = tk.Label(progress_window, image=fish_open_photo)
    fish_label.place(x=0, y=70)

    # Add transparent pie at the end
    pie_transparent_image = Image.open("pie_transparent.png")
    pie_transparent_image = pie_transparent_image.resize((35, 35),
                                                         Image.LANCZOS)
    pie_transparent_photo = ImageTk.PhotoImage(pie_transparent_image)
    pie_transparent_label = tk.Label(progress_window,
                                     image=pie_transparent_photo)
    pie_transparent_label.image = pie_transparent_photo
    pie_transparent_label.place(x=260, y=70)


# Define a function to update the progress bar
def update_progress_bar(value: float) -> None:
    progress_bar['value'] = value
    progress_window.update_idletasks()  # Refresh the GUI
    move_fish(value)  # Update fish movement


# Move the fish according to progress
def move_fish(progress_value: float) -> None:
    max_width = 250
    fish_x = int((progress_value / 100) * max_width)
    fish_label.place(x=fish_x, y=70)

    # Alternate fish images based on progress
    if int(progress_value) % 2 == 0:
        fish_label.config(image=fish_open_photo)
    else:
        fish_label.config(image=fish_closed_photo)

    # Remove pie slices as the fish reaches them
    for index in range(len(pie_slices)):
        pie_slice_x = index * 25
        if fish_x >= pie_slice_x:
            if pie_slices[index]:
                pie_slices[index].destroy()
                pie_slices[index] = None


# Create results window with plots
def create_results() -> None:
    global fig, fig_chrom, popup_window

    # Create bar chart of TE counts using PuOr colormap
    te_types = list(te_type_counts.keys())
    counts = [te_type_counts[te_type] for te_type in te_types]

    # Use PuOr colormap for bar chart colors
    colormap = plt.get_cmap('Spectral')
    num_colors = len(te_type_counts)
    colors = [to_hex(colormap(i / num_colors)) for i in range(num_colors)]

    if fig:
        fig.clear()

    fig, ax = plt.subplots()
    ax.bar(te_types, counts, color=colors)
    ax.set_xlabel("TE Classification", fontweight='bold')
    ax.set_ylabel("Count", fontweight='bold')
    ax.set_title("Complementary TE Counts by Type", fontweight='bold')
    plt.xticks(fontweight='bold')
    plt.yticks(fontweight='bold')
    plt.grid(alpha=0.5)

    # Add enrichment information to bar chart using a binomial test
    total_te_types = sum(te_type_counts.values())
    for i, te_type in enumerate(te_types):
        te_count = te_type_counts[te_type]
        p_value = binom.sf(te_count - 1, total_te_types, 1 / len(te_types))
        if p_value < 0.05:
            ax.text(i, counts[i] + 0.1, '*', color='black', fontsize=20,
                    ha='center')

    if popup_window:
        popup_window.destroy()

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

    export_chart_button = tk.Button(export_buttons_frame,
                                    text="Export Bar Chart",
                                    command=export_bar_chart,
                                    bg="orange", fg="navy", borderwidth="0")
    export_chart_button.pack(pady=5)

    export_te_list_button = tk.Button(export_buttons_frame,
                                      text="Export Complementary piRNA:TEs",
                                      command=export_complementary_te_list,
                                      bg="orange", fg="navy",
                                      borderwidth="0")
    export_te_list_button.pack(pady=5)

    export_chrom_plot_button = tk.Button(popup_window,
                                         text="Export Chromosomal Plot",
                                         command=lambda: export_chrom_plot(
                                             fig_chrom),
                                         bg="orange", fg="navy",
                                         borderwidth="0")
    export_chrom_plot_button.pack()

    plot_te_chromosomal_locations()


# Define a function to plot TE chromosomal locations with PuOr colormap
def plot_te_chromosomal_locations() -> None:
    global fig_chrom

    # Determine the selected species for TE reference file
    species = species_var.get()
    chrom_end_file = species_files[species]["chrom_end"]
    use_centromere = "chrcen" in species_files[species]
    chrcen_file = species_files[species].get("chrcen")

    # Read chromosome lengths from the selected chromosome end file
    chromosome_lengths = {}
    try:
        with open(chrom_end_file, "r") as length_file:
            for line in length_file:
                parts = line.strip().split()
                chrom = parts[0]
                length = int(parts[1])
                chromosome_lengths[chrom] = length
    except FileNotFoundError:
        messagebox.showerror("File Not Found",
                             f"The chromosome end file for {species} "
                             f"was not found.")
        return

    centromere_positions = {}
    if use_centromere:
        try:
            with open(chrcen_file, "r") as cen_file:
                for line in cen_file:
                    parts = line.strip().split()
                    chrom = parts[0]
                    position = int(parts[1])
                    centromere_positions[chrom] = position
        except FileNotFoundError:
            messagebox.showwarning("File Not Found",
                                   "Centromere data file not found, "
                                   "centromere info will be excluded.")

    if fig_chrom:
        fig_chrom.clear()

    fig_chrom, ax_chrom = plt.subplots(figsize=(10, 8))

    # Plot chromosomes
    x_pos = list(range(len(chromosome_lengths)))
    for idx, (chrom, length) in enumerate(chromosome_lengths.items()):
        ax_chrom.plot([x_pos[idx], x_pos[idx]], [0, length],
                      color='black', lw=2)

        # Plot centromeres if available
        if use_centromere and chrom in centromere_positions:
            centromere_position = centromere_positions[chrom]
            ax_chrom.plot(x_pos[idx], centromere_position, 'D',
                          color='grey', markersize=8, zorder=20)

    # Use PuOr colormap for TE locations
    colormap = plt.get_cmap('Spectral')
    num_colors = len(te_type_counts) + 1
    colors = [to_hex(colormap(i / num_colors)) for i in range(num_colors)]

    # Create a mapping from TE types to colors
    te_types = list(te_type_counts.keys())
    te_colors = dict(zip(te_types, colors))

    # Plot TE locations with colors from the te_colors dictionary
    te_types_plotted = set()
    for te_name, _ in complementary_TE_list:
        parts = te_name.split("::")
        if len(parts) < 2:
            continue

        chrom = parts[1].split(":")[0]
        try:
            positions = parts[1].split(":")[1].split("-")
            start_pos = int(positions[0])
            if chrom in chromosome_lengths:
                chrom_index = list(chromosome_lengths.keys()).index(chrom)
                te_type = classify_te_type(te_name)
                color = te_colors.get(te_type, 'black')
                ax_chrom.plot(chrom_index, start_pos, 'o', color=color,
                              markersize=5, alpha=0.75)
                te_types_plotted.add(te_type)
        except (IndexError, ValueError):
            continue

    # Set labels and title
    ax_chrom.set_xticks(x_pos)
    ax_chrom.set_xticklabels(list(chromosome_lengths.keys()),
                             rotation=45, ha='right')
    ax_chrom.set_xlabel("Chromosome", fontweight='bold')
    ax_chrom.set_ylabel("Position", fontweight='bold')
    ax_chrom.set_title("Complementary TE Genomic Locations", fontweight='bold')
    plt.grid(alpha=0.5)

    # Custom legend ensuring that colors match TE types exactly
    legend_elements = []
    if use_centromere:
        legend_elements.append(Line2D([0], [0], marker='D', color='w',
                                      markerfacecolor='grey', markersize=8,
                                      label='Centromere'))
    for te_type in te_types:
        if te_type in te_types_plotted:
            color = te_colors[te_type]
            legend_elements.append(Line2D([0], [0], marker='o', color='w',
                                          markerfacecolor=color,
                                          markersize=8, label=te_type))

    ax_chrom.legend(handles=legend_elements, loc='upper right')

    # Add to the GUI
    chrom_frame = tk.Frame(popup_window)
    chrom_frame.pack(side=tk.RIGHT, padx=5, fill=tk.Y)

    chrom_canvas = FigureCanvasTkAgg(fig_chrom, master=chrom_frame)
    chrom_canvas.draw()
    chrom_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    export_chrom_plot_button = tk.Button(chrom_frame,
                                         text="Export Chromosomal Plot",
                                         command=lambda: export_chrom_plot(
                                             fig_chrom),
                                         bg="orange", fg="navy",
                                         borderwidth="0")
    export_chrom_plot_button.pack(pady=5)


# Define a function to export the chromosomal plot as an image
def export_chrom_plot(fig_chrom):
    filename = simpledialog.askstring(
        "Input", "Enter the file name for the chromosomal plot "
        "(e.g., te_chromosomal_plot.png):"
    )
    if filename:
        fig_chrom.savefig(filename, dpi=600, bbox_inches='tight')
        messagebox.showinfo("File Saved",
                            f"Chromosomal plot saved as {filename}")


# Define a function to export the complementary TE list as a CSV file
def export_complementary_te_list() -> None:
    filename = simpledialog.askstring("Input",
                                      "Enter the file name "
                                      "(e.g., complementary_TE_list.csv):")
    if filename:
        with open(filename, "w", newline="") as csvfile:
            fieldnames = ["TE Name", "TE Sequence"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for te_name, te_sequence in complementary_TE_list:
                writer.writerow({"TE Name": te_name,
                                 "TE Sequence": te_sequence})
        messagebox.showinfo("File Saved",
                            f"Complementary TE List saved to {filename}")


# Define a function to export the bar chart as an image
def export_bar_chart() -> None:
    global fig
    filename = simpledialog.askstring("Input",
                                      "Enter the file name for the chart "
                                      "(e.g., te_counts_chart.png):")
    if filename:
        fig.savefig(filename, dpi=600, bbox_inches='tight')
        messagebox.showinfo("File Saved", f"Bar chart saved as {filename}")


# Define a function to upload and process a file containing piRNA sequences
def upload_pirna_file() -> None:
    global uploaded_pirna_sequences
    file_path = filedialog.askopenfilename(filetypes=[("FASTA Files",
                                                       "*.fasta *.fa"),
                                                      ("All Files", "*.*")])
    if file_path:
        with open(file_path, "r") as file:
            uploaded_pirna_sequences.clear()
            piRNA_sequences = file.read().strip().split('>')
            for sequence in piRNA_sequences:
                if sequence:
                    lines = sequence.splitlines()
                    piRNA_sequence = ''.join(lines[1:])
                    uploaded_pirna_sequences.append(piRNA_sequence)
        messagebox.showinfo("File Uploaded",
                            "piRNA file uploaded successfully. "
                            "Click 'Analyse Uploaded piRNA' to proceed.")


# Define a function to analyse the next piRNA sequence from the uploaded file
def analyse_uploaded_pirna() -> None:
    if uploaded_pirna_sequences:
        analyse_sequence()
    else:
        messagebox.showwarning("No piRNA Sequences",
                               "No piRNA sequences to analyse. "
                               "Please upload a file first.")


# Create the main application window FishPi
app = tk.Tk()
app.title("FishPi: piRNA:TE Complementarity Analyser")
app.geometry("600x700")

load2 = Image.open("fishpi.png")
load2 = load2.resize((250, 250))
render2 = ImageTk.PhotoImage(load2)
img_label2 = tk.Label(image=render2, bg="orange")
img_label2.image = render2
img_label2.place(x=175, y=475)

app.configure(bg="orange")

piRNA_label = tk.Label(app, text="Enter piRNA Sequence:",
                       bg="orange", fg="navy", font="bold")
piRNA_label.pack()
piRNA_label.place(x=200, y=30)

piRNA_entry = tk.Entry(app, bg="white", fg="navy", borderwidth="0")
piRNA_entry.pack()
piRNA_entry.place(x=200, y=60, width=200)

# Add text box for range input
range_label = tk.Label(app, text="Enter Base Range (e.g., 1-10):",
                       bg="orange", fg="navy", font="bold")
range_label.pack()
range_label.place(x=150, y=100)

search_range_entry = tk.Entry(app, bg="white", fg="navy", borderwidth="0")
search_range_entry.pack()
search_range_entry.place(x=200, y=130, width=100)

# Add checkbox for reverse complement search
search_reverse_complement = tk.BooleanVar()
reverse_complement_check = tk.Checkbutton(app,
                                          text="Also Search for Reverse "
                                          "Complement of Seed",
                                          variable=search_reverse_complement,
                                          bg="orange", fg="navy", font="bold")
reverse_complement_check.pack()
reverse_complement_check.place(x=150, y=160)

# Add input box for mismatches
mismatch_label = tk.Label(app,
                          text="Enter Number of Mismatches to seed "
                          "(Default is 0):", bg="orange", fg="navy",
                          font="bold")
mismatch_label.pack()
mismatch_label.place(x=150, y=200)

mismatch_entry = tk.Entry(app, bg="white", fg="navy", borderwidth="0")
mismatch_entry.pack()
mismatch_entry.place(x=200, y=230, width=100)

# Add radio buttons for species selection
species_var = tk.StringVar(value="zebrafish")
species_label = tk.Label(app, text="Select Species:",
                         bg="orange", fg="navy", font="bold")
species_label.pack()
species_label.place(x=150, y=270)

species_frame = tk.Frame(app, bg="orange")
species_frame.pack()
species_frame.place(x=150, y=300)

for species in species_files.keys():
    rb = tk.Radiobutton(species_frame, text=species.capitalize(),
                        variable=species_var, value=species,
                        bg="orange", fg="navy", font="bold")
    rb.pack(anchor="w")

analyse_button = tk.Button(app, text="Analyse piRNA sequence",
                           command=analyse_sequence,
                           bg="white", fg="navy", borderwidth="0")
analyse_button.pack()
analyse_button.place(x=200, y=380, width=200)

# Add an upload button for piRNA file
upload_button = tk.Button(app, text="Upload piRNA fasta file",
                          command=upload_pirna_file, bg="white",
                          fg="navy", borderwidth="0")
upload_button.pack()
upload_button.place(x=200, y=410, width=200)

# Add a button to analyse the uploaded piRNA file
analyse_uploaded_button = tk.Button(app,
                                    text="Analyse Uploaded piRNAs",
                                    command=analyse_uploaded_pirna,
                                    bg="white", fg="navy", borderwidth="0")
analyse_uploaded_button.pack()
analyse_uploaded_button.place(x=200, y=440, width=200)

app.mainloop()
