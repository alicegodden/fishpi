import tkinter as tk
from tkinter import messagebox, simpledialog
from PIL import Image, ImageTk
import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# test piRNA sequence should give, sequence: TACACGAAGACTGTGGTGTGATTGGGCG
# result: dna 1223, ltr 369, line 136, sine 42, rc 15, satellite 8 counts

# Create a list to store complementary TE sequences
complementary_TE_list = []

# Initialize piRNA_substring as an empty string
piRNA_substring = ""

# Create global variables for the bar chart and popup window
fig = None
popup_window = None

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
def analyze_sequence():
    global complementary_TE_list, piRNA_substring, fig, popup_window, counted_te_names

    # Clear the previous results and reset te_type_counts
    complementary_TE_list.clear()
    for te_type in te_type_counts:
        te_type_counts[te_type] = 0

    # Reset the set of counted TE names
    counted_te_names = set()

    # Get the new piRNA sequence
    piRNA_sequence = piRNA_entry.get()
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
    ax.set_xlabel("TE Class", fontweight='bold')
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
    popup_window.geometry("800x700")

    canvas = FigureCanvasTkAgg(fig, master=popup_window)
    canvas.get_tk_widget().pack()

    # Display the number of complementary TE sequences in the popup window
    te_count_text = "\n".join([f"{te_type}: {te_count}" for te_type, te_count in te_type_counts.items()])

    # Calculate the total count of all TEs
    total_te_count = sum(te_type_counts.values())
    te_count_text += f"\nTotal TE count: {total_te_count}"  # Add the total count to the label

    counts_label = tk.Label(popup_window, text=te_count_text)
    counts_label.pack()

    # Add an export button for the bar chart
    export_chart_button = tk.Button(popup_window, text="Export Bar Chart", command=export_bar_chart, bg="orange",
                                    fg="navy", borderwidth="0")
    export_chart_button.pack()

    # Add an export button for the TE list
    export_te_list_button = tk.Button(popup_window, text="Export complementary piRNA:TEs", command=export_complementary_TE_list, bg="orange",
                                      fg="navy", borderwidth="0")
    export_te_list_button.pack()

# Define a function to export the complementary TE list to a CSV file
def export_complementary_TE_list():
    global complementary_TE_list
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

# Create the main application window FishPi
app = tk.Tk()
app.title("FishPi: piRNA Sequence Analyzer")
app.geometry("600x600")

load2 = Image.open("fishpi.png")
load2 = load2.resize((400, 400))
render2 = ImageTk.PhotoImage(load2)
img_label2 = tk.Label(image=render2, bg="orange")
img_label2.image = render2
img_label2.place(x=100, y=150)

app.configure(bg="orange")

piRNA_label = tk.Label(app, text="Enter piRNA Sequence:", bg="orange", fg="navy")
piRNA_label.pack()

piRNA_entry = tk.Entry(app, bg="white", fg="navy", borderwidth="0")
piRNA_entry.pack()

analyze_button = tk.Button(app, text="Submit analysis", command=analyze_sequence, bg="orange", fg="navy",
                           borderwidth="0")
analyze_button.pack()

app.mainloop()
