# FishPi
piRNA TE analysis tool

Project Title: FishPi - piRNA Sequence Analyzer

# Description:
FishPi is a Python-based piRNA sequence analysis tool designed to identify complementary transposable element (TE) sequences in a given piRNA sequence. This project combines the power of various Python libraries, including tkinter for the graphic user interface, matplotlib for data visualization, and PIL (Pillow) for image rendering. All results and bar charts can be exported for downstream analyses and publication in hi-resolution.

# Features:

Analyze piRNA sequences for complementary TE matches.
Classify TE types based on keywords.
Generate a bar chart showing complementary TE counts by type.
Export the list of complementary TE sequences to a CSV file.
Export the bar chart as an image file.

# How to Use FishPi:

Input piRNA Sequence: Launch the FishPi application. Enter your piRNA sequence in the provided text field and click the "Submit analysis" button.

View Analysis Results: FishPi will analyze the piRNA sequence and display the results in a popup window. The results include a bar chart showing complementary TE counts by type and the total TE count.

Export Data: You can export the complementary TE sequences to a CSV file or save the bar chart as an image file using the export buttons provided in the popup window.

# Requirements:

Python 3.x
Python libraries: tkinter, matplotlib, Pillow (PIL)
TE sequence file (e.g., "GRCz11.teseqs.fasta") for TE matching

# Getting Started:

Clone or download the FishPi repository to your local machine.

# Install and use Fishpi (Linux)

```
$ git clone https://github.com/alicegodden/fishpi/tree/fishpi/
$ cd fishpi
$ python fishpi.py # opens the GUI
```

Ensure you have Python and the required libraries installed.

Prepare the TE sequence file and place it in the same directory as FishPi.

Run FishPi by executing the Python script (e.g., python fishpi.py).

# License:
This project is distributed under the GNU General Public License (GPL), which ensures that it remains open source and freely accessible for use and modification.

# Author & Contact:
Dr Alice M. Godden

Contact:

Email: alice.godden@uea.ac.uk
Feel free to contribute to this project, report issues, or suggest improvements. FishPi is designed to make piRNA sequence analysis accessible and informative for researchers and bioinformaticians.

Note:
The "fishpi.png" image file is used for the graphical interface, but you can replace it with your project's relevant image.
