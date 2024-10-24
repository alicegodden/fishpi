![GitHub Logo](https://github.com/alicegodden/fishpi/blob/main/files/fishpi.png)
# FishPi
![GitHub Actions Status](https://github.com/alicegodden/fishpi/actions/workflows/python-test.yml/badge.svg)

Created by Dr. Alice M. Godden & Dr. Benjamin Rix, 2024

# FishPi: A Powerful Tool for piRNA Analysis
FishPi is a cutting-edge, Python-powered solution for piRNA sequence analysis, designed to swiftly identify complementary transposable element (TE) sequences from individual piRNA sequences or entire fasta files. Featuring an intuitive graphical user interface built with tkinter, vivid data visualizations through matplotlib, and high-quality image rendering with Pillow, FishPi makes complex genomic analysis accessible and visually compelling. Export your results and high-resolution plots effortlessly, making FishPi ideal for downstream analysis and publication-ready outputs.

# Performance Benchmark and System Requirements
FishPi was benchmarked using GitHub Actions for continuous integration, ensuring reliable performance across different environments. The entire database of piRNA sequences from piRNADB (5000+ piRNA sequences) was tested, requiring less than 3 GB of memory and taking 1-2 minutes to generate plots and outputs on a standard machine. FishPi was developed using Python v3.11. For full requirements, see requirements.txt.

# Features
Comprehensive piRNA Analysis: 
Analyze piRNA sequences to identify complementary transposable element (TE) matches from TE reference sequences.
Customizable Seed Region Parameters: 
Define the length of the piRNA seed region, introduce mismatches for increased flexibility, and optionally search for reverse complement matches. For teleost species we recommend the defaul of 1-10 complementary base pairs, which ensure complementarity between the first 10 nucleotides of the TE. 
TE Classification:
- Automatically classify TE types based on comprehensive keyword matching to assign each TE to categories like DNA, LTR, LINE, SINE, RC, and Satellite.
Data Visualization:
- Generate bar charts to visualize complementary TE counts by type, providing immediate insights into sequence relationships.
Export Functionality:
- Export the list of complementary TE sequences, complete with details, to a CSV file for further analysis.
- Export generated bar charts as high-resolution image files for reports, presentations, and publication purposes.

# Input file preparation
Users may have piRNA sequences for their given species, if not we recommend the following databases:

SmallRNAGroup's databse of species piRNAs and piRNA clusters is linked [here](https://www.smallrnagroup.uni-mainz.de/piRNAclusterDB/) 

See Drerio_piRNAs_test_sRNAlab.fasta as example for required fasta format input for Fishpi piRNA sequences.

piRBase can be accessed [here](http://bigdata.ibp.ac.cn/piRBase/) 

If you would like to use FishPi on a reference genome other than Zebrafish (GRCz11), Medaka (oryLat2) or Tilapia (Onil_1.2), please prepare your "teseqs.fasta" file by following the same instructions written below, but by selecting appropriate organism on UCSC browser. Also check the piRNA:TE seed rules of your organism, as these can be customised pm the FishPi GUI. If you are working with the Zebrafish, Medaka or Tilapia please download all the annotated teseqs.fasta files [here](https://zenodo.org/record/13911872).

To prepare your teseqs.fasta file:
1. To generate a fasta file of the TE sequences for your reference genome, first make a .bed file. The UCSC table genome browser was used to generate the TE.bed file with these options: clade: Vertebrate, group: Variation and Repeats, genome: Zebrafish, assembly: May 2017 GRCz11, track: RepeatMasker, table: rmsk.
2. To extract DNA sequences from the reference genome based on the co-ordinates supplied in the bed file the following command was used:
   bedtools getfasta -s -name -fi Danio_rerio.GRCz11.dna.primary_assembly.fa -fo GRCz11.teseqs.use.fasta -bed GRCz11.teannotation.bed


Example output files are based on dre-piRNA-1 5'-TGGTTAGTACTTGGATGGGAGACCGCCTGGG-3', taken from piRBase. 


# Install and use FishPi (Linux)

```
# First, clone the github repo for fishpi
git clone https://github.com/alicegodden/fishpi

cd fishpi # navigate to FishPi directory
cd files # navigate to FishPi scripts and files

conda env create -f environment.yml # Create your conda environment
conda activate fishpi_environment # activate your environment

# Verify the Installation
# After activating the environment, verify that everything is installed correctly:
python --version
conda list

# To download the TE sequence fasta files, download all files
wget https://zenodo.org/records/13911872/files/GRCz11_ensembl_teseqs.fishpi.fasta.gz # For Zebrafish
wget https://zenodo.org/records/13911872/files/Onil_1.2_ensembl_teseqs.fishpi.fasta.gz # For Tilapia
wget https://zenodo.org/records/13911872/files/oryLat2_ensembl_teseqs.fishpi.fasta.gz # For Medaka

# Then uncompress fasta files
gunzip *.fasta.gz

# You should now be ready to run FishPi
python FishPi.py # opens the GUI to use FishPi


```
# Using FishPi with PyCharm IDE
To get started with FishPi in PyCharm:

Download Project Files: Make sure you download all necessary files from  [here](https://github.com/alicegodden/fishpi/tree/main/files). Place these files into your working directory within PyCharm, ensuring that they are accessible to your Python scripts.

Download TE Sequence Files: Additionally, download the required TE sequence files from Zenodo  [here](https://zenodo.org/record/13911872). Extract these files to uncompress them, and place them in the same working directory for seamless integration with FishPi. 

# How to Use FishPi

Input piRNA Sequence:
Launch the FishPi application, and enter your piRNA sequence in the provided text box. Click on the "Analyse piRNA sequence" button to start the analysis.
(Tip: You can use your own sequences or obtain them from public repositories like piRbase.)

View Analysis Results:
Once the analysis is complete, FishPi will display the results in a new popup window. You will see a bar chart representing complementary TE counts by type, as well as a chromosomal location plot for complementary TEs.

Export Data:
Export options are provided directly in the popup window for easy sharing or further analysis. 
You can save:
The list of complementary TE sequences as a CSV file.
The generated plots as an image file (high resolution, suitable for publications).
(Figure 1 shows an example of the results display and export options. Figure 2 presents a flow chart summarizing the key steps and processes undertaken by FishPi.)


*Figure 1- Graphic User interface for FishPi and results for dre-piRNA 5'-TACACGAAGACTGTGGTGTGATTGGGCG-3'*
![FishPi_GUI](https://github.com/user-attachments/assets/a7e64def-29a4-4aae-b883-5335beca54b8)


# Flow chart of FishPi functions
Below is a flow chart illustrating the key steps and processes involved in the FishPi workflow. This flow chart provides an overview of the sequence analysis, transposable element (TE) matching, and the export functionalities available in FishPi.

![fishpi_flow_diagram](https://github.com/user-attachments/assets/c998d404-6c3f-4ca2-a335-fa11581721df)
*Figure 2- Flow chart covering key steps and processes by FishPi*

# Output results
The transposable element (TE) species are grouped by order in the bar chart displayed, as shown in Figure 1 above. The TEs analyzed by FishPi are categorized based on RepeatMasker classifications, grouped as follows:

DNA Transposons: hAT, Tc1, Tc-Mar, Harbinger, Enspm, Kolobok, Merlin, Crypton, PiggyBac, Dada, Zatar, Ginger, TDR, Polinton, Maverick, Acrobat, Looper, TZF, Angel, Mariner
LTR Retrotransposons: Gypsy, DIRS, Ngaro, ERV, Pao, Copia, BEL, HERV, Bhikari
LINEs (Long Interspersed Nuclear Elements): L1, L2, L1-Tx1, Rex-Babar, RTE, Penelope, Keno, Rex
SINEs (Short Interspersed Nuclear Elements): Alu, tRNA-V-RTE
Rolling Circle (RC) Transposons: Helitron
Satellite Repeats: BRSATI, MOSAT
These groupings help to identify and quantify complementary TEs in the analyzed piRNA sequence, providing a comprehensive overview of TE relationships.

# Licence
This project is distributed under the GNU General Public License (GPL), which ensures that it remains open source and freely accessible for use.

# Citation
If you use FishPi in your research, please cite:

FishPi: A Bioinformatic Prediction Tool to Link piRNA and Transposable Elements in Zebrafish
Alice May Godden, Benjamin Rix, Simone Immler
bioRxiv 2024.09.10.612046; DOI: https://doi.org/10.1101/2024.09.10.612046


# Author & Contact
Dr. Alice M. Godden
Dr. Benjamin Rix
Prof. Simone Immler
University of East Anglia, School of Biological Sciences, Norwich, United Kingdom, NR4 7TJ

For inquiries or to report issues, feel free to reach out:

Email: alice.godden@uea.ac.uk
We encourage feedback, suggestions for improvements, and collaboration to help make FishPi a more effective tool for piRNA sequence analysis and bioinformatics research.

