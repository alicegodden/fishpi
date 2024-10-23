![GitHub Logo](https://github.com/alicegodden/fishpi/blob/main/files/fishpi.png)
# FishPi

Created by Dr. Alice M. Godden & Dr. Benjamin Rix, 2024

# Description:
FishPi is a Python-based piRNA sequence analysis tool designed to identify complementary transposable element (TE) sequences in a given piRNA sequence or fasta file of piRNAs. This project combines the power of various Python libraries, including tkinter for the graphic user interface, matplotlib for data visualization, and PIL (Pillow) for image rendering. All results and plots can be exported for downstream analyses and publication in hi-resolution.

# Memory and system requirements
FishPi was benchmarked with the entire database of piRNA sequences from piRNADB (5000+ piRNA sequences), this required < 3 GB to run, taking 1-2 minutes to generate plots and outputs for export on a standard machine. FishPi was developed on Python v 3.11.


# Features:

Analyze piRNA sequences for complementary TE matches.
Classify TE types based on keywords.
Generate a bar chart showing complementary TE counts by type.
Export the list of complementary TE sequences to a CSV file.
Export the bar chart as an image file.

# Requirements and input file preparation:
User may have piRNA sequences for their given species, if not we recommend the following databases:

SmallRNAGroup's databse of species piRNAs and piRNA clusters: https://www.smallrnagroup.uni-mainz.de/piRNAclusterDB/ 

See Drerio_piRNAs_test_sRNAlab.fasta as example for required fasta format input for Fishpi piRNA sequences.

piRBase http://bigdata.ibp.ac.cn/piRBase/ 

Python 3.11
See requirements.txt for list of requirements.
FishPi was designed to be executed from the Linux command line terminal on Mac. Windows and also Mac users may use Python tools such as PyCharm Community Edition here: https://www.jetbrains.com/pycharm/download/ 

If you would like to use FishPi on a reference genome other than Zebrafish (GRCz11), Medaka (oryLat2) or Tilapia (Onil_1.2), please prepare your "teseqs.fasta" file by following the same instructions written below, but by selecting appropriate organism on UCSC browser. Also check the piRNA:TE seed rules, as these can be customised to your model organism in the python script behind the FishPi GUI. If you are working with the Zebrafish, Medaka or Tilapia please download all the annotated teseqs.fasta files here: [https://zenodo.org/records/10656843](https://zenodo.org/records/13911872) .

To prepare your GRCz11.teseqs.fasta file:
1. To generate a fasta file of the TE sequences in Zebrafish reference genome GRCz11 first make a .bed file. The UCSC table genome browser was used to generate the TE.bed file with these options: clade: Vertebrate, group: Variation and Repeats, genome: Zebrafish, assembly: May 2017 GRCz11, track: RepeatMasker, table: rmsk.
2. To extract DNA sequences from the reference genome based on the co-ordinates supplied in the bed file the following command was used:
   bedtools getfasta -s -name -fi Danio_rerio.GRCz11.dna.primary_assembly.fa -fo GRCz11.teseqs.use.fasta -bed GRCz11.teannotation.bed



Example output files are based on dre-piRNA-1 5'-TGGTTAGTACTTGGATGGGAGACCGCCTGGG-3', taken from piRBase (http://bigdata.ibp.ac.cn/piRBase/browse.php). 


# Install and use FishPi (Linux)

```
$ git clone https://github.com/alicegodden/fishpi
$ cd fishpi # navigate to FishPi directory
$ cd files # navigate to FishPi scripts and files
$ conda env create -f environment.yml # Create your conda environment
$ conda activate fishpi_environment # activate your environment

# Verify the Installation
# After activating the environment, verify that everything is installed correctly:
$ python --version
$ conda list

# To download the TE sequence fasta files, download all files
$ wget https://zenodo.org/records/13911872/files/GRCz11_ensembl_teseqs.fishpi.fasta.gz # For Zebrafish
$ wget https://zenodo.org/records/13911872/files/Onil_1.2_ensembl_teseqs.fishpi.fasta.gz # For Tilapia
$ wget https://zenodo.org/records/13911872/files/oryLat2_ensembl_teseqs.fishpi.fasta.gz # For Medaka

# Then uncompress fasta files
$ gunzip *.fasta.gz

# You should now be ready to run FishPi
$ python FishPi.py # opens the GUI to use FishPi


```
# Using FishPi with PyCharm IDE

Make sure you download all files here https://github.com/alicegodden/fishpi/tree/main/files into your working directory in Pycharm, plus make sure you have the relevant fasta files from Zenodo linked above. 


# How to Use FishPi:

Input piRNA Sequence: Launch the FishPi application. Enter your piRNA sequence in the provided text field and click the "Analyse piRNA sequence" button.
piRNA sequences can be obtained from your own analyses or from a repository like piRbase: http://bigdata.ibp.ac.cn/piRBase/ 

View Analysis Results: FishPi will analyze the piRNA sequence and display the results in a popup window. The results include a bar chart showing complementary TE counts by type and chromosomal location of complementary TEs. 

Export Data: You can export the complementary TE sequences to a CSV file and save the bar chart as an image file using the export buttons provided in the popup window, (Figure 1).



*Figure 1- Graphic User interface for FishPi and results for dre-piRNA 5'-TACACGAAGACTGTGGTGTGATTGGGCG-3'*
![FishPi_GUI](https://github.com/user-attachments/assets/a7e64def-29a4-4aae-b883-5335beca54b8)


# Flow chart of FishPi functions
![fishpi_flow_diagram](https://github.com/user-attachments/assets/d13b522a-8e11-411b-9e76-2b8fe350511f)
*Figure 2- Flow chart covering key steps and processes by FishPi*

# Output results
Transposable element species are grouped by order in the bar chart displayed, as in Figure 1 above. 
The transposable elements covered are grouped as follows:

DNA:	hAT, Tc1, Tc-Mar, Harbinger, Enspm, Kolobok, Merlin, Crypton, PiggyBac, Dada, Zatar, Ginger, TDR, Polinton, Maverick, Acrobat, Looper, TZF, Angel, Mariner

LTR:	Gypsy, DIRS, Ngaro, ERV, Pao, Copia, BEL, HERV, Bhikari

LINE:	L1, L2, L1-Tx1, Rex-Babar, RTE, Penelope, Keno, Rex

SINE:	Alu, tRNA-V-RTE

RC:	Helitron

Satellite: 	BRSATI, MOSAT

This is based off Dfam and Repeatmasker classifications.

# Licence:
This project is distributed under the GNU General Public License (GPL), which ensures that it remains open source and freely accessible for use.

# Citation:
FishPi: a bioinformatic prediction tool to link piRNA and transposable elements in zebrafish
Alice May Godden, Benjamin Rix, Simone Immler
bioRxiv 2024.09.10.612046; doi: https://doi.org/10.1101/2024.09.10.612046

# Author & Contact:
Dr Alice M. Godden

Dr. Benjamin Rix 

Prof. Simone Immler

University of East Anglia, School of Biological Sciences, Norwich, United Kingdom, NR4 7TJ

Contact:

Email: alice.godden@uea.ac.uk
Feel free to report issues or suggest improvements. FishPi is designed to make piRNA sequence analysis accessible and informative for researchers and bioinformaticians.

Note:
The "fishpi.png" image file is used for the graphical interface, if you are running FishPi in Pycharm you will need this image in order for the software to run.
