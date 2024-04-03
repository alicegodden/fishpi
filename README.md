![GitHub Logo](https://github.com/alicegodden/fishpi/blob/main/fishpi.png)
# FishPi

Created by Dr. Alice M. Godden & Dr. Benjamin Rix, 2023

# Description:
FishPi is a Python-based piRNA sequence analysis tool designed to identify complementary transposable element (TE) sequences in a given piRNA sequence. This project combines the power of various Python libraries, including tkinter for the graphic user interface, matplotlib for data visualization, and PIL (Pillow) for image rendering. All results and bar charts can be exported for downstream analyses and publication in hi-resolution.


# Features:

Analyze piRNA sequences for complementary TE matches.
Classify TE types based on keywords.
Generate a bar chart showing complementary TE counts by type.
Export the list of complementary TE sequences to a CSV file.
Export the bar chart as an image file.

# Requirements and input file preparation:

Python 3.11
Python libraries: tkinter, matplotlib, Pillow (PIL)
FishPi was designed to be executed from the Linux command line terminal on Mac. Windows and also Mac users may use Python tools such as PyCharm Community Edition here: https://www.jetbrains.com/pycharm/download/ 

If you would like to use FishPi on a Zebrafish reference genome other than GRCz11, please prepare your "teseqs.fasta" file by following the same instructions written below, but by selecting appropriate organism on UCSC browser. Also check the piRNA:TE seed rules, as these can be customised to your model organism in the python script behind the FishPi GUI. If you are working with the Zebrafish reference GRCz11 please download the annotated GRCz11.teseqs.fasta file here: https://zenodo.org/records/10656843 .

To prepare your GRCz11.teseqs.fasta file:
1. To generate a fasta file of the TE sequences in Zebrafish reference genome GRCz11 first make a .bed file. The UCSC table genome browser was used to generate the TE.bed file with these options: clade: Vertebrate, group: Variation and Repeats, genome: Zebrafish, assembly: May 2017 GRCz11, track: RepeatMasker, table: rmsk.
2. To extract DNA sequences from the reference genome based on the co-ordinates supplied in the bed file the following command was used: bedtools getfasta -s -name -fi Danio_rerio.GRCz11.dna.primary_assembly.fa -fo GRCz11.teseqs.use.fasta -bed GRCz11.teannotation.bed.
3. To mask the sequences in the reference genome the following script was used: bedtools maskfasta -fi Danio_rerio.GRCz11.dna.primary_assembly.fa -fo GRCz11.masked.fasta -bed GRCz11.teannotation.bed.
4. To create the TE-merged reference genome: cat GRCz11.masked.fasta GRCz11.teseqs.1.fasta > GRCz11.teseqs.fishpi.fasta (GRCz11.teseqs.use.fasta). 


Clone or download the FishPi repository to your local machine, download all files here: https://github.com/alicegodden/fishpi/releases/tag/FishPi  
You need to make sure you have the GRCz11.teseqs.fasta, fishpi.png and FishPi.py in the area you execute your script, or same directory as you are running your project if using PyCharm.
Ensure you have Python v3.11 and the required libraries installed.
Prepare the TE sequence file and place it in the same directory as FishPi (only relevant if you want to use TE's not in GRCz11 annotation).
Fishpi.py could be run on the command line as described below on Linux, or with Pycharm.
Example output files are based on dre-piRNA-1 5'-TGGTTAGTACTTGGATGGGAGACCGCCTGGG-3', taken from piRBase (http://bigdata.ibp.ac.cn/piRBase/browse.php). 


# Install and use FishPi (Linux)

```
$ git clone github.com/alicegodden/fishpi
$ cd fishpi # navigate to FishPi directory
$ python fishpi.py # opens the GUI to use FishPi


```


# How to Use FishPi:

Input piRNA Sequence: Launch the FishPi application. Enter your piRNA sequence in the provided text field and click the "Submit analysis" button.
piRNA sequences can be obtained from your own analyses or from a repository like piRbase: http://bigdata.ibp.ac.cn/piRBase/ 

View Analysis Results: FishPi will analyze the piRNA sequence and display the results in a popup window. The results include a bar chart showing complementary TE counts by type and the total TE count.

Export Data: You can export the complementary TE sequences to a CSV file and save the bar chart as an image file using the export buttons provided in the popup window, (Figure 1).


*Figure 1- Graphic User interface for FishPi and results for dre-piRNA-1*
![Screenshot 2023-12-01 at 12 39 58](https://github.com/alicegodden/fishpi/assets/136358959/f122ffd1-7791-4957-b834-88072fbc7c11)


# Output results
Transposable element species are grouped by order in the bar chart displayed, as in Figure 1 above. 
The transposable elements covered are grouped as follows:

DNA:	hAT, Tc1, Tc-Mar, Harbinger, Enspm, Kolobok, Merlin, Crypton, PiggyBac, Dada, Zatar, Ginger, TDR, Polinton, Maverick, Acrobat, Looper, TZF, Angel, Mariner

LTR:	Gypsy, DIRS, Ngaro, ERV, Pao, Copia, BEL, HERV, Bhikari

LINE:	L1, L2, L1-Tx1, Rex-Babar, RTE, Penelope, Keno, Rex

SINE:	Alu, tRNA-V-RTE

RC:	Helitron

Satellite: 	BRSATI, MOSAT


# Licence:
This project is distributed under the GNU General Public License (GPL), which ensures that it remains open source and freely accessible for use.

# Citation:
Publication of Fishpi to follow, please check back for full citation or contact authors, details below. In the meantime please reference this git repository.

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
