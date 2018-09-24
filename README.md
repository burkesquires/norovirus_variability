# norovirus_variability

Gabriel I. Parra, R. Burke Squires, Consolee K. Karangwa, Jordan Johnson, Cara Lepore, Stanislav V. Sosnovtsev, and Kim Y. Green. [Static and Evolving Norovirus Genotypes: Implications for Epidemiology and Immunity](http://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1006136). PLOS Pathogens. January 19, 2017

The data for figures 1b and 2b wasgenerated with this program.

# README

norovirus_variability.py

INSTALLATION
-----------------

Installation instructions:

1. Download and install the Anaconda Python 3 distribution.

NOTE: Install for your use only in your user directory.

This installation does not require any administrative privileges.

Additionally, if you ever encounter any issues with Anaconda or python, simply drag the
anaconda3 folder or directory, that is in your home directory, to the trash and reinstall.

2. the multiple sequence alignment software muscle needs to be installed and runnable from the command line. If you have not installed muscle yet, after you have installed the Anaconda distribution you can install muscle by following the direction here:

https://bioconda.github.io/

In summary, add channels to the condo installer by typing:

$ conda config --add channels conda-forge
$ conda config --add channels defaults
$ conda config --add channels bioconda

Then install muscle:

$ conda install muscle

You can confirm that muscle is installed by typing 'which muscle' on the command line. You should get a path similiar to:

# /Users/squiresrb/anaconda/bin/muscle


3. Copy the norovirus_variability.py script to any folder in your directory

4. From the terminal or command line run the script like this:

$python python norovirus_variability.py -f ./path/to/data_file.fasta


HELP
====

To get help type: $ python norovirus_variability.py --h

You will see:

$ python norovirus_variability.py --h
usage: python norovirus_variability.py -f input_file

The norovirus variability script will plot the variability of norovirus
proteins.

optional arguments:
-h, --help
show this help message and exit

Required:
-f FILE, --input_file FILE
The nucleotide or amino acid sequence file to be analyzed

Options:
-t SEQUENCE_TYPE, --sequence_type SEQUENCE_TYPE
The format of the sequence file (options: "protein", "nucleotide").
-s SEPARATOR, --separator SEPARATOR
The separator used in the description. (default '/')?
-y YEAR_POSITION, --year_position YEAR_POSITION
The position (working backwards) of the year in the description (default:
-2).

Examples of usage:

$ python norovirus_variability.py -f ./2017_data/GII.2_FINAL2017_VP1_c_reviewed_n134_AA.fas

When running the script to create a log folder and output folder. Output is saved in the
output directory.

Cheers,

Burke Squires
