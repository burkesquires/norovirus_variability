usage: python norovirus_variability.py -f file

The norovirus variability script will plot the variability of norovirus
proteins.

optional arguments:
  -h, --help
     show this help message and exit

Required:
  -f FILE, --file FILE
     The nucleotide sequence file to be translated and analyzed.

Options:
  -t SEQUENCE_TYPE, --sequence_type SEQUENCE_TYPE
     The format of the sequence file (options: "protein", "nucleotide").
  -s SEPARATOR, --separator SEPARATOR
     The separator used in the description. (default '/')?
  -y YEAR_POSITION, --year_position YEAR_POSITION
     The position (working backwards) of the year in the description (default:
     -2).