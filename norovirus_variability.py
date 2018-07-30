#!~/anaconda/bin/ python

__author__ = 'R. Burke Squires'
__version__ = 2.1


def align_sequences(input_file_path):
    """
    Align the sequences if not already aligned.

    :param input_file_path: the path for the sequence file needed to be aligned
    :return:
    """
    import os

    # parse file name and extensions
    #muscle_exe = r"/usr/local/bin/muscle"
    muscle_exe = r"~/anaconda/bin/muscle"

    in_file = r"%s" % input_file_path

    dirname, basename, root, ext = parse_file_path(input_file_path)
    if dirname == '':
        out_file = r"%s.aligned.fasta" % root
    else:
        out_file = r"%s/%s.aligned.fasta" % (dirname, root)

    print(out_file)
    arguments = ""
    # args = "-verbose -log ~/muscle.log"
    # args = "-verbose -log ~/muscle.log -gapopen -50 -gapextend 50"

    commandline = "%s -in %s -out %s %s" % (muscle_exe, in_file, out_file, arguments)

    os.system(commandline)
    return out_file


def translate_file(input_file_path):
    """
    Translate nucleotide sequence to protein sequence

    :param input_file_path:
    :return:
    """
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    input_handle = open(input_file_path, "rU")

    dirname, basename, root, ext = parse_file_path(input_file_path)
    output_file = '%s.faa' % root
    output_handle = open(output_file, "w")

    sequences = SeqIO.parse(input_handle, "fasta")
    protein_sequences = []
    for sequence in sequences:
        protein_sequences.append(SeqRecord(sequence.seq.translate(), id=sequence.id, name=sequence.name,
                                           description=sequence.description))
    count = SeqIO.write(protein_sequences, output_handle, "fasta")

    output_handle.close()
    input_handle.close()
    print("Coverted %i records" % count)
    return output_file


def get_difference(string_1, string_2):
    """
    Count every different between two strings

    :param string_1:
    :param string_2:
    :return:
    """
    return sum(ele_x != ele_y for ele_x, ele_y in zip(string_1, string_2))


def type_of_value(text):
    """
    If given string is an int, return the original str

    :param text:
    :return:
    """
    try:
        int(text)
        return int(text)
    except ValueError:
        pass

    return str


def parse_year_from_description(description, separator, year_position):
    """
    Parse year from virus description

    :param description: The description string that will be parsed
    :param separator: What separator is used
    :param year_position: The position is the year information
    :return:
    """
    from datetime import date

    today = date.today()
    current_year = today.year

    data = description.split(separator)
    year = type_of_value(data[int(year_position)])

    if isinstance(year, str):
        year = type_of_value(data[-1])

    if isinstance(year, int):
        if (year > 1800) and (year < current_year):
            return year
        elif (year > (current_year - 2000)) and (year < 99):
            return 1900 + year
        elif (year > 0) and (year < (current_year - 2000)):
            return 2000 + year
            # else:
            # logging.info("No date found for: %s" % description)
            # else:
            # logging.info("No date found for: %s" % description)


def alignment_to_list(alignment):
    """
    Convert alignment to list of aligned sequences

    :param alignment:
    :return:
    """

    from Bio import AlignIO

    align = AlignIO.read(alignment, "fasta")

    return list(align)


def compute_differences(alignment, separator, year_position):
    """
    Compute the difference between years and sequence differences for each pair

    :param alignment: The alignment to analyze
    :param separator: What separator is used
    :param year_position: The position is the year information
    :return:
    """
    logging.info('compute_differences - starting\n')

    data_points = []

    align_list_1 = alignment_to_list(alignment)
    align_list_2 = list(align_list_1)

    logging.info('length of list: %s' % len(align_list_1))

    for i, ele_1 in enumerate(align_list_1):

        align_list_2.remove(ele_1)

        for j, ele_2 in enumerate(align_list_2):

            if ele_1.description != ele_2.description:

                logging.info('string 1: %s: (length: %s)\n%s' % (ele_1.description, len(ele_1.seq), ele_1.seq))
                logging.info('string 2: %s: (length: %s)\n%s' % (ele_2.description, len(ele_2.seq), ele_2.seq))

                year_1 = parse_year_from_description(ele_1.name, separator, year_position)
                year_2 = parse_year_from_description(ele_2.name, separator, year_position)

                if isinstance(year_1, int) and isinstance(year_2, int):

                    sequence_1 = ele_1.seq.rstrip('-')
                    sequence_2 = ele_2.seq.rstrip('-')

                    logging.info('Stripped string 1: %s' % len(sequence_1))
                    logging.info('Stripped string 2: %s' % len(sequence_2))

                    year_diff = abs(year_1 - year_2)
                    logging.info('Year differences: %i' % year_diff)

                    residue_difference = get_difference(sequence_1, sequence_2)
                    logging.info('Differences: %i\n' % residue_difference)

                    data_points.append((year_diff, residue_difference))

    logging.info('Total number of data points: %i\n' % len(data_points))

    return data_points


def convert_list_to_dataframe(data_points):
    """

    :param data_points:
    :return:
    """
    import collections
    import math
    import pandas as pd

    df = pd.DataFrame()

    if len(data_points) > 0:

        logging.debug(data_points)
        max_year_diff, max_residue_diff = map(max, zip(*data_points))

        factor = 25
        years = math.ceil(max_year_diff / factor) * factor
        residues = math.ceil(max_residue_diff / factor) * factor

        df = pd.DataFrame(0, index=range(years), columns=range(residues))
        counter = collections.Counter(data_points)

        for coordinates, frequency in counter.most_common():
            df.ix[coordinates[0], coordinates[1]] = frequency

        return df

    else:

        return df


def save_dataframe_to_file(df, input_file_path):
    """
    Save data frame to file

    :param df: The DataFrame of data
    :param input_file_path: the input file path
    :return:
    """
    if df is not df.empty:

        if not os.path.exists('output'):
            os.makedirs('output')

        dirname, basename, root, ext = parse_file_path(input_file_path)

        output_file = 'output/%s-results.csv' % root
        df.to_csv(output_file)

        return output_file, root

    else:
        logging.error("DataFrame is empty.")

        return -1


def parse_file_path(input_file_path):
    """
    Parse a file path into components

    :param input_file_path:
    :return:
    """
    import os

    dirname, basename = os.path.split(input_file_path)
    root, ext = os.path.splitext(basename)
    return dirname, basename, root, ext


def create_heatmap(dataframe, root):
    """

    :param dataframe: the dataframe to plat data from 
    :param root:
    :return:
    """
    import matplotlib.pyplot as plt

    font = {'size': 8}
    plt.rc('font', **font)

    plt.pcolor(dataframe, cmap='Blues')
    # plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    # plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns)
    plt.xlabel('Amino Acid Difference')
    plt.ylabel('Isolation Year Difference')
    plt.colorbar()
    
    output_file = 'output/%s-results.png' % root
    plt.savefig(output_file)


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg


def main():
    """
    Compute difference between years and sequences.

    :return:
    """
    input_file = args.input_file_path
    if args.sequence_type == 'nucleotide':

        logging.info('Starting analysis. Translating FASTA file.\n')
        input_file = translate_file(input_file)

    logging.info('Aligning sequences.\n')
    aligned_file = align_sequences(input_file)

    logging.info('Computing differences.\n')
    data_points = compute_differences(aligned_file, args.separator, args.year_position)

    logging.info('Computing frequencies and converting to dataframe.')
    df = convert_list_to_dataframe(data_points)

    logging.info('Saving output.')
    output_file, root = save_dataframe_to_file(df, input_file)

    create_heatmap(df, root)


if __name__ == '__main__':
    import argparse
    import logging
    import datetime
    import os, sys

    # print(sys.argv)

    log_directory = "logs"

    d = datetime.date.today()
    # d = datetime.datetime.today()
    if not os.path.exists(log_directory):
        os.makedirs(log_directory)
    fh = logging.FileHandler('logs/%s.log' % d.isoformat())
    logging.basicConfig(filename=fh.baseFilename, level=logging.DEBUG, filemode="w")

    parser = argparse.ArgumentParser(prog='norovirus_variability.py',
                                     usage="python %(prog)s -f input_file\n",
                                     description='The norovirus variability script will plot the '
                                                 'variability of norovirus proteins.',
                                     formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=5))
    required_group = parser.add_argument_group('Required')
    required_group.add_argument("-f", '--input_file', dest="input_file_path", required=True,
                                help="The nucleotide or amino acid sequence file to be analyzed", metavar="FILE",
                                type=lambda x: is_valid_file(parser, x))

    options_group = parser.add_argument_group('Options')
    options_group.add_argument('-t', '--sequence_type', required=False, default='protein',
                               help='The format of the sequence file (options: "protein", "nucleotide").')
    options_group.add_argument("-s", '--separator', required=False, default='/',
                               help="The separator used in the description. (default '/')?")
    options_group.add_argument("-y", '--year_position', required=False, default=-2,
                               help='The position (working backwards) of the year in the description (default: -2).')

    args = parser.parse_args()

    main()