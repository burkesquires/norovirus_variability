import norovirus_variability as nv

good_descriptions = [">GU930737/Hu/GII.6/E9913646/1997/USA", "> AF414410/Hu/GII.6/Miami/292/1994/US",
                     ">AB039778/Hu/GII.6/SaitamaU16/1997/JPN", ">Hu/GII.6/Spencer/1971",
                     ">AB083781/YURI32073/JP/1997"]

bad_descriptions = [">AB078337/Hu/GII.6/Ueno7k/JPN", ">Test/Hu/GGG/04/AAA"]


def test_parse_year_from_description():

    assert nv.parse_year_from_description(good_descriptions[0], '/', -2) == 1997
    assert nv.parse_year_from_description(good_descriptions[1], '/', -2) == 1994
    assert nv.parse_year_from_description(good_descriptions[2], '/', -2) == 1997
    assert nv.parse_year_from_description(good_descriptions[3], '/', -1) == 1971
    assert nv.parse_year_from_description(good_descriptions[4], '/', -1) == 1997
    assert nv.parse_year_from_description(bad_descriptions[1], '/', -2) == 2004

    # no year found
    assert nv.parse_year_from_description(bad_descriptions[0], '/', -2) == None


def test_compute_differences():
    assert len(nv.compute_differences("tests/test_sequences.fasta", '/', -2)) == 10


def test_align_sequences():
    assert nv.align_sequences("GGII.13.faa") == "GGII.13.aligned.fasta"
    assert nv.align_sequences("tests/test_sequences.fasta") == "tests/test_sequences.aligned.fasta"


def test_get_difference():
    string_1 = 'abcdefghij'
    string_2 = 'aacceeggii'

    assert nv.get_difference(string_1, string_2) == 5
    assert nv.get_difference(str.upper(string_1), str.upper(string_2)) == 5


def test_translate_file():
    translation = [">AB084071/Hu/GII.6/GIFU/1999/JPN\n", "MKMASNDAAPSNDGAANLVPEANNEVMALEPVVGASIAA\n",
                   ">DQ093064/Hu/GII.6/445/JPN\n", "MKMASNDAAPSNDGAANLVPEANNEVMALEPVVGASIAA\n",
                   ">AF414408/Hu/GII.6/Baltimore/274/1993/USA\n", "MKMASNDAAPSNDGAANLVPEANNEVMALEPVVGASIAA\n",
                   ">AB039777/Hu/GII.6/SaitamaU4/1997/JPN\n", "MKMASNDAAPSNDGAANLVPEANNEVMALEPVVGASIAA\n",
                   ">AB039776/Hu/GII.6/SaitamaU3/1997/JPN\n", "MKMASNDAAPSNDGAANLVPEANNEVMALEPVVGASIAA\n",
                   ">AB078337/Hu/GII.6/Ueno7k/JPN\n", "MKMASNDAAPSNDGAANLVPEANDEVMALEPVVGASIAA\n",
                   ">KC576910/Hu/GII.6/S9c/1976/SEN\n", "MKMASNDAAPSNDGAANLVPEANNEVMALEPVVGASIAA\n"]

    test_translation = open(nv.translate_file("tests/test_sequences.fasta")).readlines()

    assert test_translation == translation
