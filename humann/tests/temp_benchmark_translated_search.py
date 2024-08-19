import logging
import multiprocessing
import os.path
import unittest

from humann import config, utilities
from humann.search import translated
from humann.tests import cfg


class TestAdvancedHumannTranslatedSearchFunctions(unittest.TestCase):
    """
    Test the functions found in humann.search.translated
    """

    def setUp(self):
        config.unnamed_temp_dir ="/tmp/"

        # set default identity threshold
        config.identity_threshold = config.identity_threshold_uniref90_mode
        config.diamond_opts = config.diamond_opts_uniref90
        config.threads = multiprocessing.cpu_count()-3

        # set up nullhandler for logger
        logging.getLogger('humann.search.translated').addHandler(logging.NullHandler())
        logging.getLogger('humann.search.blastx_coverage').addHandler(logging.NullHandler())

    def test_diamond(self):
        alignment_file = utilities.unnamed_temp_file("diamond_results_")
        uniref = "./databases/uniref"
        unaligned_reads_file_fasta = os.path.join(cfg.data_folder, "large_demo.fasta")
        translated.diamond_alignment(alignment_file, uniref, unaligned_reads_file_fasta)

    def test_diamond_multistage(self):
        alignment_file = utilities.unnamed_temp_file("diamond_results_")
        uniref = "./databases/uniref2"
        unaligned_reads_file_fasta = os.path.join(cfg.data_folder, "large_demo.fasta")
        translated.diamond_multistage(alignment_file, uniref, unaligned_reads_file_fasta)
