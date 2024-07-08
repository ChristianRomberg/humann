import unittest
import logging

import cfg
import utils
import tempfile
import re

from humann.search import blastx_coverage
from humann import store, utilities
from humann import config


class TestBasicHumannBlastx_Coverage_StreamFunctions(unittest.TestCase):
    """
    Test the functions found in humann.search.nucleotide
    """

    def setUp(self):
        config.unnamed_temp_dir = tempfile.gettempdir()
        config.temp_dir = tempfile.gettempdir()
        config.file_basename = "HUMAnN_test"

        # set up nullhandler for logger
        logging.getLogger('humann.search.blastx_coverage').addHandler(logging.NullHandler())

    def test_blastx_coverage_gene_names_default(self):
        """
        Test the blastx_coverage function
        Test the gene names
        Test without filter
        """

        # set the coverage threshold to zero so as to not test with filter on
        current_coverage_threshold = config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold = 0
        config.identity_threshold = 0

        allowed_set, coverage_feed, coverage_finish = blastx_coverage.blastx_coverage_stream(
            config.translated_subject_coverage_threshold, log_messages=True)

        all_proteins = set()
        alignments = store.Alignments()

        for protein_name, gene_length, queryid, matches, bug, alignment_length, subject_start_index, subject_stop_index \
                in utilities.get_filtered_translated_alignments(cfg.rapsearch2_output_file_without_header_coverage,
                                                                alignments):
            identity = matches / alignment_length * 100.0

            coverage_feed(protein_name, gene_length, identity, subject_start_index, subject_stop_index)
            all_proteins.add(protein_name)

        coverage_finish()

        # reset the coverage threshold
        config.translated_subject_coverage_threshold = current_coverage_threshold

        # check the expected proteins are found
        self.assertEqual(sorted(all_proteins), sorted(allowed_set))

    def test_blastx_coverage_gene_names_custom_annotation(self):
        """
        Test the blastx_coverage function
        Test the gene names with custom annotation
        Test without filter
        """

        # set the coverage threshold to zero so as to not test with filter on
        current_coverage_threshold = config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold = 0
        config.identity_threshold = 0

        allowed_set, coverage_feed, coverage_finish = blastx_coverage.blastx_coverage_stream(
            config.translated_subject_coverage_threshold, log_messages=True)

        all_proteins = set()
        alignments = store.Alignments()

        for protein_name, gene_length, queryid, matches, bug, alignment_length, subject_start_index, subject_stop_index \
                in utilities.get_filtered_translated_alignments(cfg.rapsearch2_output_file_without_header_coverage,
                                                                alignments):
            identity = matches / alignment_length * 100.0

            coverage_feed(protein_name, gene_length, identity, subject_start_index, subject_stop_index)
            all_proteins.add(protein_name)

        coverage_finish()

        # reset the coverage threshold
        config.translated_subject_coverage_threshold = current_coverage_threshold

        # check the expected proteins are found
        self.assertEqual(sorted(all_proteins), sorted(allowed_set))


    def test_blastx_coverage_gene_names_chocophlan_annoation(self):
        """
        Test the blastx_coverage function
        Test the gene names with chocophlan annotations
        Test without filter
        """

        # set the coverage threshold to zero so as to not test with filter on
        current_coverage_threshold = config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold = 0
        config.identity_threshold = 0

        allowed_set, coverage_feed, coverage_finish = blastx_coverage.blastx_coverage_stream(
            config.translated_subject_coverage_threshold, log_messages=True)

        all_proteins = set()
        alignments = store.Alignments()

        for protein_name, gene_length, queryid, matches, bug, alignment_length, subject_start_index, subject_stop_index \
                in utilities.get_filtered_translated_alignments(cfg.rapsearch2_output_file_without_header_coverage_chocophlan_annotations,
                                                                alignments):
            identity = matches / alignment_length * 100.0

            coverage_feed(protein_name, gene_length, identity, subject_start_index, subject_stop_index)
            all_proteins.add(protein_name)

        coverage_finish()

        # reset the coverage threshold
        config.translated_subject_coverage_threshold = current_coverage_threshold

        # check the expected proteins are found
        self.assertEqual(sorted(all_proteins), sorted(allowed_set))


    def test_blastx_coverage(self):
        """
        Test the coverage filter
        Test with one protein with one alignment passing threshold
        Test with one protein with two alignments passing threshold (does not pass with only one alignment)
        Test with other proteins with one more more alignments not passing threshold
        """

        # set the coverage threshold to a small value so as to have some alignments pass
        current_coverage_threshold = config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold = 50.0
        config.identity_threshold = 0
        found_proteins = set()
        alignments = store.Alignments()

        allowed_set, coverage_feed, coverage_finish = blastx_coverage.blastx_coverage_stream(
            config.translated_subject_coverage_threshold, log_messages=True)


        for protein_name, gene_length, queryid, matches, bug, alignment_length, subject_start_index, subject_stop_index \
                in utilities.get_filtered_translated_alignments(cfg.rapsearch2_output_file_without_header_coverage,
                                                                alignments):
            identity = matches / alignment_length * 100.0

            coverage_feed(protein_name, gene_length, identity, subject_start_index, subject_stop_index)

            if "_coverage50" in protein_name:
                found_proteins.add(protein_name)

        coverage_finish()

        # reset the coverage threshold
        config.translated_subject_coverage_threshold = current_coverage_threshold

        # check the expected proteins are found
        self.assertEqual(sorted(found_proteins), sorted(allowed_set))

