#! /usr/bin/env python

"""
This is a HUMAnN utility function
* Do a first pass on blastx output
* Identify proteins that were well-covered by reads
* Return them as a dict
* When processing blastx output in HUMAnN, consider only these proteins
===============================================
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""
import bisect
import sys
import re
import logging
import argparse
from collections import defaultdict

from humann import config
from humann import utilities
from humann import store

# name global logging instance
logger = logging.getLogger(__name__)


def _get_coverage(start_positions, end_positions):
    coverage = 0
    for start, end in zip(start_positions, end_positions):
        coverage += len(range(start, end))
    return coverage


def _insert_read(start_positions, end_positions, start, end):
    if len(start_positions) == 0:
        start_positions.append(start)
        end_positions.append(end)
        return

    # insert the read to the left of anything that might be there
    # to make it safe and guarantee overlaps can be found if they exist
    insertion_index = bisect.bisect_left(start_positions, start)
    start_positions.insert(insertion_index, start)
    end_positions.insert(insertion_index, end)

    # merge all the overlapping contigs
    i = 1
    while i < len(start_positions):
        has_overlap = (start_positions[i] <= end_positions[i - 1] and
                       end_positions[i] >= start_positions[i - 1])
        if has_overlap:
            start_positions[i] = min(start_positions[i - 1], start_positions[i])
            end_positions[i] = max(end_positions[i - 1], end_positions[i])
            del start_positions[i - 1]
            del end_positions[i - 1]
            continue
        i += 1


def blastx_coverage_stream(min_coverage, log_messages=None, nucleotide=False,
                           query_coverage_threshold=config.translated_query_coverage_threshold,
                           identity_threshold=config.nucleotide_identity_threshold):
    allowed_set = set()
    protein_hits = dict()
    proteins_without_length = set()

    metadata = {
        "no_coverage_count": 0,
        "small_identity_count": 0,
        "small_query_coverage_count": 0,
    }

    def feed(protein_name, gene_length, identity, subject_start_index, subject_stop_index):
        filter = False
        if identity < identity_threshold:
            filter = True
            metadata["small_identity_count"] += 1

        # DISABLED in the other blastx coverage filter aswell, because the query_length is not known
        # if utilities.filter_based_on_query_coverage(query_length, subject_start_index, subject_stop_index,
        #                                            query_coverage_threshold):
        #    filter = True
        #    metadata["small_query_coverage_count"] += 1

        if filter:
            return False

        if subject_stop_index <= subject_start_index:
            metadata["no_coverage_count"] += 1
            return False

        if protein_name in allowed_set:
            return True

        if not nucleotide:
            gene_length /= 3

        if gene_length == 0:
            proteins_without_length.add(protein_name)
            if min_coverage <= 0:
                allowed_set.add(protein_name)
                return True
            else:
                return False

        if min_coverage <= 0:
            allowed_set.add(protein_name)
            return True

        if protein_name not in protein_hits:
            s, e = protein_hits[protein_name] = [subject_start_index], [subject_stop_index]
        else:
            s, e = protein_hits[protein_name]
            _insert_read(s, e, subject_start_index, subject_stop_index)

        coverage = _get_coverage(s, e) / float(gene_length) * 100
        if coverage >= min_coverage:
            allowed_set.add(protein_name)
            del protein_hits[protein_name]
            return True
        return False

    def finish():
        output_messages = [
            "Total alignments where percent identity is not a number: 0",
            "Total alignments where alignment length is not a number: 0",
            "Total alignments where E-value is not a number: 0",
            "Total alignments not included based on large e-value: 0",
            f"Total alignments not included based on small percent identity: {metadata["small_identity_count"]}",
            f"Total alignments not included based on small query coverage: {metadata["small_query_coverage_count"]}",
            f"Total alignments without coverage information: {metadata["no_coverage_count"]}",
            f"Total proteins in blastx output: {len(allowed_set) + len(protein_hits)}",
            f"Total proteins without lengths: {len(proteins_without_length)}",
            f"Proteins with coverage greater than threshold ({min_coverage}): {len(allowed_set)}",
        ]

        # write out informational messages to log or stdout, depending on input parameters
        if log_messages:
            for message in output_messages:
                logger.info(message)
        else:
            print("\n".join(output_messages))

    return allowed_set, feed, finish


def blastx_coverage(blast6out, min_coverage, alignments=None, log_messages=None, apply_filter=None, nucleotide=False,
                    query_coverage_threshold=config.translated_query_coverage_threshold,
                    identity_threshold=config.nucleotide_identity_threshold):
    # create alignments instance if none is passed
    if alignments is None:
        alignments = store.Alignments()

    # store protein lengths
    protein_lengths = {}
    # store unique positions hit in each protein as sets
    protein_hits = defaultdict(str)
    # track proteins with sufficient coverage
    allowed = set()
    # track alignments unable to compute coverage
    no_coverage = 0
    # parse blast6out file, applying filtering as selected
    for alignment_info in utilities.get_filtered_translated_alignments(blast6out, alignments, apply_filter=apply_filter,
                                                                       log_filter=log_messages,
                                                                       query_coverage_threshold=query_coverage_threshold,
                                                                       identity_threshold=identity_threshold):
        (protein_name, gene_length, queryid, matches, bug, alignment_length,
         subject_start_index, subject_stop_index) = alignment_info

        # divide the gene length by 3 to get protein length from nucleotide length
        if not nucleotide:
            gene_length = gene_length / 3

        # store the protein length
        protein_lengths[protein_name] = gene_length

        # add the range of the alignment to the protein hits
        protein_range = range(subject_start_index, subject_stop_index)
        if protein_range:
            # keep track of unique hit positions in this protein
            protein_hits[protein_name] += "{0}-{1};".format(subject_start_index, subject_stop_index)
        else:
            no_coverage += 1
    # track proteins without lengths
    no_length = 0
    # compute coverage
    for protein_name, hit_positions in protein_hits.items():

        # compile the hit positions
        range_hit_positions = set()
        for alignment_hit in hit_positions.split(";")[:-1]:
            start_index, stop_index = alignment_hit.split("-")
            new_range = range(int(start_index), int(stop_index))
            range_hit_positions.update(new_range)

        try:
            # compute coverage, with 50 indicating that 50% of the protein is covered
            coverage = len(range_hit_positions) / float(protein_lengths[protein_name]) * 100
        except ZeroDivisionError:
            coverage = 0
            no_length += 1

        if coverage >= min_coverage:
            allowed.add(protein_name)

    output_messages = ["Total alignments without coverage information: " + str(no_coverage)]
    output_messages += ["Total proteins in blastx output: " + str(len(protein_lengths))]
    output_messages += ["Total proteins without lengths: " + str(no_length)]
    output_messages += [
        "Proteins with coverage greater than threshold (" + str(min_coverage) + "): " + str(len(allowed))]

    # write out informational messages to log or stdout, depending on input parameters
    if log_messages:
        for message in output_messages:
            logger.info(message)
    else:
        print("\n".join(output_messages))

    return allowed


def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description="Compute blastx coverage\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i", "--input",
        help="the blastx formatted input file\n",
        required=True)
    parser.add_argument(
        "--coverage-threshold",
        type=float,
        help="the subject coverage threshold\n[ DEFAULT : " + str(config.translated_subject_coverage_threshold) + " ]",
        default=config.translated_subject_coverage_threshold)
    parser.add_argument(
        "--print-protein-list",
        action="store_true",
        help="print the list of proteins that meet the coverage threshold")

    return parser.parse_args()


def main():
    # parse the arguments from the user
    args = parse_arguments(sys.argv)

    # run coverage computation
    allowed = blastx_coverage(args.input, args.coverage_threshold)

    if args.print_protein_list:
        print("\n".join(allowed))


if __name__ == "__main__":
    main()
