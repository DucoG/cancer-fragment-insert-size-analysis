# binwise_fragmentomics_calculator.py
# this script calculates fragmentomic features for within each bin for given bam and a given bin size
# the script has a loader function to load the bam file, and a function to calculate the bin ranges
# then a fragmentomic feature function is defined, which is called in the main function and is evaluated for each bin
# there is also a save function which saves the results in the appropriate format
# the main function calls the loader function, the bin range function, the fragmentomic feature function and the save function
# the main function also calls a argparse function to parse the arguments
# arguments are: bam_files, bin_size, contigs_to_exclude, list of fragmentomic feature functions, output_file
# fragmentomic calculater functions: mean_fragment_length, gc_content, fraction_non-proper-pairs, fraction_read_start_in_nucleosome_center, mono-, di- and tri-nucleotide frequencies at fragment start and end

import time
import pysam
import argparse
import numpy as np
import math
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam_file", type=str)
    # parser.add_argument("fasta_file", type=str)
    parser.add_argument("nucleosome_list", type=str)
    parser.add_argument("bin_size", type=int)
    parser.add_argument("contigs_to_exclude", nargs="*")
    parser.add_argument("output_file")
    return parser.parse_args()


def calculate_bin_ranges(bam_file, bin_size, contigs_to_exclude=[]):
    bam = pysam.AlignmentFile(bam_file, "rb")

    contig_iterator = zip(bam.references, bam.lengths)
    # filter contigs to exclude
    if contigs_to_exclude:
        contig_iterator = filter(
            lambda x: x[0] not in contigs_to_exclude, contig_iterator)

    bin_ranges = []
    for contig, length in contig_iterator:
        start = 0
        while start < length:
            end = min(start + bin_size, length)
            bin_ranges.append((contig, start, end))
            start = end
    return bin_ranges

# fragmentomic feature functions


def calculate_mean_fragment_length(bam_file, bin_ranges):
    bam = pysam.AlignmentFile(bam_file, "rb")

    mean_fragment_lengths = []
    for contig, start, end in bin_ranges:
        # iterate over reads in bin
        read_lengths = []
        for read in bam.fetch(contig, start, end):
            # make sure both reads are mapped
            both_mapped = not read.is_unmapped and not read.mate_is_unmapped

            # make sure reads are correctly oriented
            correctly_oriented = (read.is_forward and read.mate_is_reverse) or (
                read.is_reverse and read.mate_is_forward)

            # make sure read is primary alignment
            primary_alignment = not read.is_secondary and not read.is_supplementary and not read.is_duplicate

            if both_mapped and correctly_oriented and primary_alignment and read.is_forward:

                length_check = read.template_length > 5 and read.template_length < 1000
                if length_check:
                    read_lengths.append(read.template_length)

        # calculate mean fragment length
        if len(read_lengths) > 0:
            mean_fragment_length = sum(read_lengths)/len(read_lengths)
        else:
            mean_fragment_length = 0

        mean_fragment_lengths.append(
            [(contig, start, end), mean_fragment_length])
    return mean_fragment_lengths


def calculate_binwise_read_statistics(bam_file, bin_ranges):
    bam = pysam.AlignmentFile(bam_file, "rb")

    total_reads = []

    secondary_alignments = []
    supplementary_alignments = []
    mate_unmapped_reads = []

    incorrect_orientation_reads = []
    extreme_fragment_lengths = []

    for contig, start, end in bin_ranges:
        total_reads_bin = 0

        secondary_alignments_bin = 0
        supplementary_alignments_bin = 0
        mate_unmapped_reads_bin = 0

        extreme_fragment_lengths_bin = 0
        incorrect_orientation_reads_bin = 0

        for read in bam.fetch(contig, start, end):
            if not read.is_duplicate and not read.is_unmapped:
                total_reads_bin += 1

                if read.is_secondary:
                    secondary_alignments_bin += 1
                if read.is_supplementary:
                    supplementary_alignments_bin += 1

                # fragment length checks
                if read.mate_is_unmapped:
                    mate_unmapped_reads_bin += 1
                elif not read.is_secondary and not read.is_supplementary and read.is_forward:
                    if read.mate_is_reverse:  # correct orientation, check for fragment length
                        if read.template_length < 5 or read.template_length > 1000:
                            extreme_fragment_lengths_bin += 1
                    else:  # incorrect orientation
                        incorrect_orientation_reads_bin += 1

        total_reads.append([(contig, start, end), total_reads_bin])
        mate_unmapped_reads.append(
            [(contig, start, end), mate_unmapped_reads_bin])
        secondary_alignments.append(
            [(contig, start, end), secondary_alignments_bin])
        supplementary_alignments.append(
            [(contig, start, end), supplementary_alignments_bin])
        extreme_fragment_lengths.append(
            [(contig, start, end), extreme_fragment_lengths_bin])
        incorrect_orientation_reads.append(
            [(contig, start, end), incorrect_orientation_reads_bin])

    return total_reads, mate_unmapped_reads, secondary_alignments, supplementary_alignments, extreme_fragment_lengths, incorrect_orientation_reads


# fraction of non-proper-pairs
def calculate_fraction_non_proper_pairs(bam_file, bin_ranges):
    bam = pysam.AlignmentFile(bam_file, "rb")

    total_valid_reads = []
    fractions_non_proper_pair = []
    for contig, start, end in bin_ranges:
        # iterate over reads in bin
        non_proper_pair_reads = 0
        total = 0
        for read in bam.fetch(contig, start, end):
            if not read.is_duplicate and not read.is_secondary:
                total += 1
                if not read.is_proper_pair:
                    non_proper_pair_reads += 1

        # calculate fraction of non-proper-pairs
        if total > 0:
            frac = non_proper_pair_reads/total

        fractions_non_proper_pair.append(
            [(contig, start, end), frac])

        total_valid_reads.append(
            [(contig, start, end), total])

    return fractions_non_proper_pair, total_valid_reads


# function for calculating distance
def _find_nearest_distance(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return value - array[idx-1]
    else:
        return value - array[idx]


def _find_closest_distances(A, B):
    distances = []
    index_b = 0

    for a in A:
        while index_b < len(B) - 1 and a - B[index_b] > B[index_b + 1] - a:
            index_b += 1

        closest_distance = a - B[index_b]
        distances.append(closest_distance)

    return distances


def _read_nucleosome_center_list(nucleosome_list):
    with open(nucleosome_list) as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]

    chrs = np.array([line.split('\t')[0] for line in lines])
    hist_center = np.array([line.split('\t')[1] for line in lines])
    hist_locs = np.column_stack([chrs, hist_center]).astype(object)
    hist_locs[:, 1] = hist_locs[:, 1].astype(int)

    # convert to dictionary of chromosome: nucleosome center locations
    hist_locs = {chr: hist_locs[hist_locs[:, 0] == chr, 1]
                 for chr in np.unique(hist_locs[:, 0])}
    return hist_locs

# delta score is based on the fraction of reads starting in the nucleosome center vs starting at the end of the nucleosome


def calculate_delta_scores(bam_file, bin_ranges, nucleosome_list):
    hist_locs = _read_nucleosome_center_list(nucleosome_list)

    bam = pysam.AlignmentFile(bam_file, "rb")

    delta_scores = []
    for contig, start, end in bin_ranges:
        # mask to get only current chromosome and only nucleosomes within 10000 bp of bin
        chrom_hist_locs = hist_locs[contig]

        read_locs = []
        for read in bam.fetch(contig, start, end):
            # skip read if duplicate or secondary
            if not (read.is_duplicate or read.is_secondary or read.is_unmapped):
                # extract the read position if its reverse or not
                if read.is_reverse:
                    read_locs.append(read.reference_end)
                else:
                    read_locs.append(read.reference_start)

        # find the distance to the closest nucleosome
        read_dists = _find_closest_distances(read_locs, chrom_hist_locs)

        # calculate delta = frac_reads_in_center - frac_reads_at_end
        read_dists = np.array(read_dists)
        if len(read_dists) > 0:
            # reads in center are reads with distance -50 to 50 using numpy
            reads_in_center = len(
                read_dists[(read_dists >= -50) & (read_dists <= 50)])
            fraction_reads_center = reads_in_center/len(read_dists)

            # reads at end are reads with distance -120 to -50 and 50 to 120
            reads_at_end = len(read_dists[(read_dists >= -120) & (read_dists <= -50)]) + len(
                read_dists[(read_dists >= 50) & (read_dists <= 120)])
            fraction_reads_at_end = reads_at_end/len(read_dists)

            delta = fraction_reads_center - fraction_reads_at_end
        else:
            delta = None

        delta_scores.append([(contig, start, end), delta])

    return delta_scores


def merge_features(feature_lists, feature_names):
    # Create a dictionary to store merged features
    merged = {}

    # Iterate through each feature list
    for feature_list in feature_lists:
        # Iterate through each item in the feature list
        for item in feature_list:
            # Get the contig, start, end, and feature value
            key, feature = item

            # If the key is not in the merged dictionary, initialize it with an empty list
            if key not in merged:
                merged[key] = []

            # Append the feature value to the list
            merged[key].append(feature)

    # Convert the merged dictionary into a list of lists
    merged_list = [[contig, start, end, *features]
                   for (contig, start, end), features in merged.items()]

    # Convert the merged list into a pandas DataFrame
    df = pd.DataFrame(merged_list, columns=[
                      "contig", "start", "end"] + feature_names)

    return df


# main function
if __name__ == "__main__":
    # parse args
    args = parse_arguments()
    # print all args
    print(args)

    # calculate bin ranges
    bin_ranges = calculate_bin_ranges(
        args.bam_file, args.bin_size, args.contigs_to_exclude)

    # calculate mean fragment length
    print("Calculating mean fragment length...")
    mean_fragment_lengths = calculate_mean_fragment_length(
        args.bam_file, bin_ranges)

    # # calculate fraction of non-proper pair reads
    # print("Calculating fraction of non-proper pair reads...")
    # fractions_non_proper_pair = calculate_fraction_non_proper_pairs(
    #     args.bam_file, bin_ranges)

    # calculate read statistics
    print("Calculating read statistics...")
    total_reads, mate_unmapped_reads, secondary_alignments, supplementary_alignments, extreme_fragment_lengths, incorrect_orientation_reads = calculate_binwise_read_statistics(
        args.bam_file, bin_ranges)

    # calculate delta score
    print("Calculating delta score...")
    delta_scores = calculate_delta_scores(
        args.bam_file, bin_ranges, args.nucleosome_list)

    # merge features
    feature_lists = [mean_fragment_lengths,
                     total_reads,
                     mate_unmapped_reads,
                     secondary_alignments,
                     supplementary_alignments,
                     extreme_fragment_lengths,
                     incorrect_orientation_reads,
                     delta_scores]

    feature_names = ["mean_fragment_length",
                     "total_reads",
                     "mate_unmapped_reads",
                     "secondary_alignments",
                     "supplementary_alignments",
                     "extreme_fragment_lengths",
                     "incorrect_orientation_reads",
                     "delta_score"]

    df = merge_features(feature_lists, feature_names)

    # write to file
    df.to_csv(args.output_file, sep="\t", index=False)
