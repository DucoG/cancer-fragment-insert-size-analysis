# this script calculates the delfi ratio of short reads (100<x<150) to long reads (150<x<220) for each sample for each bin for a given binsize
# to do this, for each bin, each fragment it comes across, it checks its length and adds it to the appropriate bin: 0-9, 10-19, 20-29, ...., 499, =>500
# it outputs a pandas dataframe with columns contig, start, end, GC_content, short_reads, long_reads, very_short_reads, very_long_reads

import pandas as pd
import numpy as np
import pysam
import argparse
import json


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam_file", type=str)
    parser.add_argument("fasta_file", type=str)
    parser.add_argument("bin_size", type=int)
    parser.add_argument("delfi_output_file", type=str)
    parser.add_argument("fragment_length_distribution_output_file", type=str)
    parser.add_argument("contigs_to_exclude", nargs="*")
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

# calculate the GC content of a given sequence


def get_gc_content(fasta_file, bin_ranges):
    fasta = pysam.FastaFile(fasta_file)
    gc_content = []
    for contig, start, end in bin_ranges:
        print(contig, start, end)
        seq = fasta.fetch(contig, start, end)
        gc = (seq.count("G") + seq.count("C"))/len(seq)
        gc_content.append(
            (contig, start, end, gc))
    return gc_content

# calculate binwise fragment length distribution


def get_fragment_length_distribution(bam_file, bin_ranges):
    bam = pysam.AlignmentFile(bam_file, "rb")
    fragment_length_distribution = []

    bin_edges = np.arange(0, 501, 10)
    bin_edges = np.append(bin_edges, np.inf)

    for contig, start, end in bin_ranges:
        fragment_lengths = []
        for read in bam.fetch(contig, start, end):
            # make sure both reads are mapped
            both_mapped = not read.is_unmapped and not read.mate_is_unmapped

            # make sure reads are correctly oriented
            correctly_oriented = (read.is_forward and read.mate_is_reverse) or (
                read.is_reverse and read.mate_is_forward)

            # make sure read is primary alignment
            primary_alignment = not read.is_secondary and not read.is_supplementary and not read.is_duplicate

            if both_mapped and correctly_oriented and primary_alignment and read.is_forward:
                fragment_lengths.append(abs(read.template_length))

        # calculate the binwise fragment length distribution
        bins = np.digitize(fragment_lengths, bin_edges) - 1
        bin_counts = np.bincount(bins)

        # add the binwise fragment length distribution to the list
        fragment_length_distribution.append(
            (contig, start, end, bin_counts.tolist()))

    return fragment_length_distribution

# calculate the delfi ratio for each bin


def calculate_delfi_ratios(fragment_length_distribution, fasta_file):
    delfi_ratios = []

    # Get only the (contig, start, end) from fragment_length_distribution
    bin_ranges = [(contig, start, end)
                  for contig, start, end, _ in fragment_length_distribution]
    gc_content = get_gc_content(fasta_file, bin_ranges)

    for (contig, start, end, fragment_lengths), (_, _, _, gc) in zip(fragment_length_distribution, gc_content):

        # calculate the delfi ratios
        short_reads = sum(fragment_lengths[10:15])
        long_reads = sum(fragment_lengths[15:22])
        very_short_reads = sum(fragment_lengths[:10])
        very_long_reads = sum(fragment_lengths[22:])
        total_reads = sum(fragment_lengths)

        # Append the delfi ratios as a dictionary to the list
        delfi_ratios.append({
            "contig": contig,
            "start": start,
            "end": end,
            "short_reads": short_reads,
            "long_reads": long_reads,
            "total_reads": total_reads,
            "very_short_reads": very_short_reads,
            "very_long_reads": very_long_reads,
            "GC_content": gc
        })

    delfi_ratios_df = pd.DataFrame(delfi_ratios)

    return delfi_ratios_df


# main function
if __name__ == "__main__":
    # parse args
    args = parse_arguments()
    # print all args
    print(args)

    # calculate bin ranges
    print("Calculating bin ranges...")
    bin_ranges = calculate_bin_ranges(
        args.bam_file, args.bin_size, args.contigs_to_exclude)

    # calculate fragment length distribution
    print("Calculating fragment length distribution...")
    fragment_length_distribution = get_fragment_length_distribution(
        args.bam_file, bin_ranges)

    # calculate delfi ratios
    print("Calculating delfi ratios...")
    delfi_ratio_df = calculate_delfi_ratios(
        fragment_length_distribution, args.fasta_file)

    # save delfi ratios as a tsv file
    print("Saving delfi ratios...")
    delfi_ratio_df.to_csv(args.delfi_output_file, sep='\t', index=False)

    # create a json file for the fragment length distribution
    print("Saving fragment length distribution...")
    json_data = []
    for contig, start, end, measurements in fragment_length_distribution:
        entry = {
            "contig": contig,
            "start": start,
            "end": end,
            "measurements": measurements
        }
        json_data.append(entry)

    # Save json_data as a JSON file
    with open(args.fragment_length_distribution_output_file, "w") as file:
        json.dump(json_data, file)
