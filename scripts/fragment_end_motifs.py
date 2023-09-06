import argparse
import pysam
import pandas as pd
from Bio.Seq import Seq
from collections import defaultdict
import bisect
import pickle
import time


def read_pair_gen(bam_pointer, region_str=None):
    read_dict = defaultdict(lambda: [None, None])

    for read in bam_pointer.fetch(region=region_str):
        both_mapped = not read.is_unmapped and not read.mate_is_unmapped
        correctly_oriented = (read.is_forward and read.mate_is_reverse) or (
            read.is_reverse and read.mate_is_forward)
        primary_alignment = not read.is_secondary and not read.is_supplementary and not read.is_duplicate

        if both_mapped and correctly_oriented and primary_alignment:
            query_name = read.query_name
            if query_name not in read_dict:
                if read.is_read1:
                    read_dict[query_name][0] = read
                else:
                    read_dict[query_name][1] = read
            else:
                if read.is_read1:
                    yield read, read_dict[query_name][1]
                else:
                    yield read_dict[query_name][0], read
                del read_dict[query_name]


def get_sequence(fasta_file, reference, start, end, contig_length):
    if start < 0:
        padding = 'P' * abs(start)
        sequence = fasta_file.fetch(
            reference=reference, start=0, end=end).upper()
        return padding + sequence
    elif end > contig_length:
        padding = 'P' * abs(contig_length - end)
        sequence = fasta_file.fetch(
            reference=reference, start=start, end=contig_length).upper()
        return sequence + padding
    else:
        return fasta_file.fetch(reference=reference, start=start, end=end).upper()

def load_RPRs(filename):
    print("Loading RPRs from", filename)
    RPRs = {}
    with open(filename, 'r') as file:
        for line in file:
            chrom, positions = line.strip().split(':')
            chrom = 'chr' + chrom
            start, end = map(int, positions.split('-'))
            if chrom not in RPRs:
                RPRs[chrom] = []
            RPRs[chrom].append((start, end))
        # Sort the RPRs by start location
        for chrom in RPRs:
            RPRs[chrom].sort()
    return RPRs

def array_to_subsequent_unique(arr):
    new_arr = [arr[0]]  # Start with the first element as True
    
    for i in range(1, len(arr)):
        if arr[i] != new_arr[-1]:  # Check if the current element is different from the previous element
            new_arr.append(arr[i])  # If different, add it to the new array
    
    return new_arr

def main(args):
    bam_file = args.bam_file
    ref_genome_file = args.ref_genome_file
    RPRs_bool_file = args.RPRs_bool_file
    output_file = args.output_file

    bam_pointer = pysam.AlignmentFile(bam_file, "rb")
    hg19 = pysam.FastaFile(ref_genome_file)
    
    with open(RPRs_bool_file, 'rb') as f:
        RPRs_bool_file = pickle.load(f)
    

    print('Extracting upstream and downstream sequences from BAM file')
    # print args
    print('BAM file: {}'.format(bam_file))
    print('Reference genome file: {}'.format(ref_genome_file))
    print('Output file: {}'.format(output_file))

    contigs = bam_pointer.references
    print('contigs: {}'.format(contigs))
    first_chunk = True

    for contig in contigs:
        extracted_data = []
        print('processing contig: {}'.format(contig))
        if contig in RPRs_bool_file:
            contig_RPRs = RPRs_bool_file[contig]
        else:
            contig_RPRs = None
        contig_length = bam_pointer.get_reference_length(contig)
        for read1, read2 in read_pair_gen(bam_pointer, region_str=contig):
            if read1.is_forward:  # read1 is upstream, handle as normal
                strand = '+'
                read1_end_pos = read1.reference_start
                read1_seq = get_sequence(
                    hg19, read1.reference_name, read1_end_pos - 10, read1_end_pos + 10, contig_length)

                read2_end_pos = read2.reference_end - 1
                read2_seq = get_sequence(
                    hg19, read2.reference_name, read2_end_pos - 10, read2_end_pos + 10, contig_length)

            else:  # read1 is downstream, flip everything so that upstream bases are not sequenced
                strand = '-'
                read1_end_pos = read1.reference_end - 1
                read1_seq = get_sequence(
                    hg19, read1.reference_name, read1_end_pos - 10, read1_end_pos + 10, contig_length)
                read1_seq = str(Seq(read1_seq).reverse_complement())

                read2_end_pos = read2.reference_start
                read2_seq = get_sequence(
                    hg19, read2.reference_name, read2_end_pos - 10, read2_end_pos + 10, contig_length)
                read2_seq = str(Seq(read2_seq).reverse_complement())

            # get fragment start and end positions
            fragment_start = min(read1_end_pos, read2_end_pos)
            fragment_end = max(read1_end_pos, read2_end_pos)

            # calculate fragment length
            fragment_length = abs(fragment_end - fragment_start)
            
            # calculate GC content
            fragment_seq = get_sequence(
                hg19, read1.reference_name, fragment_start, fragment_end, contig_length)
            if fragment_start == fragment_end:
                gc_content = 0
            else:
                gc_content = (fragment_seq.count('G') + fragment_seq.count('C')) / len(fragment_seq)

            # determine if fragment ends intersect an RPR
            if contig_RPRs is None or fragment_length < 5:
                read1_intersect = False
                read2_intersect = False
                RPR_overlap = False
            else:
                read1_intersect = contig_RPRs[read1_end_pos]
                read2_intersect = contig_RPRs[read2_end_pos]
                
                # Check if the fragment overlaps an RPR
                subsequent_unique_RPRs = array_to_subsequent_unique(contig_RPRs[fragment_start:fragment_end])
                # Determine if sub-array [False, True, False] is in array
                sub_array = [False, True, False]
                overlaps_RPR = any(subsequent_unique_RPRs[i:i+len(sub_array)] == sub_array for i in range(len(subsequent_unique_RPRs) - len(sub_array) + 1))
                

            extracted_data.append({
                'read_id': read1.query_name,
                'strand': strand,
                'chrom': read1.reference_name,
                'read1_end_pos': read1_end_pos,
                'read1_seq': read1_seq, 
                'read2_end_pos': read2_end_pos,
                'read2_seq': read2_seq,
                'fragment_length': fragment_length,
                'gc_content': gc_content,
                'read1_intersect': read1_intersect,
                'read2_intersect': read2_intersect,
                'RPR_overlap': overlaps_RPR
            })

        print('converting {} to dataframe...'.format(contig))
        extracted_data = pd.DataFrame(extracted_data)
        if first_chunk:
            print('creating new parquet file: {}...'.format(output_file))
            extracted_data.to_parquet(output_file, engine='fastparquet')
            first_chunk = False
        else:
            print('adding {} to parquet file...'.format(contig))
            extracted_data.to_parquet(
                output_file, engine='fastparquet', append=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract upstream and downstream sequences from BAM file")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("ref_genome_file",
                        help="Reference genome file in FASTA format")
    parser.add_argument("RPRs_bool_file", help='pickle file containing the dict of True/False arrays')
    parser.add_argument("output_file", help="Output file in parquet format")

    args = parser.parse_args()
    tic = time.time()
    main(args)
    toc = time.time()
    print('Elapsed time: {} minutes'.format((toc - tic) / 60))
