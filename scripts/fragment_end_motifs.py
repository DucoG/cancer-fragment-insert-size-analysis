import argparse
import pysam
import pandas as pd
from Bio.Seq import Seq
from collections import defaultdict


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


def main(args):
    bam_file = args.bam_file
    ref_genome_file = args.ref_genome_file
    output_file = args.output_file

    bam_pointer = pysam.AlignmentFile(bam_file, "rb")
    hg19 = pysam.FastaFile(ref_genome_file)

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

            extracted_data.append({
                'read_id': read1.query_name,
                'strand': strand,
                'chrom': read1.reference_name,
                'read1_end_pos': read1_end_pos,
                'read1_seq': read1_seq,
                'read2_end_pos': read2_end_pos,
                'read2_seq': read2_seq
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
    parser.add_argument("output_file", help="Output file in parquet format")

    args = parser.parse_args()
    main(args)
