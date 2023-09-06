import pysam
import numpy as np

def load_RPRs_as_dict(filename, fasta_file):
    print("Loading RPRs from", filename)
    RPRs = {}
    fasta = pysam.FastaFile(fasta_file)

    with open(filename, 'r') as file:
        for line in file:
            chrom, positions = line.strip().split(':')
            chrom = 'chr' + chrom
            start, end = map(int, positions.split('-'))

            if chrom not in RPRs:
                RPRs[chrom] = np.full(fasta.get_reference_length(chrom), False)

            RPRs[chrom][start:end] = True

    return RPRs

RPRs = load_RPRs_as_dict('/ssd/d.gaillard/nucleosome_list.txt', '/ssd/d.gaillard/hg19.fa')
