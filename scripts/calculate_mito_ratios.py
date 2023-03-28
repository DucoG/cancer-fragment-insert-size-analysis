import pysam
import sys
import pandas as pd


def get_mito_contig_name(bam):
    # Check for common mitochondrial contig names
    possible_names = ["chrM", "MT", "chrMT", "M", "ChrM",
                      "chrm", "mtDNA", "ChrMT", "chr-m", "chr_mt", "chrMito"]

    for name in possible_names:
        if name in bam.references:
            return name

    # If not found, return None
    return None


def calculate_mito_ratio(bam_file):
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        total_mapped_reads = 0
        mito_mapped_reads = 0

        mito_contig_name = get_mito_contig_name(bam)

        if mito_contig_name is None:
            print("Mitochondrial contig not found.")
            sys.exit(1)

        # Iterate through the mapped reads
        for read in bam.fetch():
            if not read.is_unmapped:
                total_mapped_reads += 1

                # Check if the read is mapped to the mitochondrial contig
                if bam.get_reference_name(read.reference_id) == mito_contig_name:
                    mito_mapped_reads += 1

    if total_mapped_reads == 0:
        return 0

    # Calculate and return the ratio
    return mito_mapped_reads / total_mapped_reads


if __name__ == "__main__":
    # take in a list of bam files, calculate the mito ratio for each, and save the results to a csv file using pandas
    output_file = sys.argv[1]
    bam_files = sys.argv[2:]
    mito_ratios = []
    sample_names = []
    for bam_file in bam_files:
        sample_names.append(bam_file.split('/')[-1].split('.')[0])
        mito_ratios.append(calculate_mito_ratio(bam_file))
    df = pd.DataFrame({'sample': sample_names, 'mito_ratio': mito_ratios})
    df.to_csv(output_file, index=False)
