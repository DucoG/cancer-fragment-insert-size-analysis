import pysam
import numpy as np
import math
import click
import pathlib
import json


# function for calculating distance
def find_nearest_distance(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return value - array[idx-1]
    else:
        return value - array[idx]


def read_nucleosome_center_list(nucleosome_list):
    with open(nucleosome_list) as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]

    chrs = np.array([line.split('\t')[0] for line in lines])
    hist_center = np.array([line.split('\t')[1] for line in lines])
    hist_locs = np.column_stack([chrs, hist_center]).astype(object)
    hist_locs[:, 1] = hist_locs[:, 1].astype(int)
    return hist_locs


@click.command()
@click.argument('bamfile_path', type=click.Path(exists=True, path_type=pathlib.Path, dir_okay=False))
@click.argument('output_file', type=click.Path(exists=False, writable=True, path_type=pathlib.Path, dir_okay=False))
@click.option('--nucleosome_list', required=True, type=click.Path(exists=True), help='path to list containing nucleosome locations on contigs in a tsv format')
@click.option('-u', '--unwanted_chrs', multiple=True, type=str, default=['chrM', 'chrX', 'chrY'], help='chromosomes to leave out in the analysis. Structure to be used: chr1')
def extract_nucleosome_distance_counts(nucleosome_list, bamfile_path, output_file, unwanted_chrs):
    """Extracts distance to closest nucleosome for each read from a bamfile.

    Args:

        bamfile_path (str): path of input file. If it doesnt end with .bam, it is assumed to be a text file containing the path for a bam file on each line.

        output_file (str): path to output the resulting count array

        nulceosome_list (str): path to list of nucleosome locations per contig in tsv format

        unwanted_chrs (str): unwanted chr
    """
    click.echo(f'samfile path: {bamfile_path}')
    if not bamfile_path.suffix == '.bam':
        click.echo(
            f'batch mode detected as bamfile_path is not a bam file \nAssuming {bamfile_path} contains a list of bamfiles')
        with open(bamfile_path) as file:
            bamfile_path = file.readlines()
            bamfile_path = [line.rstrip() for line in bamfile_path]
    else:
        bamfile_path = [bamfile_path]
    click.echo(f'output path: {output_file}')
    click.echo(f'nucleosome_list: {nucleosome_list}')
    click.echo(f'unwanted_chrs: {unwanted_chrs}')
    # generating nucleosome map array
    click.echo(f'reading nucleosome list...')
    hist_locs = read_nucleosome_center_list(nucleosome_list)
    click.echo(
        f'nucleosome list with {hist_locs.shape[0]} nucleosome locations loaded')

    # extract read locations from bam files
    for sample in bamfile_path:
        sample = pathlib.Path(sample)
        click.echo(f'starting sample {sample.name}')
        samfile = pysam.AlignmentFile(sample, "rb")

        # select chroms to do
        samfile_chroms = list(samfile.references)
        if not set(unwanted_chrs).issubset(set(samfile_chroms)):
            raise Exception(
                'not all unwanted chromosomes are in the samfile reference list')
        samfile_chroms = [x for x in samfile_chroms if x not in unwanted_chrs]

        # for each distance to a nucleosome (from -600 to 600), a dictionary is stored. the dictionary contains the counts per read flag for each distance
        dist_to_closest_nucleosome = {}
        for i in range(-600, 601):
            dist_to_closest_nucleosome[i] = {}

        # for each chromosome in the samfile
        for chrom in samfile_chroms:
            print(f'processing chromosome {chrom}')
            # mask to get only current chromosome
            mask = hist_locs[:, 0] == chrom
            chr_hist_locs = hist_locs[:, 1][mask]

            # make sure that it is sorted
            if not (all(chr_hist_locs[i] <= chr_hist_locs[i+1] for i in range(len(chr_hist_locs) - 1))):
                raise Exception('nucleosome list not sorted')

            for read in samfile.fetch(chrom):
                # skip read if unmapped
                if read.is_unmapped:
                    continue

                # extract the read position if its reverse or not
                if read.is_reverse:
                    read_pos = read.reference_end
                else:
                    read_pos = read.reference_start

                # find the distance to the closest nucleosome
                dist = find_nearest_distance(chr_hist_locs, read_pos)

                # add the read to the dictionary
                if dist in dist_to_closest_nucleosome.keys():
                    # add to the data
                    dist_to_closest_nucleosome[dist][read.flag] = dist_to_closest_nucleosome[dist].get(
                        read.flag, 0) + 1

        # save the sample to a json file
        with open(pathlib.Path(output_file), 'w') as f:
            json.dump(dist_to_closest_nucleosome, f)


if __name__ == '__main__':
    extract_nucleosome_distance_counts()
