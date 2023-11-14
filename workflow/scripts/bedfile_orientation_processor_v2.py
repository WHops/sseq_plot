import sys
import gzip
from collections import defaultdict

def flip_orientation(orientation):
    """
    Flips the orientation of a strand.

    Args:
        orientation (str): Strand orientation ('+' or '-').

    Returns:
        str: Flipped strand orientation.
    """
    return '-' if orientation == '+' else '+'

def process_cell(cell, orientation_count, outfile, cell_name):
    """
    Processes a cell's entries, flips orientations, appends cell name, and writes to the output file.

    Args:
        cell (list): List of entries representing a cell.
        orientation_count (defaultdict): Count of orientations in the cell.
        outfile: Output file object for writing.
        cell_name (str): Name of the cell.
    """
    # Check if the ratio between - and + counts is at least 25:75 in either direction. If not, set is_wc to True
    majority_orientation_count = orientation_count['-'] if orientation_count['-'] > orientation_count['+'] else orientation_count['+']
    minority_orientation_count = orientation_count['-'] if orientation_count['-'] < orientation_count['+'] else orientation_count['+']
    print(majority_orientation_count, minority_orientation_count)
    if (minority_orientation_count == 0):
        print('Zero reads found; adding a virtual one')
        minority_orientation_count += 0.1

    
    if majority_orientation_count / minority_orientation_count < 3:
        print("Looks like a WC cell! Skipping...")
        return()
    print('This is a proper cell')
    majority_orientation = '-' if orientation_count['-'] > orientation_count['+'] else '+'
    for entry in cell:
        if majority_orientation == '-':
            entry[5] = flip_orientation(entry[5])
        entry.append(cell_name)
        outfile.write('\t'.join(entry) + '\n')


def process_bedfile(input_bedfile, output_bedfile, ref_chr, ref_start, ref_end):
    """
    Processes the input bed file, flips orientations based on the majority orientation in a specific region on chrX,
    appends cell names, and writes the processed entries to the output bed file.

    Args:
        input_bedfile (str): Path to the input bed file.
        output_bedfile (str): Path to the output bed file.
    """
    if input_bedfile.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'

    with open_func(input_bedfile, mode) as infile, open(output_bedfile, 'w') as outfile:
        current_cell = []
        current_cell_orientation_count = defaultdict(int)
        cell_name = ""
        for line in infile:
            if line.startswith("track"):
                if current_cell:
                    process_cell(current_cell, current_cell_orientation_count, outfile, cell_name)
                    current_cell = []
                    current_cell_orientation_count.clear()

                cell_name = line.split()[1].split('=')[1]
            else:
                fields = line.strip().split('\t')
                #import pdb; pdb.set_trace()
                chrom, start, end, _, _, strand = fields
                start = int(start)
                end = int(end)
                #if chrom == "chrX" and 140000000 <= start <= 150000000:
                if chrom == ref_chr and ref_start <= start <= ref_end:
                    current_cell_orientation_count[strand] += 1
                current_cell.append(fields)

        if current_cell:
            process_cell(current_cell, current_cell_orientation_count, outfile, cell_name)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python bedfile_orientation_processor.py <input_bedfile> <output_bedfile> <ref_chr> <ref_start> <ref_end>")
        sys.exit(1)

    input_bedfile = sys.argv[1]
    output_bedfile = sys.argv[2]
    ref_chr = sys.argv[3]
    ref_start = int(sys.argv[4])
    ref_end = int(sys.argv[5])
    process_bedfile(input_bedfile, output_bedfile, ref_chr, ref_start, ref_end)
