import os
import gzip
import sys

# Input folder containing bed files
input_folder = sys.argv[1]
output_filename = sys.argv[2]
# Create the output folder
#output_folder = os.path.join(input_folder, 'upload')

# Concatenate all bed files into a single file
with gzip.open(output_filename, 'wt') as output_file:
    for filename in os.listdir(input_folder):
        if filename.endswith('_reads.bed.gz'):
            with gzip.open(os.path.join(input_folder, filename), 'rt') as input_file:
                for line in input_file:
                    if line.startswith('track'):
                        line = line.strip()
                        name = line.split('name=')[1].split('.')[0]
                        line = f"track name={name}.bed description=\"\" {line.split(' ')[-1]}\n"
                    output_file.write(line)

print("Concatenation complete!")
print(f"Output file: {output_filename}")