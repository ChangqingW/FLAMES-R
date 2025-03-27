# python add_bc_ub_tags.py <input_fastq> <output_fastq>
#"This script assumes that the read IDs are structured in the format @CB_UB#readID, where CB is the barcode, UB is the unique barcode (or UMI), and readID represents the specific read identifier."



import sys

def add_bc_ub_tags(input_fastq, output_fastq):
    with open(input_fastq, "r") as infile, open(output_fastq, "w") as outfile:
        while True:
            # Read four lines at a time (one record)
            lines = [infile.readline().strip() for _ in range(4)]
            if not lines[0]:  # End of file
                break
            
            # Extract barcode and umi from the read_id
            read_id = lines[0][1:]  # Remove the '@' symbol
            barcode = read_id.split('_')[0]  # First part before "_"
            umi = read_id.split('_')[1].split('#')[0]  # Part between "_" and "#"
            
            # Modify the read_id and description
            description = f"CB:Z:{barcode}\tUB:Z:{umi}"

            # Write the modified FASTQ record
            outfile.write(f"@{read_id}\t{description}\n")
            outfile.write(f"{lines[1]}\n")  # Sequence
            outfile.write(f"+\n")          # Plus line
            outfile.write(f"{lines[3]}\n") # Quality score

# Command-line arguments
if len(sys.argv) != 3:
    print("Usage: python add_bc_ub_tags.py <input_fastq> <output_fastq>")
    sys.exit(1)

input_fastq = sys.argv[1]
output_fastq = sys.argv[2]

add_bc_ub_tags(input_fastq, output_fastq)
