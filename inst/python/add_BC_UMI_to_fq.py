# Usage: python add_BC_UMI_to_fq.py <input_fastq[.gz]> <output_fastq[.gz]>
# This script assumes that the read IDs are structured in the format
# @CB_UB#readID, where CB is the barcode, UB is the unique barcode (or UMI),
# and readID represents the specific read identifier. The barcode and UMI are
# written into the read header as SAM tags (CB:Z:/UB:Z:) in the comment field,
# so tools such as minimap2 -y can propagate them to the resulting BAM.
#
# Equivalent awk one-liner:
#   zcat in.fq.gz | awk 'NR%4==1{split(substr($1,2),a,"#"); split(a[1],b,"_"); \
#     print $1"\tCB:Z:"b[1]"\tUB:Z:"b[2]; next}1' | gzip > out.fq.gz
# Using mawk and pigz  (pigz -dc / pigz) could be faster than the python script

import sys
import gzip


def _open(path, mode):
    # Transparently handle gzip-compressed FASTQ based on the file extension.
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def add_bc_ub_tags(input_fastq, output_fastq):
    with _open(input_fastq, "rt") as infile, _open(output_fastq, "wt") as outfile:
        for n, header in enumerate(infile):
            # Read the remaining three lines of the record.
            try:
                seq, _plus, qual = next(infile), next(infile), next(infile)
            except StopIteration:
                sys.exit(f"Truncated FASTQ record starting at line {n + 1}")

            # Drop any existing comment, keep only the read ID itself.
            read_id = header[1:].rstrip("\n").split(None, 1)[0]

            # Extract barcode and umi from the read_id (@CB_UB#readID).
            bc_ub, sep, _rest = read_id.partition("#")
            barcode, bc_sep, umi = bc_ub.partition("_")
            if not sep or not bc_sep or not umi:
                sys.exit(
                    f"Read {read_id!r} (record {n // 4}) is not in the "
                    "expected @CB_UB#readID format"
                )

            # Write the modified FASTQ record; quality line is preserved verbatim.
            outfile.write(f"@{read_id}\tCB:Z:{barcode}\tUB:Z:{umi}\n{seq}+\n{qual}")


# Command-line arguments
if len(sys.argv) != 3:
    sys.exit("Usage: python add_BC_UMI_to_fq.py <input_fastq[.gz]> <output_fastq[.gz]>")

add_bc_ub_tags(sys.argv[1], sys.argv[2])
