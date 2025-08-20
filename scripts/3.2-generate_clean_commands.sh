#!/bin/bash

# Clear any existing clean file
> clean_commands.txt

# Check if files are compressed or uncompressed
if ls ../01_merged_files/*.fastq.gz >/dev/null 2>&1; then
    echo "Files are compressed - will decompress first"
    cd ../01_merged_files
    gunzip *.fastq.gz
    cd ../03_trim
    echo "Files decompressed"
fi

# Generate commands for uncompressed files (FIXED - filename input, not stdin)
for F in ../01_merged_files/*.fastq; do
    base=$(basename "$F" .fastq)
    echo "perl /labs/Bolnick/ROL/test-samples/tag-based_RNAseq/tagseq_clipper.pl $F | cutadapt - -a AAAAAAAA -a AGATCGG -q 15 -m 25 -o ${base}.trim" >> clean_commands.txt
done

echo "Generated $(wc -l < clean_commands.txt) trimming commands"
