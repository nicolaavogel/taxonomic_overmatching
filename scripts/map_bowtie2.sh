#!/bin/bash

# Check that the correct number of arguments were provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 <fastq_list_file> <output_bam_dir> <reference_fasta>"
    exit 1
fi

FASTQ_LIST=$1
OUTPUT_DIR=$2
REFERENCE=$3

module load samtools/1.21
module load bowtie2/2.4.2

# Check if the output directory exists
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Error: Output directory '${OUTPUT_DIR}' does not exist."
    exit 1
fi

# Check if the FASTQ list file exists
if [ ! -f "$FASTQ_LIST" ]; then
    echo "Error: FASTQ list file '${FASTQ_LIST}' does not exist."
    exit 1
fi

# Check if the reference FASTA file exists
if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference FASTA file '${REFERENCE}' does not exist."
    exit 1
fi

# Run Bowtie2 index if the index does not already exist
if [ ! -f "${REFERENCE}.1.bt2l" ] || [ ! -f "${REFERENCE}.1.bt2" ]; then
    echo "Bowtie2 index not found for ${REFERENCE}, creating index..."
    bowtie2-build --threads 24 "$REFERENCE" "$REFERENCE"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create Bowtie2 index for ${REFERENCE}."
        exit 1
    fi
fi

# Read the FASTQ list file line by line
while IFS= read -r FASTQ; do
    # Check if the provided FASTQ file exists
    if [ ! -f "$FASTQ" ]; then
        echo "Error: FASTQ file '${FASTQ}' does not exist."
        continue
    fi

    OUTPUT_BAM="${OUTPUT_DIR}/$(basename "$FASTQ" .fastq)_$(basename "$REFERENCE" .fasta).bam"

    # Run Bowtie2 alignment and samtools pipeline with the specified parameters
    echo "Mapping FASTQ file '${FASTQ}' to the Bowtie2 index '${REFERENCE}' and filtering aligned reads..."
    
    bowtie2 --threads 24 -k 1000 --no-unal --mm -t -x "$REFERENCE" -U "$FASTQ" | samtools view -@ 10 -bS - | samtools sort -@ 10 -o "$OUTPUT_BAM"

    if [ $? -ne 0 ]; then
        echo "Error: Bowtie2 mapping or samtools processing failed for ${FASTQ} with reference ${REFERENCE}."
        continue
    fi

    echo "Mapping and filtering completed successfully. Output BAM file: $OUTPUT_BAM"

done < "$FASTQ_LIST"  # End of FASTQ loop

echo "All mapping tasks completed."
