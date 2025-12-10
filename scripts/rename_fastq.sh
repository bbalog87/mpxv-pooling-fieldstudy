#!/bin/bash

# Usage information
usage() {
    echo "Usage: $0 -i <input_directory> [-o <output_directory>]"
    echo ""
    echo "Options:"
    echo "  -i, --input <input_directory>    Path to the folder containing raw Illumina FASTQ files."
    echo "  -o, --output <output_directory>  Path to output directory (optional)."
    echo "                                   If not specified, the renamed files will be copied to the input directory."
    echo "  -h, --help                       Show this help message."
    echo ""
    echo "Example:"
    echo "  $0 -i /path/to/input -o /path/to/output"
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check if input directory is provided
if [ -z "$INPUT_DIR" ]; then
    echo "Error: Input directory is required."
    usage
fi

# Set default output directory to input directory if not provided
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR=$INPUT_DIR
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all .fastq.gz files in the input directory
for file in "$INPUT_DIR"/*.fastq.gz; do
    # Extract the sample number (e.g., S1) and append the appropriate _R1 or _R2 to the new name
    new_name=$(basename "$file" | sed -r 's/.*_S([0-9]+)_L[0-9]+_R([12])_.*/S\1_R\2.fastq.gz/')
    
    # Copy the file to the output directory with the new name
    cp "$file" "$OUTPUT_DIR/$new_name"
done

echo "Copying completed. Files have been copied to: $OUTPUT_DIR"
