#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <FASTA_DIR> <OUTPUT_DIR>"
    exit 1
fi

# Directory containing the fasta files
FASTA_DIR="$1"

# Output directory for HSPs
OUTPUT_DIR=$(readlink -f "$2")

# Initialize counter, total_files, and time tracking
counter=0
total_files=$(ls "$FASTA_DIR"/segment_*.fasta | wc -l)
start_time=$(date +%s)

# Iterate over each fasta file in the directory
for fasta_file in $(ls "$FASTA_DIR"/segment_*.fasta | sort -V); 
do
    if [ ! -s "$fasta_file" ]; then
        echo "Skipping empty file: $fasta_file"
        continue
    fi
    # Extract the base name of the fasta file
    base_name=$(basename "$fasta_file" .fasta)
    refPath=$(readlink -f $fasta_file)
    
    #if [ -f "$OUTPUT_DIR/${base_name}.plus.segments" ]; then
    #    echo "File $OUTPUT_DIR/${base_name}.plus.segments already exists. Skipping."
    #    continue
    #fi

    OUTPUT_FOLDER=$PWD/output_$base_name/
    mkdir -p $OUTPUT_FOLDER

    if [ ! -w $OUTPUT_FOLDER ]; then
        1>&2 echo "Cannot create data directory in $OUTPUT_FOLDER because of permissions"
        (exit 5)
    fi

    DATA_FOLDER=$OUTPUT_FOLDER/data_$RANDOM/
    mkdir -p $DATA_FOLDER


    #if [ $ng -eq 0 ]; then
    #    cd $DATA_FOLDER
    #    1>&2 echo ""
    #    1>&2 echo "Converting $base_name to 2bit format"
    #    faToTwoBit $refPath ref.2bit
    #fi



    1>&2 echo ""
    1>&2 echo "Executing: \"kegalign $refPath $refPath $DATA_FOLDER --nogapped\""

    cd $OUTPUT_FOLDER
    echo "Output Directory: $OUTPUT_DIR"
    echo "Output Folder: $OUTPUT_FOLDER"
    echo "Base Name: $base_name"

    # Run the command
    kegalign $refPath $refPath $DATA_FOLDER --nogapped --debug

    # join all files matching pattern *plus* to a single file
    # Check if the expected output files exist
    if ls *plus*.segments 1> /dev/null 2>&1; then
        # Join all files matching pattern *plus* to a single file
        cat *plus*.segments > "$OUTPUT_DIR/${base_name}.plus.segments"
        # Extract start_coord from base_name
        start_coord=$(echo "$base_name" | awk -F'_' '{print $2}')
        start_coord=$((start_coord - 1))

        # Create a temporary file to store the modified segments
        temp_file="$OUTPUT_FOLDER/temp_$(basename "$fasta_file" .fasta).segments"
        touch "$temp_file"

        # Add start_coord to the specified columns in the segments file
        echo "Adding start_coord to columns in $OUTPUT_DIR/${base_name}.plus.segments"
        awk -v start_coord="$start_coord" 'BEGIN {FS=OFS="\t"} 
        {
            $2 = $2 + start_coord;
            $3 = $3 + start_coord;
            $5 = $5 + start_coord;
            $6 = $6 + start_coord;
            print
        }' "$OUTPUT_DIR/${base_name}.plus.segments" > "$temp_file"

        # Move the temporary file to the original segments file
        mv "$temp_file" "$OUTPUT_DIR/${base_name}.plus.segments"
    else
        1>&2 echo "No files matching *plus*.segments found in $OUTPUT_FOLDER"
    fi
    # Update counter and calculate elapsed time
    counter=$((counter + 1))
    current_time=$(date +%s)
    elapsed_time=$((current_time - start_time))
    avg_time_per_file=$((elapsed_time / counter))
    remaining_files=$((total_files - counter))
    estimated_time_left=$((remaining_files * avg_time_per_file))

    echo "=> Processed $counter files out of $total_files. Elapsed time: $elapsed_time seconds. Estimated time left: $estimated_time_left seconds."

    cd ..
done

# Concatenate all segments files into a single file
result_file=all_segments.plus.segments
cat "$OUTPUT_DIR"/*.plus.segments > "$OUTPUT_DIR/${result_file}"

# Remove individual segments files except the concatenated file
find "$OUTPUT_DIR" -type f -name '*.plus.segments' ! -name "${result_file}" -exec rm {} +
echo "All segments have been concatenated into $OUTPUT_DIR/${result_file} and individual segment files have been removed."

# Cleanup temporary directories
rm -rf output_*