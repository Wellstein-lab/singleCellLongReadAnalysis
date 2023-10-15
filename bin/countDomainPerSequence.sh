#!/bin/bash
#
#  program countDomainPerSequence.sh
#
#  purpose:
#   Given open reading frames as linearized (without carriage routine)
#   as amino acid sequences from mRNA experiments (preferabbly long read mRNA)
#   search the protein domains that are included in the measurements.
#
#  why:
#    One set of data is from single cell mRNA long read sequencing.
#    we are checking to see if there are putatively likely alternative splicing
#    isoforms in our mixture.  This would point to potentially different functions
#    and would explain differential expression of the transcript isoforms
#
#
#  
#
# Define the relative directory where your experimental files are located
#experiment_dir="/path/to/experiment_directory"
experiment_dir=$1
# Define the relative location file containing protein and domain data
#  (one pair per line, in the format "ProteinName:DomainSequence")
#  FUTURE THOUGHTS: Replace this with an API.
#

#data_file="/path/to/protein_domain_data.txt"
data_file=$2

# Check if the data file exists
#if [ ! -d "$data_file" ]; then
#    echo "Error: Data file '$data_file' not found."
#    exit 1
#fi

# Loop through each line in the data file
while IFS= read -r line; do
    # Split the line into protein name and domain sequence using ":" as a delimiter
    IFS=":" read -r protein_name domain_name aa_position domain_sequence <<< "$line"
    
    # Loop through each experiment file in the directory
    for experiment_file in "$experiment_dir"/*_linear_aa.fa; do
        # Use grep to search for the domain sequence in the experiment file
        read_count=$(grep -c "$domain_sequence" "$experiment_file")
        
        # Generate the output filename based on the input experiment file
        output_file="${experiment_file##*/}_results.txt"
        
        # Write the results to the output file
        echo "Experiment: $experiment_file" >> "$output_file"
        echo "Experiment: $experiment_file"
        echo "Protein: $protein_name" >> "$output_file"
        echo "Protein: $protein_name"
        echo "Domain_Name: $domain_name" >> "$output_file"
        echo "Domain_Name: $domain_name" 
	echo "AA_position: $aa_position" >> "$output_file"
	echo "AA_position: $aa_position"
        echo "Domain: $domain_sequence" >> "$output_file"
        echo "Domain: $domain_sequence"
        echo "Read_Count: $read_count" >> "$output_file"
	echo "read_count = $read_count"
 #       echo "---" >> "$output_file"
        echo "output_file = $output_file"
	read_count=0
    done
done < "$data_file"
