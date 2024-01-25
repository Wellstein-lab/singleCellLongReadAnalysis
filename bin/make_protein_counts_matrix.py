#!/usr/bin/env python

import argparse
import sys
import os

def process_file(file_path, results):
    print(f"Processing file: {file_path}")

    with open (file_path,"r") as file:
        header = file.readline().strip().split("\t")
        # get the last element in the header
        experiment_name = header[-1]
        print(f"Experiment Name: {experiment_name}")

        for line in file:
            data = line.strip().split(' ')
#            print(f"Data Line: {data}")

            # Here, data[:4] is a list slice that extracts the first four elements (from index 0 to 3) 
            key = tuple(data[:4])
            value = int(data[4])

            if key not in results:
                results[key] = {experiment_name: value}
            else:
                results[key][experiment_name] = value
            
        
def main():
    
    parser = argparse.ArgumentParser(description="Convert input files to a matrix matched on first 4 columns matrix of 5th voclumn")
    parser.add_argument("input_directory", help="Directory containing input files to be used")
    # file output will be a csv file
    parser.add_argument("output_matrix", help="Name of the output file recommend ending with .csv")
    
    args = parser.parse_args()

    input_directory = args.input_directory
    output_matrix = args.output_matrix

    print(f"Input Directory: {input_directory}")
    print(f"Output Matrix: {output_matrix}")

    
    files = [f for f in os.listdir(input_directory) if f.endswith(".txt")]

    # Sort files for consistent order
    files.sort()
    results = {}

    for file in files:
        file_path = os.path.join(input_directory, file)
        process_file(file_path, results)

    with open(output_matrix,'w') as output_file:

        # Writing header
        header = "\t".join(["Protein", "Domain_Name", "AA_position", "Domain_sequence"] + files)
        output_file.write(header + "\n")

        # Writing data
        for key, values in results.items():
            row = ",".join(map(str, key))
            for file in files:
                row +="," + str(values.get(file,0))
            output_file.write(row + "\n")

if __name__ == "__main__":
    main()
