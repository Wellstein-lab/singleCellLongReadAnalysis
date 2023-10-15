#!/usr/bin/env python

import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Convert structured data to a wide matrix format.")
    parser.add_argument("input_file", help="Path to the input file")
    parser.add_argument("output_file", help="Path to the output file")
    parser.add_argument("-o", "--order", nargs='+', default=["Experiment", "Protein", "Domain_Name", "Domain", "AA_position", "Read_Count"],
                        help="Desired order of keys separated by spaces")
    
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    order = args.order

    # Create a dictionary to store the key-value pairs
    data = {}

    # Initialize the data dictionary with empty values for the desired keys
    for key in order:
        data[key] = []

    # Process the input file
    try:
        with open(input_file, "r") as infile:
            for line in infile:
                parts = line.strip().split(": ")
                if len(parts) == 2:
                    key, value = parts
                    if key in data:
                        data[key].append(value)
    except FileNotFoundError:
        print(f"Error: The input file '{input_file}' was not found.")
        sys.exit(1)

    # Write the data to the output file
    with open(output_file, "w") as outfile:
        for key in order:
            values = ", ".join(data[key])
            outfile.write(f"{key}, {values}\n")

if __name__ == "__main__":
    main()
