#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""JSON to CSV Converter

Usage:
  json_to_csv.py <input_directory> <output_csv>

Options:
  -h --help     Show this help message and exit.
"""

import os
import json
import csv
from docopt import docopt

def convert_json_to_csv(input_directory, output_csv):
    # Initialize a list to store extracted data
    data_list = []

    # Iterate through JSON files in the directory
    for filename in os.listdir(input_directory):
        if filename.endswith('.json'):
            file_path = os.path.join(input_directory, filename)
            
            # Read JSON file
            with open(file_path, 'r') as json_file:
                data = json.load(json_file)
            
            # Extract percent_mapped and file name
            percent_mapped = data.get('percent_mapped', None)
            
            # Append data to the list
            data_list.append({'file_name': filename, 'percent_mapped': percent_mapped})

    # Create a CSV file and write the data
    with open(output_csv, 'w', newline='') as csv_file:
        fieldnames = ['file_name', 'percent_mapped']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        
        # Write header
        writer.writeheader()
        
        # Write data
        writer.writerows(data_list)

    print(f"CSV file saved at: {output_csv}")

if __name__ == '__main__':
    arguments = docopt(__doc__)
    input_directory = arguments['<input_directory>']
    output_csv = arguments['<output_csv>']
    
    convert_json_to_csv(input_directory, output_csv)
