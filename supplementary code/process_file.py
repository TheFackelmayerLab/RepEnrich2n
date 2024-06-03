# a simple Python program to read a tab-delimited file and process it to
# replace commas and space characters with an underscore (_)
# and remove all double quotes (")

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Replace commas and spaces with underscores and remove double quotes
            modified_line = line.replace(',', '_').replace(' ', '_').replace('"', '')
            
            # Write the modified line to the output file
            outfile.write(modified_line)

# Run the code only when the script is executed directly (not imported as a module)
if __name__ == "__main__":
    input_file = 'input.txt'  		 # Replace with your input file name and path
    output_file = 'input_converted.txt'  # Replace with your desired output file name and path
    process_file(input_file, output_file)
