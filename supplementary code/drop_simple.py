def filter_repeats_with_counts(input_file, output_file):
  """
  Filters a file to a new file, dropping lines with "Simple_repeat" or "Low_complexity" in the "repClass" column, and provides counts.

  Args:
    input_file: Path to the input file.
    output_file: Path to the output file.
  """

  total_lines = 0
  low_complexity_dropped = 0
  simple_repeat_dropped = 0
  filtered_lines = []

  with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    # Read the lines from the input file
    lines = infile.readlines()

    # Skip the header line (assuming the first line is header)
    header_line = lines[0]
    outfile.write(header_line)
    total_lines = len(lines)  # Count total lines including header

    # Filter lines based on repClass value
    for line in lines[1:]:
      if "Simple_repeat" in line:
        simple_repeat_dropped += 1
      elif "Low_complexity" in line:
        low_complexity_dropped += 1
      else:
        filtered_lines.append(line)

    # Write filtered lines to the output file
    outfile.writelines(filtered_lines)

  # Print results
  final_lines_written = len(filtered_lines)
  print(f"Original lines: {total_lines}")
  print(f"Lines dropped due to Low_complexity: {low_complexity_dropped}")
  print(f"Lines dropped due to Simple_repeat: {simple_repeat_dropped}")
  print(f"Lines in final result: {final_lines_written}")

# Example usage
input_file =  "/path/to/file/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out"
output_file = "/path/to/file/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14_dropped_simple.out"
filter_repeats_with_counts(input_file, output_file)
print(f"Filtered data from {input_file} to {output_file}")


