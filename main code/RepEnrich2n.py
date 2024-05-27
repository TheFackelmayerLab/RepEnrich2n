# RepEnrich2n, an extended version of RepEnrich2
# (c)2024 Frank O. Fackelmayer
# published under MIT license

import json
import argparse
import csv
import os
import shlex
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

#
# read the command line arguments and perform checks to increase robustness
# 

parser = argparse.ArgumentParser(description='''
             Repenrich aligns reads to Repeat Elements pseudogenomes\
             and counts aligned reads. RepEnrich_setup must be run\
             before its use''')
parser.add_argument('--annotation_file', action='store',
                    metavar='annotation_file',
                    help='RepeatMasker.org annotation file for your\
                          organism. The file may be downloaded from\
                          RepeatMasker.org. E.g. hg19_repeatmasker.txt')
parser.add_argument('--alignment_bam', action='store', metavar='alignment_bam',
                    help='Bam alignments of unique mapper reads.')
parser.add_argument('--fastqfile', action='store', metavar='fastqfile',
                    help='File of fastq reads mapping to multiple\
                          locations. Example: /data/multimap.fastq')
parser.add_argument('--fastqfile2', action='store', dest='fastqfile2',
                    metavar='fastqfile2', default='',
                    help='fastqfile #2 when using paired-end option.\
                          Default none')
parser.add_argument('--cpus', action='store', dest='cpus', metavar='cpus',
                    default="1", type=int,
                    help='Number of CPUs. The more cpus the\
                          faster RepEnrich performs. Default: "1"')
parser.add_argument('--repeatlist', action='store', dest='repeatfile', type=str, 
                    default='',
                    help="List of repeats to be quantified, instead of taking them from --alignment.bam")
parser.add_argument('--outpath', action='store', dest='outpath', type=str, 
                    default='',
                    help="Base path for all output - default is current path")

args = parser.parse_args()


# parameters from arguments
repeat_names_file = None			# empty placeholder
annotation_file = args.annotation_file
unique_mapper_bam = args.alignment_bam
fastqfile_1 = args.fastqfile
fastqfile_2 = args.fastqfile2
cpus = args.cpus
repeat_names_file = args.repeatfile
outpath = args.outpath


# check that the required command line arguments have been provided
# not to waste time on bedtools analysis if the other files for bowtie do not exist 
# do not exist e.g. because their name or folder has a typo
def check_file_existence(args):
  missing_files = []
  required_files = [args.annotation_file, args.alignment_bam, args.fastqfile]

  for filename in required_files:
    if not os.path.isfile(filename):
      missing_files.append(filename)

  # Check fastqfile2 only if it's provided in the arguments
  if args.fastqfile2 and not os.path.isfile(args.fastqfile2):
    missing_files.append(args.fastqfile2)

  if missing_files:
    error_message = "Error: required file(s) missing: \n" + "\n".join(missing_files)
    print(error_message)
    raise SystemExit(1)  # Exit program with error code 1

# call the function to check files, continue only if all files exist
check_file_existence(args)


# check that the programs we need are available
try:
    subprocess.call(shlex.split("coverageBed -h"),
                    stdout=open(os.devnull, 'wb'),
                    stderr=open(os.devnull, 'wb'))
    subprocess.call(shlex.split("bowtie2 --version"),
                    stdout=open(os.devnull, 'wb'),
                    stderr=open(os.devnull, 'wb'))
except OSError:
    print("Error: Bowtie2 or bedtools not available")
    raise


# determine whether data are paired-end or not
if args.fastqfile2:
    paired_end = True
else:
    paired_end = False

#
# all startup checks are completed at this point
#

print("__________________________________________________________________\n")

#
# define helper functions 
#

# function to check if the first element in a list is a numerical value
def starts_with_numerical(list):
    try:
        if len(list) == 0:
            return False
        int(list[0])
        return True
    except ValueError:
        return False

# define a text importer for .out/.txt format of repbase
def import_text(filename, separator):
    csv.field_size_limit(sys.maxsize)
    file = csv.reader(open(filename), delimiter=separator,
                      skipinitialspace=True)
    return [line for line in file if starts_with_numerical(line)]

# define an importer for a file with a user-defined list of repeats 
# (unless they are generated from the --annotation_file)
def read_repeat_names(file_path):
    with open(file_path, 'r') as file:
        repeat_list = [line.strip() for line in file]
    return sorted(list(set(repeat_list)))

# set a reference repeat list for the script
# the list is either passed in --repeatlist
# or generated from the repnames.bed file 
# change the next line if the renames.bed file is in a different location 
# than the current folder/directory
input_file = 'repnames.bed'	# add path if necessary
if repeat_names_file:
    output_file = os.path.join(outpath, 'filtered_repeats.bed')
    print(f"Using your list of repeats in {repeat_names_file}")
    repeat_list = read_repeat_names(repeat_names_file)
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            columns = line.strip().split()
            if len(columns) >= 4 and columns[3] in repeat_list:
                outfile.write(line)
    use_these = output_file		# use the shortened list
else:
    print("Generating list of repeats, please wait...")
    repeat_list = [listline[9].translate(
       str.maketrans('()/', '___')) for listline in import_text(annotation_file, '\t')]
    repeat_list = sorted(list(set(repeat_list)))
    use_these = input_file		# use the original list


# we will save bedtools analysis in a cache file so that bedrolls must run only once
# and can later retrieve results from the cached file in all subsequent runs with the same data
#
# this allows the program to be stopped later (in the time-consuming step of
# determining fractional counts) and be re-started.
# NOTE: we will delete the cache file when the analysis has completed successfully 

cache_file = os.path.join(outpath, 'repeat_counts_cache.json')

def generate_data():
    # unique mapper counting
    cmd = f"bedtools bamtobed -i {unique_mapper_bam} | bedtools coverage -b stdin -a {use_these}"
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    bedtools_counts = p.communicate()[0].decode().rstrip('\r\n').split('\n')

    # parse bedtools output
    counts = defaultdict(int)  # key: repeat names, value: unique mapper counts
    sumofrepeatreads = 0
    for line in bedtools_counts:
        line = line.split('\t')
        counts[line[3]] += int(line[4])
        sumofrepeatreads += int(line[4])

    # save unique mapper counts to file
    output_tsv = os.path.join(outpath, 'unique_mapper_counts.tsv')
    with open(output_tsv, 'w') as fout:
        fout.write("#element\tcount\n")
        for count in sorted(counts):
            fout.write(f"{count}\t{counts[count]}\n")
    
    return dict(counts)

# Define `counts` to ensure it's always available
counts = defaultdict(int)

try:
    # Try to load data from the cache file
    with open(cache_file, 'r') as f:
        data = json.load(f)
        print("Loaded bedtools results from cache file:", cache_file)
        sumofrepeatreads = sum(data.values())
        # Convert data back to defaultdict(int)
        for key, value in data.items():
            counts[key] = value
except FileNotFoundError:
    # If cache file is not found, generate data and save it
    print("Generating data with bedtools, please wait...")
    data = generate_data()
    with open(cache_file, 'w') as f:
        json.dump(data, f)
#        print("Saved data to cache file:", cache_file)
    sumofrepeatreads = sum(data.values())
    # Use the generated data directly
    counts = defaultdict(int, data)


# print the result of uniquely mapped repeats
print(f"Identified {sumofrepeatreads} unique reads that mapped to repeats.") 


# counting uniquely mapped reads is finished here


# go on with determining fractional counts for every repeat  
# ----------------------------------------------------------
# modified to handle potential bowtie errors gracefully
def run_bowtie(args):
  """
  Writes to files to save memory.
  """
  metagenome = args
  b_opt = "-k 1 -p 2 --quiet --no-hd --no-unal"
  if paired_end is True:
    command = shlex.split(f"bowtie2 {b_opt} -x {metagenome}"
                          f" -1 {fastqfile_1} -2 {fastqfile_2}")
  else:
    command = shlex.split(f"bowtie2 {b_opt} -x {metagenome} {fastqfile_1}")

# Print metagenome and constructed command for debugging
#  print(f"Processing metagenome: {metagenome}")
#  print(f"Command: {' '.join(command)}")

  # Capture both stdout and stderr
  result = subprocess.run(command, check=True, capture_output=True, text=True)
  bowtie_align = result.stdout.rstrip('\r\n').split('\n')

  # Check for errors in stderr
  if result.stderr:
    print(f"Error running bowtie2 for {metagenome}:\n{result.stderr}")
    return  # Indicate error and return

  # Open in write mode to create the file if missing, instead of append mode
  with open(f"{metagenome}.reads", "w") as readfile:
    print(f"processed: {metagenome}")
    for line in bowtie_align:
      read = line.split()[0].split("/")[0]
      if read:
        readfile.write(f"{read}\n")


# function to create a new args_list that ONLY contain the entries of repeat_list 
# IF NOT the {metagenome}.reads file does not exist
# This handles second (and later) runs when the first run failed to produce a .reads file 
# for one or more repeats, because bowtie failed or was interrupted.
# This also allows the program to be stopped at any point during the analysis (e.g. by pressing
# ctrl-c), and continue later without having to start from the beginning again.

def check_and_add_to_args(metagenome):
  """Checks if the .reads file exists and adds metagenome to args_list if it does NOT exist.
  Args:
      metagenome: The name of the metagenome.

  Returns:
      The metagenome if the corresponding .reads file does NOT exist, None otherwise.
  """

  # Construct the file path
  file_path = f"{metagenome}.reads"

  # Check if the file exists using os.path.exists() for efficiency
  if not os.path.exists(file_path):
#    print(f"progress: '{file_path}' missing, creating it now.")   # can be un-commented for second or later runs
    return metagenome
  else:
    # Skip entries with existing files
    # print(f"File '{file_path}' exists, skipping {metagenome}.")
    return None  # Indicates skipped entry

# Filter repeat_list to include only entries with missing .reads files
# this skips analysis of already completed repeats, which is important when re-running the program
args_list = [metagenome for metagenome in repeat_list if check_and_add_to_args(metagenome) is not None]

# multimapper parsing with ThreadPoolExecutor
# Use ThreadPoolExecutor with the filtered args_list (repeats with missing .reads files)
# because it works for multiprocessing also on macOS with Python >3.8
with ThreadPoolExecutor(max_workers=cpus) as executor:
  results = executor.map(run_bowtie, args_list)


# Aggregate results (avoiding race conditions)
metagenome_reads = defaultdict(list)  # metagenome: list of multimap reads

# Now we read .reads files to populate metagenomes_reads
# Iterate through metagenomes and gracefully handle those that do not exist
is_complete = "ANALYSIS COMPLETED SUCCESSFULLY"
for metagenome in repeat_list:
  # Construct the file path
  file_path = f"{metagenome}.reads"

  # Check if the file exists using os.path.exists() for efficiency
  if os.path.exists(file_path):
    # Open the file in read mode
    with open(file_path, "r") as readfile:
      for read in readfile:
            metagenome_reads[metagenome].append(read)
      # read are only once in list
      metagenome_reads[metagenome] = list(set(metagenome_reads[metagenome]))
  else:
     print(f"File '{file_path}' not found, skipping {metagenome}.")
     is_complete = " >>>> WARNING: ANALYSIS IS INCOMPLETE, RE-RUN PROGRAM TO COMPLETE" 

# Now 'metagenome_reads' dictionary only contains entries for processed metagenomes 
# with existing .reads file, and a flag for complete analysis is set appropriately 


# implement repeats_by_reads from the inverse dictionary metagenome_reads
repeats_by_reads = defaultdict(list)  # readids: list of repeats names
for repname in metagenome_reads:
    for read in metagenome_reads[repname]:
        repeats_by_reads[read].append(repname)
for repname in repeats_by_reads:
    repeats_by_reads[repname] = list(set(repeats_by_reads[repname]))
    # this repeats_by_reads dictionary is far too big

# 3 dictionaries and 1 pointer variable to be populated
fractionalcounts = defaultdict(float)
familyfractionalcounts = defaultdict(float)
classfractionalcounts = defaultdict(float)
sumofrepeatreads = 0

# Update counts dictionary with sets of repeats (was "subfamilies")
# matched by multimappers
for repeat_set in repeats_by_reads.values():
    repeat_set_string = ','.join(repeat_set)
    counts[repeat_set_string] += 1
    sumofrepeatreads += 1

print(f'Identified {sumofrepeatreads} multimapper repeat reads')

# Populate fractionalcounts
for key, count in counts.items():
    key_list = key.split(',')
    for i in key_list:
        fractionalcounts[i] += count / len(key_list)

# build repeat_ref for easy access to rep class and rep families
repeat_ref = defaultdict(dict)
repeats = import_text(annotation_file, '\t')
for repeat in repeats:
    repeat_name = repeat[9].translate(str.maketrans('()/', '___'))
    try:
        repclass = repeat[10].split('/')[0]
        repfamily = repeat[10].split('/')[1]
    except IndexError:
        repclass, repfamily = repeat[10], repeat[10]
    repeat_ref[repeat_name]['class'] = repclass
    repeat_ref[repeat_name]['family'] = repfamily

# Populate classfractionalcounts and familyfractionalcounts
for key, value in fractionalcounts.items():
    classfractionalcounts[repeat_ref[key]['class']] += value
    familyfractionalcounts[repeat_ref[key]['family']] += value


# we're done with the analysis here

# now we only need to write the final results to files
# first, assemble file names for output at outpath
class_output  = os.path.join(outpath, 'class_fraction_counts.tsv')
family_output = os.path.join(outpath, 'family_fraction_counts.tsv')
repeat_output = os.path.join(outpath, 'fraction_counts.tsv')

# second, print class-, family- and fraction-repeats counts to files
with open(class_output, 'w') as fout:
    for key in sorted(classfractionalcounts):
        fout.write(f"{key}\t{round(classfractionalcounts[key], 2)}\n")

with open(family_output, 'w') as fout:
    for key in sorted(familyfractionalcounts):
        fout.write(f"{key}\t{round(familyfractionalcounts[key], 2)}\n")

with open(repeat_output, 'w') as fout:
    for key in sorted(fractionalcounts):
        fout.write(f"{key}\t{repeat_ref[key]['class']}\t"
                   f"{repeat_ref[key]['family']}\t"
                   f"{round(fractionalcounts[key], 2)}\n")

# report whether analysis was completed successfully or program should be re-run 
print(f"\n{is_complete}\n")

# final clean-up, deleting the bedtools cache file (only after successfully completing the analysis)
if is_complete == "ANALYSIS COMPLETED SUCCESSFULLY":
    try:
        os.remove(cache_file)
    except:
        pass