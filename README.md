# RepEnrich2n
An extended version of RepEnrich2 for the analysis of repetitive elements in eukaryotic genomes<br>
Our work on RepEnrich2n was funded by the <br>**Hellenic Foundation for Research and Innovation (HFRI) through grant #3757, Acronym SARNAC**, to Frank O. Fackelmayer

## **Standing on the Shoulders of Giants**<br>
RepEnrich2n was built on original software provided by Nicholas Skvir (nicholas_skvir@brown.edu) at https://github.com/nerettilab/RepEnrich2 as a continuation of RepEnrich1 by Steven Criscione.<br>
The original publication to cite is: <br>
Criscione, S. W., Zhang, Y., Thompson, W., Sedivy, J. M. and Neretti, N. (2014), Transcriptional landscape of repetitive elements in normal and cancer human cells. BMC Genomics 15, 583<br>

Please make sure you have read the paper before attempting to install the program and try to run it, as it is important to understand the logic behind the analysis (especially regarding the idea of using pseudogenomes for alignments). <br>

Much of the program logic was adapted from the modernized version of RepEnrich2, provided at the repository for Galaxy tools developed in the ARTbio platform at the IBPS CNRS & UPMC, at https://github.com/ARTbio/tools-artbio/tree/main/tools/repenrich2<br>

## **Our contribution**<br>
Our extended version aims at enhancing user experience, ease of use, and improvements in the code to make it more robust, run faster, and provide improved documentation. It is written entirely in Python 3, and we run it successfully on a Apple Macbook pro with M1 Apple silicon. It should, however, run identically in any other operating system and computer setup that has functional Python 3 installed.<br>

# **Tutorial**<br>
In the same way as the original software, RepEnrich2n is run from the command line of your operating system, with arguments that provide the program with additional information. If you are not familiar with running programs from the command line, fear not! It isn’t hard if you follow our tutorial below. And certainly, there is a computer-savvy person somewhere around in your workspace who can help you to get started, or when you run into a problem. Otherwise, there are online tutorials for using the command line (e.g. using the Windows Powershell program, or the Mac Terminal app); or ask the AI of your preference, such as ChatGPT. 

Briefly, the analysis with RepEnrich2n will follow these steps:

1. Prepare your computer setup and install required software
2. Download the genome of choice, and its RepeatMasker annotation file
3. Prepare additional files for the RepEnrich2n analysis
4. Run the RepEnrich2n analysis to determine read counts for every repeat
5. Perform statistical evaluation of the RepEnrich2n count data 

Step 1 must be done only once for as many analyses you plan to run. Step 2 and a good part of step 3 must be performed once per genome. Parts of step 3 and step 4 are performed for every individual sample (e.g. separately for each biological replicate and condition). Step 5 is done once per desired statistical analysis (e.g. all replicates and conditions together). <br><br>
The following documentation and detailed tutorial is under construction and may thus NOT YET BE READY in every aspect. Please let us know if you encounter problems or need clarifications.

## **1. Preparing your computer setup**
Importantly, the RepEnrich2n program requires additional software to be installed. First, as the analysis is done with Python 3, it is important that your computer has a functional Python 3 installation. The deprecated older Python 2.7 is not compatible with Python 3 code, and cannot be used. This won’t be a problem if you install Python now, but if your computer is older and has pre-installed Python, make sure the version is 3.8 or higher. To learn more about Python 3, see the official website https://www.python.org/ from where you can also download installers for almost every operating system under the sun. We have coded RepEnrich2n in Python on a Mac, but it should run on all other systems, too, maybe with minor adaptations.<br>

We also need other software: 
1. [BioPython](https://biopython.org/). Install this with the pip command of Python, ```pip install BioPython```
2. [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2). Install this from [here](https://bowtie-bio.sourceforge.net/bowtie2)
3. [samtools](https://www.htslib.org/) 
4. [bedtools](https://bedtools.readthedocs.io/)

Samtools and Bedtools are a little bit harder to install. Both for Windows and Mac, you can install it with Conda (Anaconda or Miniconda), after first installing either one from the official websites at https://www.anaconda.com/products/distribution#download-section 
or https://docs.conda.io/en/latest/miniconda.html. If you want to also use Python 3 for other projects, consider creating a separate virtual environment for RepEnrich2n; see the official Python documentation or ask your favorite AI to learn about Python environments. 
<br>Finally, Samtools and Bedtools can be installed with the commands <br>

```conda install -c bioconda samtools``` and ```conda install -c bioconda bedtools```

For Mac, both programs can alternatively also be installed with [Homebrew](https://brew.sh/), if that software package management system is what you are more familiar with.<br>

**For reference: We currently use Python 3 (version 3.11.4), Biopython (version 1.83), Bowtie2 (version 2.4.3), samtools (version 1.19.2), and bedtools (v2.31.1)**

## **2. Selecting the suitable reference genome**
*The better the reference genome is chosen, the better the results of the analysis will be*. For our analysis of repetitive elements in the human genome, we have used the latest available version, T2T-CHM13v2.0, instead of the standard hg38 genome. The T2T genome includes gapless telomere-to-telomere assemblies for all 22 human autosomes and chromosome X, comprising 3,054,815,472 bp of total sequence. The entire genome sequence, documentation and additional files can be found at the GitHub page of the [Telomere-to-telomere consortium CHM13 project](https://github.com/marbl/CHM13). At the same GitHub page, you also find the required RepeatMasker annotation file. These files are in the public domain and can be directly downloaded through the following links (but you don’t have to do this now, read on!):

1. [T2Tchm13v2.0 genome](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz)
2. [RepeatMasker file for T2Tchm13v2.0](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out)

For RepEnrich2n, a special format of the genome is required, called Bowtie2 indexes. So, you do not necessarily need to download the whole genome sequence, if you find the genome’s representation as Bowtie2 indexes. For the latest T2T genome, we provide these indexes for [download here](https://1drv.ms/u/c/a671e173671c80ce/ERvEhTNVkwlFlLw5oaOmdlQBwC_NORsKc2bgpgfjqvc2hg?e=jRexUJ). If you do not want to use these, feel free to create them yourself with the command (which will take a few hours): 

```bowtie2-build /path/to/file/chm13v2.0.fa name```<br>

This is also the way to go if you work with a different organism or genome version, for which you cannot find pre-made Bowtie2 indexes. Simply replace the “chm13v2.0.fa” part of the command with the filename of your downloaded raw genome sequence. The “name” part of the above command will be used as the base name for the generated indexes; you can choose it freely. Depending on the power of your computer, generating the Bowtie2 indexes will take several hours, but it must be done only once for each genome you want to analyze. The result will be six files (e.g. name.1.bt2, name.2.bt2, name.3.bt2, name.4.bt2, name.rev.1.bt2 and name.rev.2.bt2). We used T2Tchm13v2 as the base name, so the generated files were called T2Tchm13v2.1.bt2 etc. For convenience, it is advisable to copy these six files to a dedicated folder such as “bowtie2_indexes”.<br>

Regarding annotation files, the T2T consortium also supplies more specialized annotations, e.g. for a more comprehensive centromere/satellite repeat annotation, or for telomers. These can also be used for RepEnrich2n, as outlined below. However, we stick to the regular RepeatMasker annotation for the standard analysis described here.<br>

Similar resources are available for other organisms. In any case, it is important to use the soft-masked version of the genome, in which repetitive regions are shown in lower-case letters, NOT the hard-masked version, in which repetitive regions are replaced by runs of N. 

## **3. Preparations to create additional resources for the analysis**
The RepEnrich2n analysis requires a number of additional files that depend on the reference genome of your choice, and on your raw data. <br>

**First**, the analysis requires a repeat annotation file for the genome, such as the RepeatMasker file described in point 2 above. This original file contains a large number of “simple” and “low complexity” repeats, which will not be analyzed with RepEnrich2n anyway, and should be removed. For the T2T-CHM13v2.0 genome, we provide a cleaned-up version ```chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14_dropped_simple.out``` for [download here]( https://1drv.ms/u/c/a671e173671c80ce/EST-PelBaW9AlCLxG0pYjpcB9CgI4m12Qpl02Cqrh5R2Mw?e=sMekIM), which can be used directly.<br>
For a different genome, search the internet for the RepeatMasker annotation file (in “native out” format) for that genome. Creating it yourself is possible, too, but the procedure is outside the scope of this tutorial. 
After downloading, clean the RepeatMasker annotation file with the Python program drop_simple.py, which we provide in the supplementary code folder (see the in-program documentation on how to use it).<br>

**Using different annotation files:** Here, we use the typical annotation file for the T2T genome, which is the native output of a RepeatMasker analysis. But that is not the only annotation file that can be used, and other repeat annotations are available for the more prominent genomes, too. Unfortunately, these annotation files were prepared by different consortia using different programs, and thus do not have a unified format. So, the first step to utilize these alternative annotation files is to re-format the content so that it can be used by ```RepEnrich2n_setup.py```. This is most easily done by a short program in Python, which extracts the required information from the original file and creates a new file with the format required by our programs. Basically, all information should be in the original file, but in different columns, which simply must be re-ordered. So, the code would read the annotation file line-by-line, split the content of the line into separate junks of information depending on the delimiter (usually, but not always, a tab character ´\t´), and write a new file with these junks of information in a different order. An (arbitrary) example for code which would keep the first three column, then use original column 10 as new column 4, create two empty columns, and finally use original column 4 as new columns 7 and 8 would be: 
```
def reorder_columns(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # change this for a different delimiter
            columns = line.strip().split('\t')	
            # Extract the required columns
            new_columns = columns[:3] + [columns[9], '', '', columns[3], columns[3]]
            # Join the new columns with tab delimiter and write to the output file
            outfile.write('\t'.join(new_columns) + '\n')

# run the code only when the script is executed directly (not imported as a module)
if __name__ == "__main__": 
input_file = 'input.txt' # Replace with your input file name and path
output_file = 'output.txt' # Replace with your desired output file name and path reorder_columns(input_file, output_file)
```
<br>

To adapt this sample code to your specific case, the first step will be to open the alternative annotation file and check in which column the required information is stored. Then, change the code to re-order the columns as described below (it is possible that additional code is needed e.g. to re-format text). Run it on your annotation file, and check that the order of the columns in the generated output file is identical to the example file we provide for download, ```chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14_dropped_simple.out```. Most of the columns in this file will not be used by RepEnrich2, and can be filled with a zero (0). In the end, the file must have the respective chromosome in column 5, the start and end of the sequence in columns 6 and 7, the name of the repeat type in column 10, and its class and family in column 11 (separated by a slash “/”). It will look something like this:<br>
0 | 0 | 0 | 0 | chr1 | 116796047 | 121405145 | 0 | 0 | name | class/family
<br>
We plan to soon provide a ready-to-use version of the annotation file for a more comprehensive centromere/satellite repeat annotation, adapted from the one published by the UCSC Genomics Institute [here](https://genome.ucsc.edu/cgi-bin/hgTables?db=hub_3671779_hs1&hgta_group=map&hgta_track=hub_3671779_censat&hgta_table=hub_3671779_censat&hgta_doSchema=describe+table+schema). Others may also follow later - stay tuned!

**Second**, RepEnrich2n requires a list of all repeats in the reference genome, and the pseudogenomes for every one of them. The original versions of RepEnrich1 and RepEnrich2 use a setup program that must be run before RepEnrich. This program reads information from the annotation file to create the repeats list in a file (```repnames.bed```), as well as repeat pseudogenomes and their Bowtie2 representations that will be needed for RepEnrich2n (the intermediate .fa files are not needed after the bt2 files have been generated, and can be deleted). We also provide an updated version of this program, RepEnrich2n_setup.py in the supplementary code folder, but if you want to continue with the T2T-CHM13v2.0 genome, there is no need to use this program. Instead, [download our pre-generated files from here]( https://1drv.ms/u/c/a671e173671c80ce/EX3T0z0xy1FDhIfF1Nw3bi4BIj532CPXiqHg-M--EzLcuw?e=l7UWPN). <br>
If, however, you want to analyze repeats in any other genome (version), or utilize an alternative annotation file, running our setup program is definitely necessary!<br>

No matter whether you downloaded the pre-generated files or generated them with the setup program, it is highly recommended to copy all the .bt2 files (1373 files in the case of the T2T genome) into a separate folder (e.g. called “bt_files”), and also place the ```repnames.bed``` file into the same folder. In the downloaded version, the repnames.bed file is already in the correct location (i.e. together with the .bt2 files).

**Finally**, we need to run Bowtie2 to align the raw reads of each of your sample to the genome (so, e.g., if you have two controls and two experiments, you run the alignments four times). This is done with the following commands: <br>

```bowtie2 -q -p 16 -x /path/to/bowtie2_indexes/T2Tchm13v2 -1 /path/to/name_1.fq -2 /path/to/name_2.fq -S /path/to/mapped_name.sam```<br>

Here, you provide the path to the bowtie2_indexes (downloaded above), while the arguments -1 and -2 have to point to the raw read data of your end-paired RNAseq experiment, in fq format. 
Choose “name” for each of your sample, e.g. “ctrl_1” or “treated_1”. The result of the alignment will be a .sam file, which we have to convert/translate to a different format, .bam, by the following command <br>

```samtools view -bS /path/to/mapped_name.sam > /path/to/mapped_name.bam```<br>

This will generate a .bam file for your data. After this you will have (in our example), the four files “ctrl_1.bam”, “ctrl_2.bam”, “treated_1.bam” and “treated_2.bam”. The intermediary .sam files are now no longer needed and can be deleted.

## **4. Running the RepEnrich2n program**
The RepEnrich2n program is run from the computer’s command line, using the Python 3 interpreter. The program requires additional information passed through arguments, as follows: <br>

--annotation_file [path to annotation file for the chosen genome] <br>
--alignment_bam [path to the .bam file of the condition] <br>
--fastqfile [path to the first file with raw data from RNAseq] <br>
--fastqfile2 [path to the second file with raw data from RNAseq, for paired-end sequencing] <br>
--repeatlist [path to the custom repeat list if wanted] <br>
--cpus [number of cpu cores to run the analysis on; we use 8, but default is 1] <br>
--outpath [basepath for the output] <br>

The command to run RepEnrich2n is simply ```python repenrich2n.py```, followed by the arguments and their values, for example: 

```python repenrich2n.py --annotation_file /path/to/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14_dropped_simple.out --alignment_bam /path/to/mapped_name.bam --fastqfile /path/to/name_1.fq --fastqfile2 /path/to/name_2.fq --repeatlist /path/to/custom/repeat_list.txt --cpus n --outpath /basepath/for/output```<br>

Of these, the arguments repeatlist, cpus and outpath are not strictly required, but are implemented to provide extra convenience. Also, if your RNAseq raw data are not end-paired (which is only rarely the case), the fastqfile2 is not required and should not be declared.<br>

**Watch out:** In its current version, the program assumes that it is run from in the folder/directory that contains the .bt2 files and the repnames.bed file. So, move to this folder/directory with ```cd /path/to/bt_files```<br> 

Note that running RepEnrich2n requires a lot of processing power, and - even on a reasonably fast MacBook - can take 50 hours or more. As you have to do this separately for each of your experimental conditions, we are talking a week or more for a full analysis. To overcome this issue, we have implemented a function to allow the program to be run with a subset of repeats that you are interested in. For example, if you are interested in analyzing L1 elements only, there is no need to wait many hours for the program to compute Alu or LTR elements. To limit the analysis to only the repeats of your interest, RepEnrich2n has the command line argument --repeatlist, where you can supply the path to a tab-delimited txt file with your personal choice of repeats. <br>
For convenience, we provide ready-to-use lists for several abundant classes and families of repeats in the human T2T genome; these lists can be used directly as arguments for --repeatlist, e.g. ```--repeatlist LINE_class.txt```, or serve as examples on how to create your own list. If you are interested in ALL repeats, you do not have to declare the --repeatlist argument, or you can pass the value allrepeats.txt, as ```--repeatlist allrepeats.txt```; both options are equivalent, but explicitly passing allrepeats.txt is significantly faster. 

We know that it is very inconvenient to dedicate a computer full-time to the analysis, and (especially for a notebook) not be able to shut it down and move it somewhere else. Our RepEnrich2n therefore allows you to stop the execution at any time, e.g. by pressing ctrl-C in the shell. This may give you an ugly error message, but the program will understand and remember where it was stopped, and will continue with the analysis at this point after being re-started. The re-startability of the program also allows you to kill/interrupt single bowtie2 processes in case they get stuck - which happens quite regularly in our hands. When the program has finished, it will tell whether the analysis was completed successfully, or must be re-run to process the (hopefully few) repeats that were missing due to interrupting a stuck bowtie2 process. <br> 

After successful completion, the result of every RepEnrich2n analysis is saved in four files, which contain the count data in tab-delimited txt file format as follows: <br>

1. **unique_mapper_counts.tsv**: uniquely mapped counts per repeat type
2. **fraction_counts.tsv**: fractionally mapped counts per repeat type
3. **family_fraction_counts.tsv**: fractionally mapped counts per repeat family
4. **class_fraction_counts.tsv** : fractionally mapped counts per repeat class

The most important of these is the fraction_counts.tsv file, which contain the final results as read counts per individual repeat type (e.g. L1PA5); these counts are the sum of the uniquely mapped reads in unique_mapper_counts.tsv, and the fractional counts for the same repeat type. The last two files contain count data aggregated for repeat family (e.g. L1) and repeat class (e.g. LINE). These allow statistical evaluation (see below) on the level of whole families and classes. 


## **5. Statistical Evaluation and Interpretation of the RepEnrich2n results**
The last step of the analysis is to statistically evaluate the counts data generated by RepEnrich2n. There are many ways to do this, and it is beyond the scope of this tutorial to give detailed protocols. Most likely, if you are doing RNAseq, you already have your own analysis pipeline. Most conveniently, the analysis can be performed with a script in R (a powerful open-source programming language and environment for statistical computing and graphics (https://www.r-project.org/). R is available for every operating system, and easy to install. Using R can be quite a challenge, though, and it is highly recommended to follow a basic tutorial (e.g. on Youtube). For a more in-depth introduction, see [this official site]( https://cran.r-project.org/manuals.html).
<br>

For the purpose of this tutorial, we provide a sample script (R_samplescript.txt) for data from an RNA-seq experiment with two conditions (e.g. “control” and “treated”), each in duplicate. For different data, the R script can be adapted accordingly. Before the script can be run, the count data from the (in our case four) different RepEnrich2n runs must be consolidated into one tab-delimited file, e.g. using Excel. The sample script assumes that the data is a tab-delimited .txt file with five columns, where the first one has the individual repeat name, the second and third column have the count data from “control” samples, and the fourth and fifth column have the data from “treated” samples. The first row must provide the titles “repeat_name”, and the names of the samples, e.g. “ctrl1”, “ctrl2” etc. in this way:<br>

|repeat_name |   ctrl1   |   ctrl2     | treated1  | treated2   |
|-------------------|-------------|---------------|---------------|-----------------|
|5S                     | 1140       | 1891        | 1122        | 1289            |
|7SK                   | 820304  | 197764    | 795431   | 685558       |
|7SLRNA           | 7145018 | 1585139 | 6442318 | 6248160      |
|ACRO1             | 64            |82               | 53            | 29                |
.<br>
.<br>
.<br>

This file can be created by copying the count data of any of the results files of RepEnrich2n. The most comprehensive “full” analysis is that of data from the fraction_counts.tsv files, but you can also analyze the higher-order clusters of families or classes, using data from the family_fraction_counts.tsv or class_fraction_counts.tsv files, respectively. <br>
Our R_samplescript.txt will read the data from this consolidated file, but in order to find it, the path to the file must be edited in the fourth line of the code ```counts_data <- read.delim("/path/to/the/file/RepEnrichResults.txt", row.names="repeat_name")```, to point to your file.
The script then uses [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) for statistical evaluation, and [ggplot2](https://ggplot2.tidyverse.org/) for visualization of the final result. Both packages must be installed in R, as per the instructions on their respective website, before the script is run. <br>
In the analysis, we set a significance threshold to 0.05 after correction for multiple hypothesis test, print the 25 most significant hits to a file “results.txt”, and visualized the differential expression of repetitive elements with a volcano plot. We have added comments for documentation into the R script, to make it easy to adapt it to other demands, e.g. to save more that 25 hits, or change the appearance of the volcano plot. <br>
As said above, this is only a sample script, and much much more can be done in R, depending on the scientific question. It is definitely worth investing some time into learning the basics of R!
