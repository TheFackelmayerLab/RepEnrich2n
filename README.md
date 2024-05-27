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
4. Run the RepEnrich2n analysis to determine read counts of every repeat
5. Perform statistical evaluation of the RepEnrich2n count data 

Step 1 must be done only once for as many analyses you plan to run. Step 2 and a good part of step 3 must be performed once per genome. Parts of step 3 and step 4 are performed for every individual sample (e.g. separately for each biological replicate and condition). Step 5 is done once per desired statistical analysis (e.g. all replicates and conditions together). <br><br>
The following documentation and detailed tutorial is under construction and is NOT YET READY.

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
The better the reference genome is chosen, the better the results of the analysis will be. For our analysis of repetitive elements in the human genome, we have used the latest available version, T2T-CHM13v2.0, instead of the standard hg38 genome. The T2T genome includes gapless telomere-to-telomere assemblies for all 22 human autosomes and chromosome X, comprising 3,054,815,472 bp of total sequence. The entire genome sequence, documentation and additional files can be found at the GitHub page of the Telomere-to-telomere consortium CHM13 project. At the same GitHub page, you also find the required RepeatMasker annotation file. These files are in the public domain and can be directly downloaded through the following links (but you don’t have to do this now, read on!):

1. [T2Tchm13v2.0 genome](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz)
2. [RepeatMasker file for T2Tchm13v2.0](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out)

For RepEnrich2n, a special format of the genome is required, called Bowtie2 indexes. So, you do not necessarily need to download the whole genome sequence, if you find the genome’s representation as Bowtie2 indexes. For the latest T2T genome, we provide these indexes for [download here](https://1drv.ms/u/c/a671e173671c80ce/ERvEhTNVkwlFlLw5oaOmdlQBwC_NORsKc2bgpgfjqvc2hg?e=jRexUJ). If you do not want to use these, feel free to create them yourself with the command (which will take a few hours): 

```bowtie2-build /path/to/file/chm13v2.0.fa name```<br>

This is also the way to go if you work with a different organism or genome version, for which you cannot find pre-made Bowtie2 indexes. Simply replace the “chm13v2.0.fa” part of the command with the filename of your downloaded raw genome sequence. The “name” part of the above command will be used as the base name for the generated indexes; you can choose it freely. Depending on the power of your computer, generating the Bowtie2 indexes will take several hours, but it must be done only once for each genome you want to analyze. The result will be six files (e.g. name.1.bt2, name.2.bt2, name.3.bt2, name.4.bt2, name.rev.1.bt2 and name.rev.2.bt2). We used T2Tchm13v2 as the base name, so the generated files were called T2Tchm13v2.1.bt2 etc. For convenience, it is advisable to copy these six files to a dedicated folder such as “bowtie2_indexes”.<br>

Regarding annotation files, the T2T consortium also supplies more specialized annotations, e.g. for a more comprehensive centromere/satellite repeat annotation, or for telomers. These can also be used for RepEnrich2n, as outlined below. However, we stick to the regular RepeatMasker annotation for the standard analysis described here.<br>

Similar resources are available for other organisms. In any case, it is important to use the soft-masked version of the genome, in which repetitive regions are shown in lower-case letters, NOT the hard-masked version, in which repetitive regions are replaced by runs of N. 

## **3. Preparations to create additional resources for the analysis**
The RepEnrich2n analysis requires a number of additional files that depend on the reference genome of your choice, and on your raw data. <br>

**First**, the analysis requires a repeat annotation file for the genome, such as the RepeatMasker file described in point 2 above. This original file contains a large number of “simple” and “low complexity” repeats, which will not be analyzed with RepEnrich2n anyway, and should be removed. For the T2T-CHM13v2.0 genome, we provide a cleaned-up version ```chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14_dropped_simple.out``` that can be used directly for [download here]( https://1drv.ms/u/c/a671e173671c80ce/EST-PelBaW9AlCLxG0pYjpcB9CgI4m12Qpl02Cqrh5R2Mw?e=sMekIM).<br>
For a different genome, search the internet for the RepeatMasker annotation file (in “native out” format) for that genome. Creating it yourself is possible, too, but the procedure is outside the scope of this tutorial. 
After downloading, clean the RepeatMasker annotation file with the Python program drop_simple.py, which we provide in the supplementary code folder (see the in-program documentation on how to use it).<br>

**Second**, RepEnrich2n requires a list of all repeats in the reference genome, and the pseudogenomes for every one of them. The original versions of RepEnrich1 and RepEnrich2 provide a setup program () that must be run before RepEnrich. This program creates the repeats list in a file (```repnames.bed```), as well as repeat pseudogenomes and their Bowtie2 representations that will be needed for RepEnrich2n (the intermediate .fa files are not needed after the bt2 files have been generated, and can be deleted). We also provide an updated version of this program, but if you want to continue with the T2T-CHM13v2.0 genome, there is no need to use this program. Instead, [download our pre-generated files from here]( https://1drv.ms/u/c/a671e173671c80ce/EX3T0z0xy1FDhIfF1Nw3bi4BIj532CPXiqHg-M--EzLcuw?e=l7UWPN). If you use any other genome (version), using the setup program is definitely necessary!<br>

No matter whether you downloaded the pre-generated files or generated them with the setup program, it is highly recommended to copy all the .bt2 files (1373 files in the case of the T2T genome) into a separate folder, and also place the ```repnames.bed``` file into the same folder. In the downloaded version, the repnames.bed file is already in the correct location (i.e. together with the .bt2 files).

**Finally**, we need to run Bowtie2 to align the raw reads of each of your sample to the genome (so, e.g., if you have two controls and two experiments, you run the alignments four times). These initial runs of Bowtie2 will create count data for reads that can be uniquely mapped to any of the repeats. This is done with the following commands: <br>

```bowtie2 -q -p 16 -x /path/to/bowtie2_indexes/T2Tchm13v2 -1 /path/to/name_1.fq -2 /path/to/ name_2.fq -S /path/to/mapped_name.sam```<br>

Here, you provide the path to the bowtie2_indexes (downloaded above), while the arguments -1 and -2 have to point to the raw read data of your end-paired RNAseq experiment, in fq format. The result of the alignment will be a .sam file, which we have to convert/translate to a different format, .bam, by the following command <br>

```samtools view -bS /path/to/mapped_name.sam > /path/to/mapped_name.bam```<br>

After this, you will have one .sam and one .bam file for each of your samples, of which the .sam file was only an intermediate and can now be deleted. 

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

```python repenrich2n.py --annotation_file /path/to/annotation/file/for/genome --alignment_bam /path/to/bam/file/of/condition --fastqfile /path/to/first/file/with/raw/reads –fastqfile2 /path/to/first/file/with/raw/reads --repeatlist /path/to/custom/repeat/list --cpus n --outpath /basepath/for/output```<br>

As the RepEnrich2n analysis is time-consuming (and can take several days of computing power for the complete set of human repeats when you use a slower computer setup), the program can also be run with a subset of repeats that you are interested in. For example, if you are interested in analysing L1 elements only, there is no need to wait many hours for the program to compute Alu or LTR elements. To limit the analysis to only those repeat of your interest, RepEnrich2n has the command line argument --repeatlist, where you can supply the path to a tab-delimited txt file with your personal choice of repeats. <br>
