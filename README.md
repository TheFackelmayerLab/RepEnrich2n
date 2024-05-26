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
Identical to the original software, RepEnrich2n is run from the command line of your operating system, with arguments that provide the program with additional information. If you are not familiar with running programs from the command line, fear not! It isn’t hard if you follow our tutorial below. And certainly, there is a computer-savvy person somewhere around in your workspace who can help you to get started, or when you run into a problem. Otherwise, there are online tutorials for using the command line (e.g. using the Windows Powershell program, or the Mac Terminal app); or ask the AI of your preference, such as ChatGPT. 

Briefly, the analysis with RepEnrich2n will follow these steps:

1. Prepare your computer setup and install required software
2. Download the genome of choice, and its RepeatMasker annotation file
3. Prepare additional files for the RepEnrich2n analysis
4. Run the RepEnrich2n analysis to determine expression of repeats
5. Perform statistical evaluation of the RepEnrich2n 

Step 1 must be done only once for as many analyses you plan to run. Step 2 and a good part of step 3 must be performed once per genome. Parts of step 3 and step 4 are performed for every individual sample (e.g. separately for each biological replicate and condition). Step 5 is done once per desired statistical analysis (e.g. all replicates and conditions together). <br><br>

The documentation and detailed tutorial is under construction and will be published here soon - stay tuned!
