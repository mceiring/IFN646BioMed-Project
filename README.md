# IFN646BioMed-Project
Project 3 - Reproducibility and robustness of genomics results

Data
For this project, you are in charge of generating the data, so the only thing we are setting is the context, in the sense of which genome to work with.

We recommend using the fruit fly (Drosophila melanogaster), as it is widely used as an experimental model organism and its genome is small enough that you should not run into any problem in terms of scalability. You can access the data by going to the fruit fly page on the NCBI websiteLinks to an external site. and scrolling down to the download button. Select "RefSeq only" and download the genome sequences (FASTA), the annotation (GTF) and the transcripts (FASTA). You will get a zip file of just under 68MB.

If you want to test larger genomes, you could use the mouse (mus musculus) genome. The latest version is GRCm39, which you can access here: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/Links to an external site.. You can follow the same step as above, and get a zip file just over 975MB.

Tasks
The overall objective of this project is to understand how robust and how reproducible genomic analyses are. You will consider the whole pipeline, starting from a set of reads and ending with a list of differentially expressed genes.

To generate sets reads for the genome of your choice, we recommend using an existing tool rather than implementing your own approach. There are several options available, such as polyesterLinks to an external site.. The choice of tool is up to you. Using these tools, you can simulate having multiple experiments for the same conditions. 

You will also consider noise in the data. Here, we are particularly interested in what happens if you shuffle the reads in the sets you have generated. For instance, if a file contains read A, read B, read C and read D, will it lead to the same results as a file that contains read D, read B, read A and read C? You could also look at introducing some sequencing errors in some reads.

For each set of reads that you obtain, you will use two analysis methods. The first one follows the same pipeline as earlier in the unit: alignment to the reference genome (with Bowtie2), quantification at the gene level (with featureCounts), and differential expression analysis (with DESeq2). The second is performing quantification using a pseudoalignment (using your choice of KallistoLinks to an external site. or SalmonLinks to an external site.), followed by differential expression analysis with DESeq2.

Please pay close attention to the user manuals to know which files to use with each tool.

The overall goal is to evaluate the level of consistency in the results. Hints for your analysis are given below.

IFN646_reproducibility_project.drawio.png

How consistent are the results within the same method, for different variations on the input? (green lines)
How consistent are the results for the same input, for different methods? (purple lines)
(note that not all lines are shown)
