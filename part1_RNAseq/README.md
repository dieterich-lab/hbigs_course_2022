# HBIGS course 2021/22: Introduction to computational RNA biology
## Part 1 - From raw sequencing reads to differential expression

## Description

This part of the lecture contains material for the processing and analysis of the *RNA-seq* data from the example dataset.

Participants will gain practical experience and skills to be able to:

* Perform command-line Linux-based analyses using an HPC cluster
* Assess the quality of RNA-seq data
* Pre-process (adapter removal, *etc.* ) and align RNA-seq data to a reference genome
* Perform differential expression analysis
* Understand the different file formats used in a standard RNA-seq workflow (FASTQ, FASTA, BAM, GTF, *etc.* )
* Understand essential aspects of alignment, abundance estimation, experimental design, *etc.*

**Note:** Aspects such as *(i)* library preparation and sequencing, or *(ii)* more advanced topics such as performing a reference guided or a *de novo* assembly of transcripts will NOT be covered in this introductory lecture.

## Dependencies

All the core software that will be used in this part of the lecture is already available in the HPC cluster using 
the [environment modules](http://modules.sourceforge.net/). 
Some scripts require R-Package dependencies, and must be installed first.

```bash
# load the R module
module load R
# start R-Console
R
```

```python
# once in the R-console
install.packages("optparse")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")
```

**Note:** It is possible that you are prompted with the following message, if you install libraries for the first time in your account

```bash
Warning in install.packages("optparse") :
  'lib = "/beegfs/biosw/R/4.0.5_deb10/lib/R/library"' is not writable
Would you like to use a personal library instead? (yes/No/cancel) 
```

You can safely answer `yes` to these questions. The R packages will be installed locally in your $HOME directory.
After installing packages, to quit the R-Console, just type `quit()`. You don't need to save your workspace.


## RNA-seq workflow

RNA sequencing (RNA-seq) is used to interrogate in an unbiased manner the transcriptome, by quantifying RNA transcript abundance and diversity.
In this first part of the course, we will go through a standard RNA-seq workflow using `Make`. The make utility requires a file called `Makefile` (or `makefile`), which defines a set of tasks to be executed. These tasks can be grouped into one or several rules, in the form of

```
target: prerequisites
<TAB> recipe
```

The target, prerequisites, and recipes together make a rule. For this part of the course, you don't need to know more about makefiles! 
This is sufficient to have an idea of what the `Makefile` provided here is doing for each step of the workflow.
Although all rules can be called at once, we will call them one by one, and explain the different steps of the standard RNA-seq workflow.

### Workflow in detail

1. Quality assessment with [**MultiQC**](https://multiqc.info/)
2. Removing sequencing adapters with [**FlexBar**](https://github.com/seqan/flexbar)
3. Removal of ribosomal RNA with  [**bowtie2**](https://github.com/BenLangmead/bowtie2)
4. Read mapping with  [**STAR**](https://github.com/alexdobin/STAR)
5. Read-to-gene assignment (abundance estimation, or read counting) with [**subread/featurecount**](http://subread.sourceforge.net/)


### Data

The data for this part of the course is available under `/pub/hbigs_course_2022/part1_RNAseq`.

