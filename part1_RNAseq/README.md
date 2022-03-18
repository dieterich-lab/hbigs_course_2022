# HBIGS course 2021/22: Introduction to computational RNA biology: Part 1 - From raw sequencing reads to count tables and beyond

Participants will gain practical experience and skills to be able to:

* Perform command-line Linux-based analyses using an HPC cluster
* Understand the different file formats used in a standard RNA-seq workflow (FASTQ, FASTA, BAM, GTF, *etc.* )
* Assess quality of RNA-seq data
* Pre-process (adapter removal, *etc.* ) and align RNA-seq data to a reference genome
* Perform differential expression analysis

**Note:** Aspects such as *(i)* library preparation and sequencing, or *(ii)* more advanced topics such as performing a reference guided or a de novo assembly of transcripts will NOT be covered in this introductory lecture.

## Description

RNA sequencing (RNA-seq) is used to interrogate in an unbiased manner the transcriptome, by quantifying RNA transcript abundance and diversity.
In this first part of the course, we will go through a standard RNA-seq workflow using `Make`. The make utility requires a file called `Makefile` (or `makefile`), which defines a set of tasks to be executed. These tasks can be grouped into one or several rules, in the form of

```
target: prerequisites
<TAB> recipe
```

The target, prerequisites, and recipes together make a rule. For this part of the course, you don't need to know more about makefiles. 
This is sufficient to have an idea of what the `Makefile` provided here is doing for each step of the workflow.
Although all rules can be called at once, we will call them one by one, and explain the different steps of the standard RNA-seq workflow.

### Workflow in detail

1. Quality assessment with [**MultiQC**](https://multiqc.info/)
1. Removing sequencing adapters with [**FlexBar**](https://github.com/seqan/flexbar)
1. Removal of ribosomal RNA with  [**bowtie2**](https://github.com/BenLangmead/bowtie2)
1. Principal read mapping with  [**STAR**](https://github.com/alexdobin/STAR)
2. Read-to-gene assignment with [**subread/featurecount**](http://subread.sourceforge.net/)


### Data

The data for this part of the course is available under `/pub/hbigs_course_2022/part1_RNAseq`.

